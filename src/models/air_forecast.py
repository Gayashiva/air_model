import pandas as pd
import numpy as np
import math
import logging
from src.data.config import fountain, surface, option

pd.options.mode.chained_assignment = None  # Suppress Setting with warning


def albedo(df, surface):

    surface["decay_t_d"] = (
        surface["decay_t_d"] * 24 * 60 / 5
    )  # convert to 5 minute time steps
    surface["decay_t_w"] = (
        surface["decay_t_w"] * 24 * 60 / 5
    )  # convert to 5 minute time steps
    s = 0
    f = 0

    for i in range(1, df.shape[0]):

        if df.loc[i, "T_s"] > -2:  # Wet ice
            ti = surface["decay_t_w"]
            a_min = surface["a_mw"]
        else:
            ti = surface["decay_t_d"]
            a_min = surface["a_md"]

        # Precipitation
        if (df.loc[i, "Fountain"] == 0) & (df.loc[i, "Prec"] > 0):
            if df.loc[i, "T_a"] < surface["rain_temp"]:  # Snow
                s = 0
                f = 0
            else:  # Rainfall
                ti = surface["decay_t_w"]
                a_min = surface["a_mw"]

        if df.loc[i, "Fountain"] > 0:
            f = 1
            s = 0
            ti = surface["decay_t_w"]
            a_min = surface["a_mw"]

        if f == 0:  # last snowed
            df.loc[i, "a"] = a_min + (surface["a_s"] - a_min) * math.exp(-s / ti)
            s = s + 1
        else:  # last sprayed
            df.loc[i, "a"] = a_min + (surface["a_i"] - a_min) * math.exp(-s / ti)
            s = s + 1

    return df["a"]


def projectile_xy(v, theta_f, hs=0.0, g=9.8):
    """
    calculate a list of (x, y) projectile motion data points
    where:
    x axis is distance (or range) in meters
    y axis is height in meters
    v is muzzle velocity of the projectile (meter/second)
    theta_f is the firing angle with repsect to ground (degrees)
    hs is starting height with respect to ground (meters)
    g is the gravitational pull (meters/second_square)
    """
    data_xy = []
    t = 0.0
    theta_f = math.radians(theta_f)
    while True:
        # now calculate the height y
        y = hs + (t * v * math.sin(theta_f)) - (g * t * t) / 2
        # projectile has hit ground level
        if y < 0:
            break
        # calculate the distance x
        x = v * math.cos(theta_f) * t
        # append the (x, y) tuple to the list
        data_xy.append((x, y))
        # use the time in increments of 0.1 seconds
        t += 0.01
    return x, t


def icestupa(df, fountain, surface):

    logger = logging.getLogger(__name__)
    logger.debug("This is a debug message")
    logger.info("This is for temp")
    logger.warning("This is for solid")
    logger.error("This is for melted")
    logger.critical("This is a critical message")

    """Constants"""
    Ls = 2848 * 1000  # J/kg Sublimation
    Le = 2514 * 1000  # J/kg Evaporation
    Lf = 334 * 1000  #  J/kg Fusion
    cw = 4.186 * 1000  # J/kg Specific heat water
    ci = 2.108 * 1000  # J/kgC Specific heat ice
    rho_w = 1000  # Density of water
    rho_i = 916  # Density of Ice rho_i
    rho_a = 1.29  # kg/m3 air density at mean sea level
    k = 0.4  # Van Karman constant
    bc = 5.670367 * math.pow(10, -8)  # Stefan Boltzman constant
    g = 9.8  # gravity

    """Miscellaneous"""
    z = 2  # m height of AWS
    c = 0.5  # Cloudiness c
    time_steps = 5 * 60  # s Model time steps
    p0 = 1013  # Standard air pressure hPa
    theta_f = 45  # Fountain aperture angle
    ftl = 0  # Fountain flight time loss ftl
    dx = 0.001  # Ice layer thickness dx

    """Initialise"""
    fountain_height = fountain["h_f"]
    df["h_f"] = fountain_height
    df["e_s"] = surface["ie"]  # Surface Emissivity
    prec = 0  # Precipitation
    start = 0  # model start step
    state = 0
    ice_layer = 0

    l = [
        "T_s",  # Surface Temperature
        "delta_T_s",  # Temperature Change
        "ice",
        "iceV",
        "solid",
        "liquid",
        "vapour",
        "melted",
        "gas",
        "water",
        "sprayed",
        "TotalE",
        "SW",
        "LW",
        "Qs",
        "Ql",
        "meltwater",
        "SA",
        "h_ice",
        "r_ice",
        "SRf",
        "t_droplet",
        "vpa",
        "vp_ice",
        "vp_w",
    ]
    for col in l:
        df[col] = 0

    """ Estimating Albedo """
    df["a"] = albedo(df, surface)

    """ Estimating Fountain Spray radius """
    Area = math.pi * math.pow(fountain["aperture_f"], 2) / 4

    for j in range(1, df.shape[0]):
        if option != "schwarzsee":
            df.loc[j, "Discharge"] = fountain["discharge"] * df.loc[j, "Fountain"]
        df.loc[j, "v_f"] = df.loc[j, "Discharge"] / (60 * 1000 * Area)
        df.loc[j, "r_f"], df.loc[j, "t_droplet"] = projectile_xy(
            df.loc[j, "v_f"], theta_f, fountain["h_f"]
        )
    R_f = (
        df["r_f"].replace(0, np.NaN).mean()
    )  # todo implement variable spray radius for variable discharge
    D_t = df["t_droplet"].replace(0, np.NaN).mean()

    """ Simulation """
    for i in range(1, df.shape[0]):

        # Ice Melted
        if (df.loc[i - 1, "r_ice"] <= 0) & (df.Discharge[i:].sum() == 0):
            if df.Discharge[i:].sum() == 0:  # If ice melted after fountain run
                df.loc[i - 1, "solid"] = 0
                df.loc[i - 1, "ice"] = 0
                df.loc[i - 1, "iceV"] = 0
                break

            if start != 0:  # If ice melted in between fountain run
                df.loc[i - 1, "solid"] = 0
                df.loc[i - 1, "ice"] = 0
                df.loc[i - 1, "iceV"] = 0
                state = 0

        # Initiate ice formation
        if (df.loc[i, "Discharge"] > 0) & (start == 0):
            state = 1
            start = i - 1  # Set Model start time
            df.loc[i - 1, "r_ice"] = R_f
            df.loc[i - 1, "h_ice"] = 0

        if state == 1:

            """ Keeping r_ice = R_f to determine SA """
            if (df.Discharge[i] > 0) & (df.loc[i - 1, "r_ice"] >= R_f):
                # Ice Radius
                df.loc[i, "r_ice"] = df.loc[
                    i - 1, "r_ice"
                ]  # Ice radius same as Initial Fountain Spray Radius

                # Ice Height
                df.loc[i, "h_ice"] = (
                    3 * df.loc[i - 1, "iceV"] / (math.pi * df.loc[i, "r_ice"] ** 2)
                )

                # Height by Radius ratio
                df.loc[i, "h_r"] = df.loc[i - 1, "h_ice"] / df.loc[i - 1, "r_ice"]

                # Area of Conical Ice Surface
                df.loc[i, "SA"] = (
                    math.pi
                    * df.loc[i, "r_ice"]
                    * math.pow(
                        (
                            math.pow(df.loc[i, "r_ice"], 2)
                            + math.pow((df.loc[i, "h_ice"]), 2)
                        ),
                        1 / 2,
                    )
                )

            else:

                if df.loc[i - 1, "r_ice"] > R_f:
                    logger.warning(
                        "Ice radius increased to %s thick at %s",
                        df.loc[i - 1, "r_ice"],
                        df.loc[i, "When"],
                    )

                """ Keeping h_r constant to determine SA """
                # Height to radius ratio
                df.loc[i, "h_r"] = df.loc[i - 1, "h_ice"] / df.loc[i - 1, "r_ice"]

                # Ice Radius
                df.loc[i, "r_ice"] = math.pow(
                    df.loc[i - 1, "iceV"] / math.pi * (3 / df.loc[i, "h_r"]), 1 / 3
                )

                # Ice Height
                df.loc[i, "h_ice"] = df.loc[i, "h_r"] * df.loc[i, "r_ice"]

                # Area of Conical Ice Surface
                df.loc[i, "SA"] = (
                    math.pi
                    * df.loc[i, "r_ice"]
                    * math.pow(
                        (
                            math.pow(df.loc[i, "r_ice"], 2)
                            + math.pow(df.loc[i, "r_ice"] * df.loc[i, "h_r"], 2)
                        ),
                        1 / 2,
                    )
                )

            logger.info(
                "Ice radius is %s and height is %s at %s",
                df.loc[i, "r_ice"],
                df.loc[i, "h_ice"],
                df.loc[i, "When"],
            )

            # Initialize AIR ice layer and update
            if ice_layer == 0:

                ice_layer = dx * df.loc[i, "SA"] * rho_i
                logger.warning(
                    "Ice layer initialised %s thick at %s", ice_layer, df.loc[i, "When"]
                )
                df.loc[i - 1, "ice"] = ice_layer

            else:
                ice_layer = dx * df.loc[i, "SA"] * rho_i
                logger.info("Ice layer is %s thick at %s", ice_layer, df.loc[i, "When"])

            # Precipitation to ice quantity
            if df.loc[i, "T_a"] < surface["rain_temp"]:
                prec = prec + surface["snow_fall_density"] * df.loc[
                    i, "Prec"
                ] * math.pi * math.pow(df.loc[i, "r_ice"], 2)
                df.loc[i - 1, "ice"] = df.loc[i - 1, "ice"] + surface[
                    "snow_fall_density"
                ] * df.loc[i, "Prec"] * math.pi * math.pow(df.loc[i, "r_ice"], 2)

            # Vapor Pressure empirical relations
            if "vp_a" not in list(df.columns):
                df.loc[i, "vpa"] = (
                    6.11
                    * math.pow(
                        10, 7.5 * df.loc[i - 1, "T_a"] / (df.loc[i - 1, "T_a"] + 237.3)
                    )
                    * df.loc[i, "RH"]
                    / 100
                )
            else:
                df.loc[i, "vpa"] = df.loc[i, "vp_a"]

            df.loc[i, "vp_w"] = 6.112 * np.exp(
                17.62 * surface["T_f"] / (surface["T_f"] + 243.12)
            )
            df.loc[i, "vp_ice"] = 6.112 * np.exp(
                22.46 * (df.loc[i - 1, "T_s"]) / ((df.loc[i - 1, "T_s"]) + 243.12)
            )

            # atmospheric emissivity
            if df.loc[i, "Prec"] > 0:  # c = 1
                df.loc[i, "e_a"] = (
                    1.24
                    * math.pow(
                        abs(df.loc[i, "vpa"] / (df.loc[i, "T_a"] + 273.15)), 1 / 7
                    )
                    * 1.22
                )
            else:
                df.loc[i, "e_a"] = (
                    1.24
                    * math.pow(
                        abs(df.loc[i, "vpa"] / (df.loc[i, "T_a"] + 273.15)), 1 / 7
                    )
                    * (1 + 0.22 * math.pow(c, 2))
                )

            # Fountain water output
            df.loc[i, "liquid"] = df.loc[i, "Discharge"] * (1 - ftl) * time_steps / 60

            """ When fountain run """
            if df.loc[i, "liquid"] > 0:

                # Initial freeze up due to cold ice layer
                if df.loc[i - 1, "T_s"] < 0:

                    df.loc[i, "solid"] += (ice_layer * ci * (-df.loc[i - 1, "T_s"])) / (
                        Lf
                    )

                    if df.loc[i, "solid"] > df.loc[i, "liquid"]:
                        df.loc[i, "solid"] = df.loc[i, "liquid"]
                        df.loc[i, "liquid"] = 0
                    else:
                        df.loc[i, "liquid"] -= (
                            ice_layer * ci * (-df.loc[i - 1, "T_s"])
                        ) / Lf

                    logger.error(
                        "Ice layer made %s thick ice at %s",
                        df.loc[i, "solid"],
                        df.loc[i, "When"],
                    )
                    df.loc[i, "delta_T_s"] = -df.loc[i - 1, "T_s"]

            # Ice surface assumed
            df.loc[i, "e_s"] = surface["ie"]

            # Sublimation
            df.loc[i, "Ql"] = (
                0.623
                * Ls
                * rho_a
                / p0
                * math.pow(k, 2)
                * df.loc[i, "v_a"]
                * (df.loc[i, "vpa"] - df.loc[i, "vp_ice"])
                / (np.log(z / surface["z0mi"]) * np.log(z / surface["z0hi"]))
            )

            df.loc[i, "gas"] -= (df.loc[i, "Ql"] * (df.loc[i, "SA"]) * time_steps) / Ls

            # Removing gas quantity generated from previous time step
            df.loc[i, "solid"] += (
                df.loc[i, "Ql"] * (df.loc[i, "SA"]) * time_steps
            ) / Ls

            # Ice Temperature
            df.loc[i, "delta_T_s"] += (
                df.loc[i, "Ql"] * (df.loc[i, "SA"]) * time_steps
            ) / (ice_layer * ci)

            """Hot Ice"""
            if df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"] > 0:

                # Melting Ice by Temperature
                df.loc[i, "solid"] -= (
                    (ice_layer * ci)
                    * (-(df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
                    / (-Lf)
                )
                df.loc[i, "melted"] += (
                    (ice_layer * ci)
                    * (-(df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
                    / (-Lf)
                )

                df.loc[i, "delta_T_s"] = -df.loc[i - 1, "T_s"]

            logger.info(
                "Ice made after sublimation is %s thick at %s",
                round(df.loc[i, "solid"]),
                df.loc[i, "When"],
            )

            # Estimating Solar Area fraction if theta_s 45
            df.loc[i, "SRf"] = (
                0.5
                * (df.loc[i, "h_r"] + math.pi)
                / (
                    math.pi
                    * math.pow((1 + math.pow(df.loc[i, "h_r"], 2)), 1 / 2)
                    * math.pow(2, 1 / 2)
                )
            )

            # Short Wave Radiation SW
            df.loc[i, "SW"] = (1 - df.loc[i, "a"]) * (
                df.loc[i, "Rad"] * df.loc[i, "SRf"] + df.loc[i, "DRad"]
            )

            # Long Wave Radiation LW
            if "oli000z0" not in list(df.columns):

                df.loc[i, "LW"] = df.loc[i, "e_a"] * bc * math.pow(
                    df.loc[i, "T_a"] + 273.15, 4
                ) - df.loc[i, "e_s"] * bc * math.pow(df.loc[i - 1, "T_s"] + 273.15, 4)
            else:
                df.loc[i, "LW"] = df.loc[i, "oli000z0"] - df.loc[
                    i, "e_s"
                ] * bc * math.pow(df.loc[i - 1, "T_s"] + 273.15, 4)

            # Sensible Heat Qs
            df.loc[i, "Qs"] = (
                ci
                * rho_a
                * df.loc[i, "p_a"]
                / p0
                * math.pow(k, 2)
                * df.loc[i, "v_a"]
                * (df.loc[i, "T_a"] - df.loc[i - 1, "T_s"])
                / (np.log(z / surface["z0mi"]) * np.log(z / surface["z0hi"]))
            )

            # Total Energy W/m2
            df.loc[i, "TotalE"] = df.loc[i, "SW"] + df.loc[i, "LW"] + df.loc[i, "Qs"]

            # Total Energy Joules
            df.loc[i, "EJoules"] = df.loc[i, "TotalE"] * time_steps * df.loc[i, "SA"]

            """ Energy negative"""
            if df.loc[i, "EJoules"] < 0:

                """ And fountain on """
                if df.loc[i - 1, "liquid"] > 0:
                    """Freezing water"""

                    df.loc[i, "liquid"] -= (df.loc[i, "EJoules"]) / (-Lf)

                    if df.loc[i, "liquid"] < 0:
                        df.loc[i, "liquid"] += (df.loc[i, "EJoules"]) / (-Lf)
                        df.loc[i, "solid"] += df.loc[i, "liquid"]
                        df.loc[i, "liquid"] = 0
                    else:
                        df.loc[i, "solid"] += (df.loc[i, "EJoules"]) / (-Lf)

                    logger.warning(
                        "Ice made after energy neg is %s thick at temp %s",
                        round(df.loc[i, "solid"]),
                        df.loc[i - 1, "T_s"],
                    )

                else:
                    """ When fountain off and energy negative """

                    if df.loc[i - 1, "liquid"] < 0:
                        logger.error(
                            "Liquid is %s at %s",
                            round(df.loc[i, "delta_T_s"]),
                            df.loc[i, "When"],
                        )
                        df.loc[i - 1, "liquid"] = 0

                    # Cooling Ice
                    df.loc[i, "delta_T_s"] += (df.loc[i, "EJoules"]) / (ice_layer * ci)

                logger.info(
                    "Ice made after energy neg is %s thick at %s",
                    round(df.loc[i, "solid"]),
                    df.loc[i, "When"],
                )

            else:
                """ Energy Positive """

                # Heating Ice
                df.loc[i, "delta_T_s"] += (df.loc[i, "EJoules"]) / (ice_layer * ci)

                """Hot Ice"""
                if (df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]) > 0:

                    # Melting Ice by Temperature
                    df.loc[i, "solid"] -= (
                        (ice_layer * ci)
                        * (-(df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
                        / (-Lf)
                    )

                    df.loc[i, "melted"] += (
                        (ice_layer * ci)
                        * (-(df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
                        / (-Lf)
                    )

                    df.loc[i - 1, "T_s"] = 0
                    df.loc[i, "delta_T_s"] = 0

                logger.info(
                    "Ice melted is %s thick at %s",
                    round(df.loc[i, "solid"]),
                    df.loc[i, "When"],
                )

            if df.loc[i, "delta_T_s"] < -50:
                logger.error(
                    "Temperature change is %s at %s",
                    round(df.loc[i, "delta_T_s"]),
                    df.loc[i, "When"],
                )

            """ Quantities of all phases """
            df.loc[i, "T_s"] = df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]
            df.loc[i, "meltwater"] = df.loc[i - 1, "meltwater"] + df.loc[i, "melted"]
            df.loc[i, "ice"] = df.loc[i - 1, "ice"] + df.loc[i, "solid"]
            df.loc[i, "vapour"] = df.loc[i - 1, "vapour"] + df.loc[i, "gas"]
            df.loc[i, "sprayed"] = (
                df.loc[i - 1, "sprayed"] + df.loc[i, "Discharge"] * time_steps / 60
            )
            df.loc[i, "water"] = df.loc[i - 1, "water"] + df.loc[i, "liquid"]
            df.loc[i, "iceV"] = df.loc[i, "ice"] / rho_i

            logger.debug(
                "Surface temp. %s, is Ice is %s at %s",
                round(df.loc[i, "T_s"]),
                round(df.loc[i, "ice"]),
                df.loc[i, "When"],
            )

    df = df[start:i]

    print("Fountain sprayed", float(df["sprayed"].tail(1)))
    print("Ice Mass Remaining", float(df["ice"].tail(1)))
    print("Meltwater", float(df["meltwater"].tail(1)))
    print("Evaporated/sublimated", float(df["vapour"].tail(1)))
    print("Model ended", df.loc[i - 1, "When"])
    print("Model runtime", df.loc[i - 1, "When"] - df.loc[start, "When"])
    print(
        "Fountain efficiency",
        float((df["meltwater"].tail(1) + df["ice"].tail(1)) / df["sprayed"].tail(1))
        * 100,
    )
    print("Ice Volume Max", float(df["iceV"].max()))
    print("Ppt", prec)

    return df
