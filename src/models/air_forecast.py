import pandas as pd
import numpy as np
import math
import time
import logging
from src.data.config import fountain, surface, site, option, dates, folders
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pd.options.mode.chained_assignment = None  # Suppress Setting with warning


def conduct(T_initial, Q_dot_in = 0):
    A = 1  # cross sectional area of wall element in m^2
    rho = 916.0  # density of wall material in kg / m^3
    kk = 2.1  # thermal conductivity of wall material in W / (m*K)
    c = 2.097 * 1000  # specific heat capacity in J / (kg*K)
    sigma = 5.6704e-08  # Stefan-Boltzmann constant in W * m^-2 * K^-4
    h = 0.0  # convective heat transfer coefficient in W / (m^2 * K)

    T_initial = T_initial  # initial temperature in Kelvin
    T_inf = 273.15

    L = 0.01  # thickness of the entire wall in meters
    N = 10  # number of discrete wall segments
    ddx = L / N  # length of each wall segment in meters

    total_time = 5 * 60.0  # total duration of simulation in seconds
    nsteps = 5000  # number of timesteps
    dt = total_time / nsteps  # duration of timestep in seconds

    # The size of this nondimensional factor gives a rough idea
    # of the stability of the crude numerical integration we are using.
    # If this is too big there will be problems.
    simfac = (kk * dt) / (c * rho * ddx * ddx)

    # this is the factor by which to multiply Q_dot_in and Q_dot_out
    heatfac = ddx / (kk * A)

    # initialize volume element coordinates and time samples
    x = np.linspace(0, ddx * (N - 1), N)

    timesamps = np.linspace(0, dt * nsteps, nsteps + 1)

    # "In the 2-D case with inputs of length M and N, the outputs are of
    # shape (N, M) for 'xy' indexing and (M, N) for 'ij' indexing."
    X, TIME = np.meshgrid(x, timesamps, indexing="ij")

    # initialize a big 2D array to store temperature values
    T = np.zeros((X.shape))

    # set the initial temperature profile of the wall
    for ctr in range(len(x)):
        T[ctr, 0] = T_initial[ctr] + 273.15

    # T[len(x) - 1, :] = 273.15

    for j in range(len(timesamps) - 1):
        # get the outside wall temperature and heat flow at current time
        T_out = T[len(x) - 1, j]
        Q_dot_out = sigma * A * (pow(T_out, 4) - pow(T_inf, 4)) + h * A * (T_out - T_inf)

        # now compute temperature at the outside boundary for the next time step
        T[len(x)-1, j+1] = T_out + simfac * (T[len(x)-2, j] - T_out - heatfac * Q_dot_out)

        # and now compute temperature at the inside boundary for the next time step
        T[0, j + 1] = T[0, j] + simfac * (T[1, j] - T[0, j] + heatfac * Q_dot_in)
        # now loop through the interior elements to get their temp for the next time
        for ctr in range(len(x) - 2):
            T[ctr + 1, j + 1] = T[ctr + 1, j] + simfac * (
                T[ctr, j] - 2 * T[ctr + 1, j] + T[ctr + 2, j]
            )

    return T[:,-1] - 273.15


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
    cw = 4.186 * 1000  # J/kgC Specific heat water
    ci = 2.097 * 1000  # J/kgC Specific heat ice
    rho_w = 1000  # Density of water
    rho_i = 916  # Density of Ice rho_i
    rho_a = 1.29  # kg/m3 air density at mean sea level
    k = 0.4  # Van Karman constant
    bc = 5.670367 * math.pow(10, -8)  # Stefan Boltzman constant

    """Miscellaneous"""
    time_steps = 10 * 60  # s Model time steps
    p0 = 1013  # Standard air pressure hPa
    ftl = 0  # Fountain flight time loss ftl
    dx = 0.001  # Ice layer thickness dx #todo temperature gradient required to model

    """Initialise"""
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
        "Qc",
        "meltwater",
        "SA",
        "h_ice",
        "r_ice",
        "SRf",
        "vp_ice",
        "ppt",
        "deposition",
    ]
    for col in l:
        df[col] = 0

    T_layers = np.zeros(10)

    """ Estimating Fountain Spray radius """
    R_f = (
        df["r_f"].replace(0, np.NaN).mean()
    )  # todo implement variable spray radius for variable discharge

    # df = df.set_index("When").resample("15M").mean().reset_index()
    # print(df)

    """ Simulation """
    for i in tqdm(range(1, df.shape[0])):

        # Ice Melted
        if df.loc[i - 1, "iceV"] <= 0:
            df.loc[i - 1, "solid"] = 0
            df.loc[i - 1, "ice"] = 0
            df.loc[i - 1, "iceV"] = 0
            if df.Discharge[i:].sum() == 0:  # If ice melted after fountain run
                break
            else:  # If ice melted in between fountain run
                state = 0

        # Initiate ice formation
        if (df.loc[i, "Discharge"] > 0) & (state == 0):
            state = 1
            start = i - 1  # Set Model start time
            df.loc[i - 1, "r_ice"] = R_f
            df.loc[i - 1, "h_ice"] = dx
            ice_layer = dx * math.pi * R_f ** 2 * rho_i
            df.loc[i - 1, "iceV"] = dx * math.pi * R_f ** 2

            logger.debug(
                "Ice layer initialised %s thick at %s", ice_layer, df.loc[i, "When"]
            )

        if state == 1:

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

                # Height to radius ratio
                df.loc[i, "h_r"] = df.loc[i - 1, "h_r"]

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

            df.loc[i, "SRf"] = (
                0.5
                * df.loc[i, "h_ice"]
                * df.loc[i, "r_ice"]
                * math.cos(df.loc[i, "SEA"])
                + math.pi
                * math.pow(df.loc[i, "r_ice"], 2)
                * 0.5
                * math.sin(df.loc[i, "SEA"])
            ) / (
                math.pi
                * math.pow(
                    (math.pow(df.loc[i, "h_ice"], 2) + math.pow(df.loc[i, "r_ice"], 2)),
                    1 / 2,
                )
                * df.loc[i, "r_ice"]
            )

            logger.debug(
                "Ice radius is %s and ice is %s at %s",
                df.loc[i, "r_ice"],
                df.loc[i, "h_ice"],
                df.loc[i, "When"],
            )

            # Update AIR ice layer
            ice_layer = dx * df.loc[i, "SA"] * rho_i
            logger.debug("Ice layer is %s thick at %s", ice_layer, df.loc[i, "When"])

            # Precipitation to ice quantity
            if df.loc[i, "T_a"] < surface["rain_temp"]:
                df.loc[i, "ppt"] = (
                    surface["snow_fall_density"]
                    * df.loc[i, "Prec"]
                    * math.pi
                    * math.pow(df.loc[i, "r_ice"], 2)
                )

            # Fountain water output
            df.loc[i, "liquid"] = df.loc[i, "Discharge"] * (1 - ftl) * time_steps / 60

            """ Energy Balance starts """

            df.loc[i, "vp_ice"] = 6.112 * np.exp(
                22.46 * (df.loc[i - 1, "T_s"]) / ((df.loc[i - 1, "T_s"]) + 272.62)
            )

            df.loc[i, "vp_w"] = 6.112

            # Water Boundary
            if df.Discharge[i] > 0:
                df.loc[i, "vp_s"] = df.loc[i, "vp_w"]
                L = Le
                c_s = cw

            else:
                df.loc[i, "vp_s"] = df.loc[i, "vp_ice"]
                L = Ls
                c_s = ci

            df.loc[i, "Ql"] = (
                0.623
                * L
                * rho_a
                / p0
                * math.pow(k, 2)
                * df.loc[i, "v_a"]
                * (df.loc[i, "vp_a"] - df.loc[i, "vp_s"])
                / (
                    np.log(surface["h_aws"] / surface["z0mi"])
                    * np.log(surface["h_aws"] / surface["z0hi"])
                )
            )

            if df.loc[i, "Ql"] < 0:
                df.loc[i, "gas"] -= (df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps) / L

                # Removing gas quantity generated from previous ice
                df.loc[i, "solid"] += (
                    df.loc[i, "Ql"] * (df.loc[i, "SA"]) * time_steps
                ) / L

                # Ice Temperature

                T_layers = conduct(T_layers, df.loc[i, "Ql"])


                logger.info(
                    "Temperature gradient after Latent %s is %s at %s",
                    df.loc[i, "Ql"],
                    T_layers,
                    df.loc[i, "When"],
                )


                logger.debug(
                    "Gas made after sublimation is %s ", round(df.loc[i, "gas"])
                )

            else:  # Deposition

                df.loc[i, "deposition"] += (
                    df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps
                ) / L

                logger.debug(
                    "Ice made after deposition is %s thick",
                    round(df.loc[i, "deposition"]),
                )

            # Sensible Heat Qs
            df.loc[i, "Qs"] = (
                c_s
                * rho_a
                * df.loc[i, "p_a"]
                / p0
                * math.pow(k, 2)
                * df.loc[i, "v_a"]
                * (df.loc[i, "T_a"] - df.loc[i - 1, "T_s"])
                / (
                    np.log(surface["h_aws"] / surface["z0mi"])
                    * np.log(surface["h_aws"] / surface["z0hi"])
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
                ) - surface["ie"] * bc * math.pow(df.loc[i - 1, "T_s"] + 273.15, 4)
            else:
                df.loc[i, "LW"] = df.loc[i, "oli000z0"] - surface["ie"] * bc * math.pow(
                    df.loc[i - 1, "T_s"] + 273.15, 4
                )

            # Conduction Freezing
            if (df.loc[i, "liquid"] > 0) & (df.loc[i - 1, "T_s"] < 0):
                df.loc[i, "Qc"] = (
                    ice_layer
                    * ci
                    * (-df.loc[i - 1, "T_s"])
                    / (df.loc[i, "SA"] * time_steps)
                )

                T_layers[0] = 0

                logger.debug(
                    "Ice layer made %s thick ice at %s",
                    df.loc[i, "solid"],
                    df.loc[i, "When"],
                )

            # Total Energy W/m2
            df.loc[i, "TotalE"] = (
                df.loc[i, "SW"] + df.loc[i, "LW"] + df.loc[i, "Qs"] + df.loc[i, "Qc"]
            )

            logger.debug(
                "Energy is %s thick at %s", round(df.loc[i, "Qc"]), df.loc[i, "When"]
            )

            # Total Energy Joules
            df.loc[i, "EJoules"] = df.loc[i, "TotalE"] * time_steps * df.loc[i, "SA"]

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

                else:
                    """ When fountain off and energy negative """
                    # Cooling Ice
                    T_layers = conduct(T_layers, df.loc[i, "TotalE"])

                    logger.info(
                        "Temperature gradient after Energy %s is %s at %s",
                        df.loc[i, "TotalE"],
                        T_layers,
                        df.loc[i, "When"],
                    )

                logger.debug(
                    "Ice made after energy neg is %s thick at temp %s",
                    round(df.loc[i, "solid"]),
                    df.loc[i - 1, "T_s"],
                )

            else:

                """Hot Ice"""
                if T_layers[0] > 0:

                    # Melting Ice by Temperature
                    df.loc[i, "solid"] -= (
                        (df.loc[i, "EJoules"])
                        / (Lf*rho_i*df.loc[i, "SA"])
                    )

                    df.loc[i, "melted"] += (
                        (df.loc[i, "EJoules"])
                        / (Lf*rho_i*df.loc[i, "SA"])
                    )

                    T_layers = np.zeros(10)

                else:

                    T_layers = conduct(T_layers, df.loc[i, "TotalE"])

                logger.info(
                    "Ice melted is %s thick at %s",
                    round(df.loc[i, "melted"]),
                    df.loc[i, "When"],
                )

            """ Quantities of all phases """

            logger.debug(
                "Temperature gradient is %s at %s",
                df["T_s"][i-9:i+1].values,
                df.loc[i, "When"],
            )

            for ctr in range(0,10):
                df.loc[i-ctr, "T_s"] = T_layers[ctr]

            logger.info(
                "Temperature gradient is %s at %s",
                T_layers,
                df.loc[i, "When"],
            )

            df.loc[i, "meltwater"] = df.loc[i - 1, "meltwater"] + df.loc[i, "melted"]
            df.loc[i, "ice"] = (
                df.loc[i - 1, "ice"]
                + df.loc[i, "solid"]
                + df.loc[i, "ppt"]
                + df.loc[i, "deposition"]
            )
            df.loc[i, "vapour"] = df.loc[i - 1, "vapour"] + df.loc[i, "gas"]
            df.loc[i, "sprayed"] = (
                df.loc[i - 1, "sprayed"] + df.loc[i, "Discharge"] * time_steps / 60
            )
            df.loc[i, "water"] = df.loc[i - 1, "water"] + df.loc[i, "liquid"]
            df.loc[i, "iceV"] = (df.loc[i, "ice"] - df.loc[i, "ppt"]) / rho_i + df.loc[
                i, "ppt"
            ] / surface["snow_fall_density"]

            logger.info(
                "Ice volume is %s and temperature is %s at %s",
                df.loc[i, "ice"],
                df.loc[i, "T_s"],
                df.loc[i, "When"],
            )

    # df = df[start:i]

    print("Ice Volume Max", float(df["iceV"].max()))
    print(
        "Fountain efficiency",
        float(
            (df["meltwater"].tail(1) + df["ice"].tail(1))
            / (df["sprayed"].tail(1) + df["ppt"].sum() + df["deposition"].sum())
            * 100
        ),
    )

    print("Ice Mass Remaining", float(df["ice"].tail(1)))
    print("Meltwater", float(df["meltwater"].tail(1)))
    print("Sublimated", float(df["vapour"].tail(1)))
    print("Fountain sprayed", float(df["sprayed"].tail(1)))
    print("Ppt", df["ppt"].sum())
    print("Deposition", df["deposition"].sum())
    print("Water not frozen", float(df["water"].tail(1)))
    print("Model ended", df.loc[i - 1, "When"])
    print("Model runtime", df.loc[i - 1, "When"] - df.loc[start, "When"])
    print("Max growth rate", float(df["solid"].max() / 5))

    return df
