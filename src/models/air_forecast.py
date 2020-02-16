import pandas as pd
import numpy as np
import math
import logging
from src.data.config import fountain, surface, site, option, dates, folders

pd.options.mode.chained_assignment = None  # Suppress Setting with warning


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
    tc_i = 1.6 # W/m K
    rho_w = 1000  # Density of water
    rho_i = 916  # Density of Ice rho_i
    rho_a = 1.29  # kg/m3 air density at mean sea level
    k = 0.4  # Van Karman constant
    bc = 5.670367 * math.pow(10, -8)  # Stefan Boltzman constant

    """Miscellaneous"""
    time_steps = 5 * 60  # s Model time steps
    p0 = 1013  # Standard air pressure hPa
    ftl = 0  # Fountain flight time loss ftl
    dx = 0.001  # Ice layer thickness dx

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

    """ Estimating Fountain Spray radius """
    R_f = (
        df["r_f"].replace(0, np.NaN).mean()
    )  # todo implement variable spray radius for variable discharge

    H = fountain["h_f"]

    """ Simulation """
    for i in range(1, df.shape[0]):

        # Ice Melted
        if df.loc[i - 1, "iceV"] <  0:
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
            df["r_ice"] = R_f
            ice_layer = dx * math.pi * math.pow(df.loc[i, "r_ice"], 2) * rho_i
            df.loc[i - 1, "ice"] = ice_layer / rho_i
            df.loc[i - 1, "iceV"] = ice_layer/rho_i
            logger.critical(
                "Ice layer initialised with volume %s thick at %s", df.loc[i - 1, "iceV"], df.loc[ i - 1, "When"]
            )

        if state == 1:

            """ Truncated cone SA """
            # Ice Height
            coefficients = [1, -3 * H + 2 * H**2,  H**2, -3 * df.loc[i - 1, "iceV"] * H**2/ (math.pi * df.loc[i, "r_ice"]**2)]

            p = np.poly1d(coefficients)
            r = np.roots(p)

            for j in range(len(r)):
                if np.isreal(r[j]):
                    df.loc[i, "h_ice"] = np.real(r[j])

            df.loc[i, "r_top"] = R_f/H * ( H - df.loc[i, "h_ice"])

            df.loc[i, "SA"] = (
                    math.pi * df.loc[i, "r_top"]**2
                    + math.pi * df.loc[i, "r_ice"] * math.pow( H**2 + df.loc[i, "r_ice"]**2, 1/2 )
                    - math.pi * df.loc[i, "r_top"] * math.pow( (H-df.loc[i, "h_ice"])** 2 + df.loc[i, "r_top"] ** 2, 1 / 2)
                )

            logger.info(
                "Ice area is %s and height is %s at %s",
                df.loc[i, "SA"],
                df.loc[i, "h_ice"],
                df.loc[i, "When"],
            )


            # update AIR Ice layer
            ice_layer = dx * df.loc[i, "SA"] * rho_i
            logger.debug(
                "Ice layer is %s thick at %s", ice_layer, df.loc[i, "When"]
            )

            # Precipitation to ice quantity
            if df.loc[i, "T_a"] < surface["rain_temp"]:
                df.loc[i, "ppt"] = (
                    surface["snow_fall_density"]
                    * df.loc[i, "Prec"]
                    * math.pi
                    * math.pow(df.loc[i, "r_ice"], 2)
                )
                df.loc[i, "solid"] = df.loc[i, "solid"] + df.loc[i, "ppt"]

            # Fountain water output
            df.loc[i, "liquid"] = df.loc[i, "Discharge"] * (1 - ftl) * time_steps / 60


            # """ When fountain run """
            # if df.loc[i, "liquid"] > 0:
            #
            #     # Conduction Freezing
            #     if df.loc[i - 1, "T_s"] < 0:
            #
            #         df.loc[i, "solid"] += (ice_layer * ci * (-df.loc[i - 1, "T_s"])) / (
            #             Lf
            #         )
            #
            #         if df.loc[i, "solid"] > df.loc[i, "liquid"]:
            #             df.loc[i, "solid"] = df.loc[i, "liquid"]
            #             df.loc[i, "liquid"] = 0
            #         else:
            #             df.loc[i, "liquid"] -= (
            #                 ice_layer * ci * (-df.loc[i - 1, "T_s"])
            #             ) / Lf
            #
            #         logger.debug(
            #             "Conduction made Ice layer %s thick", df.loc[i, "solid"],
            #         )
            #         df.loc[i, "delta_T_s"] = -df.loc[i - 1, "T_s"]

            # # Conduction Heating and Melting
            # df.loc[i, "delta_T_s"] = tc_i * df.loc[i, "SA"] * (df.loc[i - 1, "T_a"] - df.loc[i - 1, "T_s"]) * time_steps/dx

            # if (df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]) > 0:
            #
            #     # Melting Ice by Temperature
            #     df.loc[i, "solid"] -= (
            #         (ice_layer * ci)
            #         * (-(df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
            #         / (-Lf)
            #     )
            #
            #     df.loc[i, "melted"] += (
            #         (ice_layer * ci)
            #         * (-(df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
            #         / (-Lf)
            #     )
            #
            #     df.loc[i - 1, "T_s"] = 0
            #     df.loc[i, "delta_T_s"] = 0
            #
            #     logger.debug(
            #         "Conduction melted ice is %s ", round(df.loc[i, "melted"]),
            #     )

            """ Energy Balance starts """

            df.loc[i - 1, "vp_ice"] = 6.112 * np.exp(
                22.46 * (df.loc[i - 1, "T_s"]) / ((df.loc[i - 1, "T_s"]) + 243.12)
            )

            # Sublimation only
            df.loc[i, "Ql"] = (
                0.623
                * Ls
                * rho_a
                / p0
                * math.pow(k, 2)
                * df.loc[i - 1, "v_a"]
                * (df.loc[i - 1, "vp_a"] - df.loc[i - 1, "vp_ice"])
                / (
                    np.log(surface["h_aws"] / surface["z0mi"])
                    * np.log(surface["h_aws"] / surface["z0hi"])
                )
            )

            if df.loc[i, "Ql"] < 0:
                df.loc[i, "gas"] -= (
                    df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps
                ) / Ls

                # Removing gas quantity generated from previous ice
                df.loc[i, "solid"] += (
                    df.loc[i, "Ql"] * (df.loc[i, "SA"]) * time_steps
                ) / Ls

                # Ice Temperature
                df.loc[i, "delta_T_s"] += (
                    df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps
                ) / (ice_layer * ci)

                logger.debug(
                    "Gas made after sublimation is %s ", df.loc[i, "gas"],
                )
                logger.debug(
                    "Solid left after sublimation is %s ", df.loc[i, "solid"],
                )

            else:  # Deposition

                df.loc[i, "deposition"] += (
                    df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps
                ) / Ls

                # Adding new deposit
                df.loc[i, "solid"] += (
                    df.loc[i, "Ql"] * (df.loc[i, "SA"]) * time_steps
                ) / Ls

                logger.debug(
                    "Ice made after deposition is %s thick",
                    round(df.loc[i, "deposition"]),
                )


            df.loc[i, "SRf"] = (
                ((df.loc[i, "r_ice"] - df.loc[i, "r_top"]) * H + df.loc[i,"r_top"] * df.loc[i,"h_ice"])
                * math.cos(df.loc[i, "SEA"])
                + math.pi
                * math.pow(df.loc[i, "r_top"], 2)
                * math.sin(df.loc[i, "SEA"])
            ) / df.loc[i, "SA"]

            # df.loc[i, "SRf"] = 1

            # Short Wave Radiation SW
            df.loc[i, "SW"] = (1 - df.loc[i, "a"]) * (
                df.loc[i, "Rad"] * df.loc[i, "SRf"] + df.loc[i, "DRad"]
            )

            # Long Wave Radiation LW
            if "oli000z0" not in list(df.columns):

                df.loc[i, "LW"] = df.loc[i, "e_a"] * bc * math.pow(
                    df.loc[i - 1, "T_a"] + 273.15, 4
                ) - surface["ie"] * bc * math.pow(df.loc[i - 1, "T_s"] + 273.15, 4)
            else:
                df.loc[i, "LW"] = df.loc[i, "oli000z0"] - surface["ie"] * bc * math.pow(
                    df.loc[i - 1, "T_s"] + 273.15, 4
                )

            # Sensible Heat Qs
            df.loc[i, "Qs"] = (
                ci
                * rho_a
                * df.loc[i - 1, "p_a"]
                / p0
                * math.pow(k, 2)
                * df.loc[i - 1, "v_a"]
                * (df.loc[i - 1, "T_a"] - df.loc[i - 1, "T_s"])
                / (
                    np.log(surface["h_aws"] / surface["z0mi"])
                    * np.log(surface["h_aws"] / surface["z0hi"])
                )
            )

            # Total Energy W/m2
            df.loc[i, "TotalE"] = df.loc[i, "SW"] + df.loc[i, "LW"] + df.loc[i, "Qs"]

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

                    logger.debug(
                        "Ice made after energy neg is %s thick at temp %s",
                        round(df.loc[i, "solid"]),
                        df.loc[i - 1, "T_s"],
                    )

                else:
                    """ When fountain off and energy negative """
                    # Cooling Ice
                    df.loc[i, "delta_T_s"] += (df.loc[i, "EJoules"]) / (ice_layer * ci)

                    logger.debug(
                        "Ice cooled after energy neg is %s C at %s",
                        round(df.loc[i, "delta_T_s"]),
                        df.loc[i, "When"],
                    )

            else:

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

                    logger.debug(
                        "Ice melted is %s thick", round(df.loc[i, "melted"]),
                    )

            logger.debug(
                "Ice is %s, Solid is %s and meltwater is %s at %s",
                df.loc[i - 1, "ice"],
                df.loc[i, "solid"],
                df.loc[i, "meltwater"],
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
            df.loc[i, "iceV"] = df.loc[i - 1, "ice"] / rho_i

            logger.info(
                "Ice volume is %s and meltwater is %s at %s",
                df.loc[i, "iceV"],
                df.loc[i, "meltwater"],
                df.loc[i, "When"],
            )

    df = df[start:i]

    print("Ice Mass Remaining", float(df["ice"].tail(1)))
    print("Meltwater", float(df["meltwater"].tail(1)))
    print("Ice Volume Max", float(df["iceV"].max()))
    print("Fountain sprayed", float(df["sprayed"].tail(1)))
    print("Ppt", df["ppt"].sum())
    print("Sublimated", float(df["vapour"].tail(1)))
    print("Model ended", df.loc[i - 1, "When"])
    print("Model runtime", df.loc[i - 1, "When"] - df.loc[start, "When"])
    print("Max growth rate", float(df["solid"].max() / 5))
    print(
        "Fountain efficiency",
        float((df["meltwater"].tail(1) + df["ice"].tail(1)) / df["sprayed"].tail(1))
        * 100,
    )

    return df
