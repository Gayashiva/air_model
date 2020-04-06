import pandas as pd
import numpy as np
import math
import time
import logging
from src.data.config import fountain, surface, site, option, dates, folders
from tqdm import tqdm
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
    # T_inf = 273.15

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

        # Q_dot_out = sigma * A * (pow(T_out, 4) - pow(T_inf, 4)) + h * A * (T_out - T_inf)
        Q_dot_out = 0

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

def stefan(Q_dot_in = 0):

    # -----PHYSICAL PARAMETERS------#
    l = 334e3  # specific latent heat [J/kg]
    rho = 1e-6  # density [kg/mm^3]
    K = 0.6e-3  # thermal conductivity [W/(mm*K)]
    c = 4.2e3  # specific heat [J/(kg*K)]
    alpha = K / (rho * c)  # thermal diffusivity [mm^/s]
    beta = c / l  # [1/K]
    L = 20  # [mm] -Length of domain
    t_max = 5 * 60  # [s] -maximum time

    # -----NUMERICAL VALUES------ #
    n = 1e2  # number of iterations
    dx = 0.1  # Steplength in x
    dt = dx ** 2 / (2 * alpha)  # Steplength in t
    ds = dx / (n * beta)  # Increment size for moving boundary when absorbing one walker.
    N_x = math.ceil(L / dx)  # Number of points in spatial domain.
    N_t = math.ceil(t_max / dt)  # Number of points in time domain.
    T = np.zeros((N_x, N_t))  # Matrix representing T(x,t)
    s_vector = np.zeros((N_t))  # Vector representing moving boundary s(t)
    t_j = 0
    s_i = 1
    s = dx

    Q = Q_dot_in * 1e-6  # W/mm2
    T_0 = Q * dx / K  # [degree C] -Temperature at fixed boundary


    while t_j < N_t - 1 and s_i < N_x - 1:

        # -----Boundary condition for fixed boundary-----#
        T[0, t_j] = T_0 * n
        # T(1,t_j)=(exp(t_j*dt)-1)*n; #Time-dependent exponential BC
        # T(1,t_j)=(sin(t_j*dt))*n; #Time-dependent oscillating BC
        s_vector[t_j] = s

        for x_i in np.arange(0, N_x):
            if T[x_i, t_j] < 0:
                sign = - 1
            else:
                sign = 1

            for k in np.arange(0, sign * T[x_i, t_j]):
                p = (- 1) ** round(np.random.uniform(0, 1))

                if x_i + p > 0 and x_i + p <= s_i and x_i <= N_x - 1:

                    T[x_i + p, t_j + 1] = T[x_i + p, t_j + 1] + sign

                else:

                    if x_i + p == s_i + 1:
                        s = s + ds * sign

                        s_i = math.floor(s / dx)

        t_j = t_j + 1

    T = T / n

    s_t_end = np.nonzero(s_vector)[-1][-1]

    s_t_vector = s_vector[0: s_t_end]

    return s_t_vector[-1]/1000


def icestupa(df, fountain, surface):

    logger = logging.getLogger(__name__)
    logger.debug("This is a debug message")
    logger.info("This is for temp")
    logger.warning("This is for solid")
    logger.error("This is for melted")
    logger.critical("This is a critical message")

    """Constants"""
    L_s = 2848 * 1000  # J/kg Sublimation
    L_e = 2514 * 1000  # J/kg Evaporation
    L_f = 334 * 1000  #  J/kg Fusion
    c_w = 4.186 * 1000  # J/kgC Specific heat water
    c_i = 2.108 * 1000  # J/kgC Specific heat ice
    rho_w = 1000  # Density of water
    rho_i = 916  # Density of Ice rho_i
    rho_a = 1.29  # kg/m3 air density at mean sea level
    k = 0.4  # Van Karman constant
    bc = 5.670367 * math.pow(10, -8)  # Stefan Boltzman constant

    """Miscellaneous"""
    time_steps = 5 * 60  # s Model time steps
    p0 = 1013  # Standard air pressure hPa
    ftl = 0  # Fountain flight time loss ftl

    """Initialise"""
    start = 0  # model start step
    state = 0
    ice_thickness = 0

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

    """ Estimating Fountain Spray radius """
    R_f = (
        df["r_f"].replace(0, np.NaN).mean()
    )  # todo implement variable spray radius for variable discharge

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
            df.loc[i - 1, "h_ice"] = surface['dx']
            df.loc[i - 1, "iceV"] = surface['dx'] * math.pi * R_f**2

            logger.debug(
                "Ice layer initialised %s thick at %s", df.loc[i - 1, "iceV"], df.loc[i, "When"]
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
            ice_thickness = surface['dx'] * rho_i
            logger.debug("Ice layer is %s thick at %s", ice_thickness, df.loc[i, "When"])

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
                L = L_e
                c_s = c_w

            else:
                df.loc[i, "vp_s"] = df.loc[i, "vp_ice"]
                L = L_s
                c_s = c_i

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



            if df.loc[i, "Ql"] < 0 :
                df.loc[i, "gas"] -= (df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps) / L

                # Removing gas quantity generated from previous ice
                df.loc[i, "solid"] += (
                    df.loc[i, "Ql"] * (df.loc[i, "SA"]) * time_steps
                ) / L

                # Ice Temperature
                df.loc[i, "delta_T_s"] += (
                    df.loc[i, "Ql"] * time_steps
                ) / (ice_thickness * c_i)

                logger.debug(
                    "Gas made after sublimation is %s ",
                    round(df.loc[i, "gas"]),
                )

            else: # Deposition

                df.loc[i, "deposition"] += (df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps) / L

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
                df.loc[i, "LW"] = df.loc[i, "oli000z0"] - surface["ie"] * bc * math.pow(df.loc[i - 1, "T_s"] + 273.15, 4)

            # Conduction Freezing
            if (df.loc[i, "liquid"] > 0) & (df.loc[i - 1, "T_s"] < 0):
                df.loc[i, "Qc"] = ice_thickness * c_i * (-df.loc[i - 1, "T_s"]) / time_steps
                df.loc[i, "delta_T_s"] = -df.loc[i - 1, "T_s"]

                logger.debug(
                    "Ice layer made %s thick ice at %s",
                    df.loc[i, "solid"],
                    df.loc[i, "When"],
                )

            # Total Energy W/m2
            df.loc[i, "TotalE"] = df.loc[i, "SW"] + df.loc[i, "LW"] + df.loc[i, "Qs"] + df.loc[i, "Qc"]

            logger.debug(
                "Energy is %s thick at %s",
                round(df.loc[i, "Qc"]),
                df.loc[i, "When"],
            )

            # Total Energy Joules
            df.loc[i, "EJoules"] = df.loc[i, "TotalE"] * time_steps * df.loc[i, "SA"]

            if df.loc[i, "EJoules"] < 0:

                """ And fountain on """
                if df.loc[i - 1, "liquid"] > 0:

                    """Freezing water"""

                    df.loc[i, "liquid"] -= (df.loc[i, "EJoules"]) / (-L_f)

                    if df.loc[i, "liquid"] < 0:
                        df.loc[i, "liquid"] += (df.loc[i, "EJoules"]) / (-L_f)
                        df.loc[i, "solid"] += df.loc[i, "liquid"]
                        df.loc[i, "liquid"] = 0
                    else:
                        df.loc[i, "solid"] += (df.loc[i, "EJoules"]) / (-L_f)

                else:
                    """ When fountain off and energy negative """
                    # Cooling Ice
                    df.loc[i, "delta_T_s"] += (df.loc[i, "TotalE"] * time_steps) / (ice_thickness * c_i)

                logger.debug(
                    "Ice made after energy neg is %s thick at temp %s",
                    round(df.loc[i, "solid"]),
                    df.loc[i - 1, "T_s"],
                )

            else:

                # Heating Ice
                df.loc[i, "delta_T_s"] += (df.loc[i, "TotalE"] * time_steps ) / (ice_thickness * c_i)

                """Hot Ice"""
                if (df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]) > 0:

                    Q = (ice_thickness * c_i )* ((df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))/ time_steps

                    melted = stefan(Q)

                    # Melting Ice by Temperature
                    df.loc[i, "solid"] -= (
                        (rho_i * melted * df.loc[i, "SA"])
                    )

                    df.loc[i, "melted"] += (
                        (rho_i * melted * df.loc[i, "SA"])
                    )

                    new = (
                        ( melted * 1000 )
                    )

                    old = (
                            (surface['dx'] * c_i)
                            * ((df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
                            / (L_f)
                            * 1000
                    )

                    print(old, new)

                    # # Melting Ice by Temperature
                    # df.loc[i, "solid"] -= (
                    #         (ice_thickness * c_i * df.loc[i, "SA"])
                    #         * (-(df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
                    #         / (-L_f)
                    # )
                    #
                    # df.loc[i, "melted"] += (
                    #         (ice_thickness * c_i * df.loc[i, "SA"])
                    #         * (-(df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]))
                    #         / (-L_f)
                    # )

                    df.loc[i - 1, "T_s"] = 0
                    df.loc[i, "delta_T_s"] = 0

                logger.debug(
                    "Ice melted because %s is %s thick at %s",
                    Q,
                    df.loc[i, "solid"],
                    df.loc[i, "When"],
                )

            """ Quantities of all phases """
            df.loc[i, "T_s"] = df.loc[i - 1, "T_s"] + df.loc[i, "delta_T_s"]
            df.loc[i, "meltwater"] = df.loc[i - 1, "meltwater"] + df.loc[i, "melted"]
            df.loc[i, "ice"] = df.loc[i - 1, "ice"] + df.loc[i, "solid"] + df.loc[i, "ppt"] + df.loc[i, "deposition"]
            df.loc[i, "vapour"] = df.loc[i - 1, "vapour"] + df.loc[i, "gas"]
            df.loc[i, "sprayed"] = (
                df.loc[i - 1, "sprayed"] + df.loc[i, "Discharge"] * time_steps / 60
            )
            df.loc[i, "water"] = df.loc[i - 1, "water"] + df.loc[i, "liquid"]
            df.loc[i, "iceV"] = (df.loc[i, "ice"] - df.loc[i, "ppt"]) / rho_i + df.loc[i, "ppt"] / surface["snow_fall_density"]

            logger.debug(
                "Ice volume is %s and meltwater is %s at %s",
                df.loc[i, "ice"],
                df.loc[i, "meltwater"],
                df.loc[i, "When"],
            )

    print("Ice Volume Max", float(df["iceV"].max()))
    print(
        "Fountain efficiency",
        float((df["meltwater"].tail(1) + df["ice"].tail(1)) / (
                    df["sprayed"].tail(1) + df["ppt"].sum() + df["deposition"].sum()) * 100
              ))

    print("\n Ice Mass Remaining", float(df["ice"].tail(1)))
    print("Meltwater", float(df["meltwater"].tail(1)))
    print("Fountain sprayed", float(df["sprayed"].tail(1)))
    print("Ppt", df["ppt"].sum())
    print("Model runtime", df.loc[i - 1, "When"] - df.loc[start, "When"])
    print("Max growth rate", float(df["solid"].max()/5))


    return df
