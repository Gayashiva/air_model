import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LightSource
import math
import sys

np.seterr(all="raise")


def euler(T, L=1, d=916, T_S=253.15, dx = 0.0001):


    """Spatial Grid"""
    nx = int(L / dx) + 1
    x = np.linspace(0, L, nx)

    """Time Grid"""
    t_sim = 5 * 60  # 5 minute time step
    dt = 1
    nt = int(t_sim / dt)

    """Material Properties"""
    ci = 2.108 * 1000  # J/kg/K Specific heat ice
    d = 916
    lamda = 2.09  # W/m/K Thermal Conductivity of Ice at 0 C

    """Simulation"""
    for n in range(1, nt):
        Tn = T.copy()

        T[1:-1] = (
            Tn[1:-1]
            + dt * lamda / (d * ci) * (Tn[2:] - 2 * Tn[1:-1] + Tn[0:-2]) / dx ** 2
        )

        T[0] = T_S
        T[-1] = T[-2]

    return x, T


def albedo(df, a_i, a_s, a_w, t):

    t = t * 24 * 60 / 5  # convert to 5 minute time steps
    ti = t
    tw = t
    s = 0  # Initialised
    w = 0
    j = 0  # Account for decay rate after rain
    rf = 2  # Rain decay factor
    Ts = -1  # Solid Ppt

    for i in range(1, df.shape[0]):

        # Wind Sensor Snow covered and fountain off
        if (df.loc[i, "v_a"] == 0) & (df.loc[i, "Discharge"] == 0):
            if df.loc[i, "Prec"] > 0 & (
                df.loc[i, "T_a"] < Ts
            ):  # Assumes snow ppt because of wind sensor
                s = 0
                w = 0
                j = 0
                tw = t  # Decay rate reset after snowfall
            df.loc[i, "a"] = a_i + (a_s - a_i) * math.exp(-s / ti)
            s = s + 1

        # No snow and fountain off
        if (df.loc[i, "v_a"] != 0) & (df.loc[i, "Discharge"] == 0):
            df.loc[i, "a"] = a_w + (a_i - a_w) * math.exp(-w / tw)
            w = w + 1

        # Fountain on
        if df.loc[i, "Discharge"] > 0:
            df.loc[i, "a"] = a_w
            w = 0
            s = 0
            j = 0
            tw = t  # Decay rate reset after fountain

        # Liquid Ppt
        if (df.loc[i, "Prec"] > 0) & (df.loc[i, "T_a"] > Ts):
            if j == 0:
                tw = tw / rf  # Decay rate speeds up after rain
                j = 1
            df.loc[i, "a"] = a_w

    return df["a"]


def projectile_xy(v, theta_f, hs=0.0, g=9.8):
    """
    calculate a list of (x, y) projectile motion data points
    where:
    x axis is distance (or range) in meters
    y axis is height in meters
    v is muzzle velocity of the projectile (meter/second)
    theta_f is the firing angle with repsect to ground (radians)
    hs is starting height with respect to ground (meters)
    g is the gravitational pull (meters/second_square)
    """
    data_xy = []
    t = 0.0
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
    return max(data_xy)[0], t


def icestupa(
    df,
    T_nucleation=-120,
    T_f=5,
    ftl=0.5,
    ie=0.96,
    we=0.95,
    rho_i=916,
    c=0.5,
    a_i=0.6,
    a_s=0.75,
    a_w=0.1,
    t_d=21.9,
    dp=70,
    z0mi=0.001,
    z0ms=0.0015,
    z0hi=0.0001,
    s=0,
    h_f=3,
):

    """Fountain Charecteristics"""
    d_f = 0.005  # Fountain hole diameter
    theta_f = 45  # Fountain aperture angle
    # h_f = 1.35  # Schwarzsee Fountain height
    # ftl = 0.5 # Flight time loss from 0 to 0.9
    # T_f = 5 # c Fountain water temperature from 4 to 9

    """Settings"""
    time_steps = 5 * 60  # s # Model time steps
    z = 2  # m height of AWS
    theta_s = 45  # Solar Angle
    dx = 0.001  # Layers for temp profile

    """Material Properties"""
    Ls = 2848 * 1000  # J/kg Sublimation
    Le = 2514 * 1000  # J/kg Evaporation
    Lf = 334 * 1000  #  J/kg Fusion
    cw = 4.186 * 1000  # J/kg Specific heat water
    ci = 2.108 * 1000  # J/kgC Specific heat ice
    rho_w = 1000  # Density of water
    rho_a = 1.29  # kg/m3 air density at mean sea level
    p0 = 1013  # hPa
    k = 0.4  # Van Karman constant
    bc = 5.670367 * math.pow(10, -8)  # Stefan Boltzman constant

    """Initialise"""
    df["T_s"] = 0  # Surface Temperature
    df["e_s"] = ie  # Surface Emissivity
    df["ice"] = 0
    df["iceV"] = 0
    df["solid"] = 0
    df["liquid"] = 0
    df["vapour"] = 0
    df["melted"] = 0
    df["gas"] = 0
    df["water"] = 0
    df["TotalE"] = 0
    df["SW"] = 0
    df["LW"] = 0
    df["Qs"] = 0
    df["Ql"] = 0
    df["meltwater"] = 0
    df["SA"] = 0  # Surface Area
    df["h_ice"] = 0
    df["h_f"] = h_f
    y = h_f
    df["d_t"] = 0
    prec = 0  # Precipitation
    start = 0  # model start step
    stop = 0  # model end step
    state = 0
    df["SRf"] = 0

    theta_f = math.radians(theta_f)  # Angle of Spray
    theta_s = math.radians(theta_s)  # solar angle

    ''' Estimating Albedo '''
    df["a"] = albedo(df, a_i, a_s, a_w, t_d)

    ''' Estimating Fountain Spray radius '''
    Area = 3.14 * math.pow(d_f, 2) / 4
    df["v"] = 0
    df["r"] = 0

    for j in range(1, df.shape[0]):
        df.loc[j, "v_f"] = (
            df.loc[j, "Discharge"] / (60 * 1000 * Area)
        )
        df.loc[j, "r_f"], df.loc[j, "d_t"] = projectile_xy(
            df.loc[j, "v_f"], theta_f, h_f
        )
    R_f = df["r_f"].replace(0, np.NaN).mean()
    D_t = df["d_t"].replace(0, np.NaN).mean()

    """ Simulation """
    for i in range(1, df.shape[0]):

        if df.loc[i - 1, "iceV"] < 0.001:
            df.loc[i - 1, "solid"] = 0
            df.loc[i - 1, "ice"] = 0
            df.loc[i - 1, "iceV"] = 0
            stop = i - 1

        if (df.loc[i, "When"] > datetime(2019, 2, 14, 18)) & (stop == i - 1):
            break

        if (df.loc[i, "Discharge"]) > 0:
            state = 1
            # Set Model start time
            if df.loc[i - 1, "ice"] > 0:

                if start == 0:
                    start = i - 1

        if state == 1:

            if (
                df.Discharge[i:].sum() > 0
            ):  # Keeping Radius for the given fountain parameters

                # Ice radius same as Fountain Spray Radius
                df.loc[i, "r_ice"] = R_f

                # Ice Height
                df.loc[i, "h_ice"] = (
                    3 * df.loc[i - 1, "iceV"] / (math.pi * df.loc[i, "r_ice"] ** 2)
                )

                df.loc[i, "h_r"] = (
                    df.loc[i, "h_ice"] / df.loc[i, "r_ice"]
                )  # Height to radius ratio

                df.loc[i, "SRf"] = (
                    0.5
                    * (df.loc[i, "h_r"] + math.pi)
                    / (
                        math.pi
                        * math.pow((1 + math.pow(df.loc[i, "h_r"], 2)), 1 / 2)
                        * math.pow(2, 1 / 2)
                    )
                )  # Estimating Solar Area fraction if theta_s 45

                if df.loc[i, "h_f"] > df.loc[i, "h_ice"]:  # Fountain height Check

                    df.loc[i, "SA"] = (
                        math.pi
                        * df.loc[i, "r_ice"]
                        * math.pow(
                            (
                                math.pow(df.loc[i, "r_ice"], 2)
                                + math.pow(df.loc[i, "h_ice"], 2)
                            ),
                            1 / 2,
                        )
                    )  # Area of Conical Ice Surface

                else:

                    y = y + h_f
                    df.loc[i:, "h_f"] = y  # Reset height

                    for j in range(1, df.shape[0]):
                        df.loc[j, "v_f"] = (
                            df.loc[j, "Discharge"]
                            / (60 * 1000 * Area)
                        )
                        df.loc[j, "r_f"], df.loc[j, "d_t"] = projectile_xy(
                            df.loc[j, "v_f"], theta_f, df.loc[i, "h_f"]
                        )
                    R_f = df["r_f"].replace(0, np.NaN).mean()
                    D_t = df["d_t"].replace(0, np.NaN).mean()

                    df.loc[i, "r_ice"] = R_f

                    df.loc[i, "SA"] = (
                        math.pi
                        * df.loc[i, "r_ice"]
                        * math.pow(
                            (
                                math.pow(df.loc[i, "r_ice"], 2)
                                + math.pow(df.loc[i, "h_ice"], 2)
                            ),
                            1 / 2,
                        )
                    )  # Area of Conical Ice Surface

            else:  # Keeping h_r constant to determine SA

                df.loc[i, "h_r"] = (
                    df.loc[i - 1, "h_ice"] / df.loc[i - 1, "r_ice"]
                )  # Height to radius ratio

                # Ice Radius
                df.loc[i, "r_ice"] = math.pow(
                    df.loc[i - 1, "iceV"] / math.pi * (3 / df.loc[i, "h_r"]), 1 / 3
                )  # Assumption height h_r and conical

                df.loc[i, "h_ice"] = df.loc[i, "h_r"] * df.loc[i, "r_ice"]

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
                )  # Area of Conical Ice Surface

            # Precipitation to ice quantity
            prec = prec + dp * df.loc[i, "Prec"] * math.pi * math.pow(R_f, 2)
            df.loc[i - 1, "ice"] = df.loc[i - 1, "ice"] + dp * df.loc[
                i, "Prec"
            ] * math.pi * math.pow(R_f, 2)

            # Vapor Pressure empirical relations
            if "vp_a" not in list(df.columns):
                Ea = (
                    6.11
                    * math.pow(
                        10, 7.5 * df.loc[i - 1, "T_a"] / (df.loc[i - 1, "T_a"] + 237.3)
                    )
                    * df.loc[i, "RH"]
                    / 100
                )
            else:
                Ea = df.loc[i, "vp_a"]

            Ew = 6.112 * np.exp(17.62 * T_f / (T_f + 243.12))
            Eice = 6.112 * np.exp(
                22.46 * (df.loc[i - 1, "T_s"]) / ((df.loc[i - 1, "T_s"]) + 243.12)
            )

            # atmospheric emissivity
            if df.loc[i, "Prec"] > 0:  # c = 1
                df.loc[i, "e_a"] = (
                    1.24 * math.pow(abs(Ea / (df.loc[i, "T_a"] + 273.15)), 1 / 7) * 1.22
                )
            else:
                df.loc[i, "e_a"] = (
                    1.24
                    * math.pow(abs(Ea / (df.loc[i, "T_a"] + 273.15)), 1 / 7)
                    * (1 + 0.22 * math.pow(c, 2))
                )

            # Fountain water output
            df.loc[i, "liquid"] = (
                df.loc[i, "Discharge"]
                * (1 - ftl)
                * time_steps
                / (60)
            )

            ''' When fountain run '''
            if df.loc[i, "liquid"] > 0:

                # Does water droplet cool down enough to nucleate
                m_d = rho_w * math.pi * d_f ** 3 / 6  # mass of droplet
                ar_d = math.pi * d_f ** 2  # SA of droplet
                E_d = (
                    cw
                    * rho_w
                    * df.loc[i, "p_a"]
                    / p0
                    * math.pow(k, 2)
                    * df.loc[i, "v_a"]
                    * (df.loc[i, "T_a"] - T_f)
                    / (np.log(z / z0mi) * np.log(z / z0hi))
                )
                df.loc[i - 1, "T_s"] = (E_d * ar_d * df.loc[i, "d_t"]) / (
                    m_d * cw
                ) + T_f

                if df.loc[i - 1, "T_s"] < 0:
                    if df.loc[i - 1, "T_s"] < T_nucleation:
                        df.loc[i, "solid"] += df.loc[i, "liquid"]
                        df.loc[i, "liquid"] = 0
                        df.loc[i - 1, "T_s"] = 0
                    else:
                        df.loc[i - 1, "T_s"] = 0  # No Nucleation

                # Evaporation or condensation
                df.loc[i, "Ql"] = (
                    0.623
                    * Le
                    * rho_a
                    / p0
                    * math.pow(k, 2)
                    * df.loc[i, "v_a"]
                    * (Ea - Ew)
                    / (np.log(z / z0mi) * np.log(z / z0hi))
                )

                df.loc[i, "gas"] -= (
                     (df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps) / Le
                )

                # Removing water quantity generated from previous time step
                df.loc[i, "liquid"] += (
                    (df.loc[i, "Ql"] * df.loc[i, "SA"] * time_steps) / Le
                )

                df.loc[i, "e_s"] = we

            else: # When fountain off
                # Sublimation or condensation
                df.loc[i, "Ql"] = (
                    0.623
                    * Ls
                    * rho_a
                    / p0
                    * math.pow(k, 2)
                    * df.loc[i, "v_a"]
                    * (Ea - Eice)
                    / (np.log(z / z0mi) * np.log(z / z0hi))
                )

                df.loc[i, 'gas'] -= (
                     (df.loc[i, 'Ql'] * (df.loc[i, 'SA']) * time_steps) / Ls
                )

                # Removing gas quantity generated from previous time step
                df.loc[i - 1, 'ice'] += (
                     (df.loc[i, 'Ql'] * (df.loc[i, 'SA']) * time_steps) / Ls
                )

                if df.loc[i - 1, "ice"] > 0:

                    '''Determine Ice Layers'''
                    if df.loc[i - 1, "ice"] < dx * df.loc[i, "SA"] * rho_i:
                        ice_layer = df.loc[i - 1, "ice"]
                    else:
                        ice_layer = dx * df.loc[i, "SA"] * rho_i

                    '''Melting Ice by layers'''
                    Ql = df.loc[i, "Ql"] * (df.loc[i, "SA"]) * time_steps
                    ctr = 0
                    #melted layers
                    ctr = int(abs(Ql)/  (ice_layer  * Lf))

                    if ctr != 0 :

                        # Subtract layers meltwater
                        Ql -= ctr * ice_layer  * Lf

                        df.loc[i, "solid"] -= (
                            (ctr * ice_layer)
                        )
                        df.loc[i, "melted"] += (
                            (ctr * ice_layer)
                        )

                    # Ice Temperature
                    df.loc[i, "T_s"] = (Ql) / (
                        ice_layer * ci
                    ) + df.loc[i - 1, "T_s"]

                    '''Hot Ice'''
                    if df.loc[i, "T_s"] > 0:
                        # Melting Ice by Temperature
                        df.loc[i, "solid"] -= (
                            (ice_layer * ci)
                            * (-df.loc[i, "T_s"])
                            / (-Lf)
                        )
                        df.loc[i, "melted"] += (
                            (ice_layer * ci)
                            * (-df.loc[i, "T_s"])
                            / (-Lf)
                        )
                        df.loc[i, "T_s"] = 0

                else:
                    df.loc[i, "ice"] = 0
                    df.loc[i, "gas"] = df.loc[i - 1, "ice"]

            # Short Wave Radiation SW
            df.loc[i, "SW"] = (1 - df.loc[i, "a"]) * (
                df.loc[i, "Rad"] * df.loc[i, "SRf"] + df.loc[i, "DRad"]
            )

            # Long Wave Radiation LW
            if "oli000z0" not in list(df.columns):
                if df.loc[i - 1, "T_s"] < -100 :
                    df.loc[i - 1, "T_s"] = -100
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
                / (np.log(z / z0mi) * np.log(z / z0hi))
            )

            # Total Energy W/m2
            df.loc[i, "TotalE"] = df.loc[i, "SW"] + df.loc[i, "LW"] + df.loc[i, "Qs"]

            # Total Energy Joules
            df.loc[i, "EJoules"] = df.loc[i, "TotalE"] * time_steps * df.loc[i, "SA"]

            if df.loc[i - 1, "ice"] > 0:
                '''Determine Ice Layers'''
                if df.loc[i - 1, "ice"] < dx * df.loc[i, "SA"] * rho_i:
                    ice_layer = df.loc[i - 1, "ice"]
                else:
                    ice_layer = dx * df.loc[i, "SA"] * rho_i

                if df.loc[i, "EJoules"] < 0:  # Energy Negative
                    if (
                        df.loc[i - 1, "liquid"] > 0
                    ):
                        '''Freezing water'''
                        df.loc[i, "liquid"] -= (df.loc[i, "EJoules"]) / (
                            -Lf + (df.loc[i - 1, "liquid"] * cw) * (-df.loc[i - 1, "T_s"])
                        )

                        # if df.loc[i, "liquid"] < 0:
                        #     df.loc[i, "liquid"] += (df.loc[i, "EJoules"]) / (
                        #         -Lf + (df.loc[i - 1, "liquid"] * cw) * (-df.loc[i - 1, "T_s"])
                        #     )
                        #     df.loc[i, "solid"] += df.loc[i, "liquid"]
                        #     df.loc[i, "liquid"] = 0
                        # else :
                        df.loc[i, "solid"] += (df.loc[i, "EJoules"]) / (
                            -Lf + (df.loc[i - 1, "liquid"] * cw) * (-df.loc[i - 1, "T_s"])
                        )

                    else:
                        # Cooling Ice
                        df.loc[i, "T_s"] = (df.loc[i, "EJoules"]) / (
                            ice_layer * ci
                        ) + df.loc[i - 1, "T_s"]

                else: # Energy Positive

                    if df.loc[i - 1, "ice"] > 0:

                        '''Melting Ice by layers'''
                        ctr = 0
                        #melted layers
                        ctr = int(abs(df.loc[i, "EJoules"])/  (ice_layer  * Lf))

                        if ctr != 0 :

                            # Subtract layers meltwater
                            df.loc[i, "EJoules"] -= ctr * ice_layer  * Lf

                            df.loc[i, "solid"] -= (
                                (ctr * ice_layer)
                            )
                            df.loc[i, "melted"] += (
                                (ctr * ice_layer)
                            )

                        # Heating Ice
                        df.loc[i, "T_s"] = (df.loc[i, "EJoules"]) / (
                            ice_layer * ci
                        ) + df.loc[i - 1, "T_s"]


                        if df.loc[i, "T_s"] > 0 :
                            # Melting Ice by Temperature
                            df.loc[i, "solid"] -= (
                                (ice_layer * ci)
                                * (-df.loc[i, "T_s"])
                                / (-Lf)
                            )
                            df.loc[i, "melted"] += (
                                (ice_layer * ci)
                                * (-df.loc[i, "T_s"])
                                / (-Lf)
                            )
                            df.loc[i, "T_s"] = 0

            ''' Quantities of all phases '''
            df.loc[i, "water"] = df.loc[i - 1, "water"] + df.loc[i, "liquid"]
            df.loc[i, "meltwater"] = df.loc[i - 1, "meltwater"] + df.loc[i, "melted"]
            df.loc[i, "ice"] = df.loc[i - 1, "ice"] + df.loc[i, "solid"]
            df.loc[i, "vapour"] = df.loc[i - 1, "vapour"] + df.loc[i, "gas"]
            df.loc[i, "iceV"] = df.loc[i, "ice"] / rho_i

            #print(df.loc[i, "When"], df.loc[i, "solid"], df.loc[i, "melted"])

    df = df[start:i]

    print("Ice Mass Remaining", float(df["ice"].tail(1)))
    print("Meltwater", float(df["meltwater"].tail(1)))
    print("Evaporated/sublimated", float(df["vapour"].tail(1)))
    print("Ice Volume Max", float(df["iceV"].max()))
    print("Ppt", prec)
    print("Model ended", df.loc[i - 1, "When"])
    print("Model runtime", df.loc[i - 1, "When"] - df.loc[start, "When"])

    return df
