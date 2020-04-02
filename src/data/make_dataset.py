import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
from pandas.plotting import register_matplotlib_converters
import math
import time
from pathlib import Path
import os
import logging
from tqdm import tqdm
from src.data.config import site, dates, option, folders, fountain, surface

start = time.time()

def projectile_xy(v, hs=0.0, g=9.8):
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
    theta_f = math.radians(45)
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
    return x

def getSEA(date, latitude, longitude, utc_offset):
    hour = date.hour
    minute = date.minute
    # Check your timezone to add the offset
    hour_minute = (hour + minute / 60) - utc_offset
    day_of_year = date.timetuple().tm_yday

    g = (360 / 365.25) * (day_of_year + hour_minute / 24)

    g_radians = math.radians(g)

    declination = (
        0.396372
        - 22.91327 * math.cos(g_radians)
        + 4.02543 * math.sin(g_radians)
        - 0.387205 * math.cos(2 * g_radians)
        + 0.051967 * math.sin(2 * g_radians)
        - 0.154527 * math.cos(3 * g_radians)
        + 0.084798 * math.sin(3 * g_radians)
    )

    time_correction = (
        0.004297
        + 0.107029 * math.cos(g_radians)
        - 1.837877 * math.sin(g_radians)
        - 0.837378 * math.cos(2 * g_radians)
        - 2.340475 * math.sin(2 * g_radians)
    )

    SHA = (hour_minute - 12) * 15 + longitude + time_correction

    if SHA > 180:
        SHA_corrected = SHA - 360
    elif SHA < -180:
        SHA_corrected = SHA + 360
    else:
        SHA_corrected = SHA

    lat_radians = math.radians(latitude)
    d_radians = math.radians(declination)
    SHA_radians = math.radians(SHA)

    SZA_radians = math.acos(
        math.sin(lat_radians) * math.sin(d_radians)
        + math.cos(lat_radians) * math.cos(d_radians) * math.cos(SHA_radians)
    )

    SZA = math.degrees(SZA_radians)

    SEA = 90 - SZA

    if SEA < 0:  # Before Sunrise or after sunset
        SEA = 0

    return math.radians(SEA)

def discharge_rate(df, fountain):

    if option == 'schwarzsee':
        df["Fountain"] = 0  # Fountain run time

        df_nights = pd.read_csv(
            os.path.join(folders["dirname"], "data/raw/schwarzsee_fountain_time.txt"),
            sep="\s+",
        )

        df_nights["Start"] = pd.to_datetime(
            df_nights["Date"] + " " + df_nights["start"]
        )
        df_nights["End"] = pd.to_datetime(df_nights["Date"] + " " + df_nights["end"])
        df_nights["Start"] = pd.to_datetime(
            df_nights["Start"], format="%Y-%m-%d %H:%M:%S"
        )
        df_nights["End"] = pd.to_datetime(df_nights["End"], format="%Y-%m-%d %H:%M:%S")

        df_nights["Date"] = pd.to_datetime(df_nights["Date"], format="%Y-%m-%d")
        mask = (df_nights["Date"] >= dates["start_date"]) & (
                df_nights["Date"] <= dates["end_date"]
        )
        df_nights = df_nights.loc[mask]
        df_nights = df_nights.reset_index()

        for i in range(0, df_nights.shape[0]):
            df_nights.loc[i, "Start"] = df_nights.loc[i, "Start"] - pd.Timedelta(days=1)
            df.loc[
                (df["When"] >= df_nights.loc[i, "Start"])
                & (df["When"] <= df_nights.loc[i, "End"]),
                "Fountain",
            ] = 1

    if option == 'temperature':
        mask = df["T_a"] < fountain["crit_temp"]
        mask_index = df[mask].index
        df.loc[mask_index, "Fountain"] = 1
        mask = df["When"] >= dates["fountain_off_date"]
        mask_index = df[mask].index
        df.loc[mask_index, "Fountain"] = 0

    if option == "energy":

        """Constants"""
        Ls = 2848 * 1000  # J/kg Sublimation
        Le = 2514 * 1000  # J/kg Evaporation
        Lf = 334 * 1000  # J/kg Fusion
        cw = 4.186 * 1000  # J/kg Specific heat water
        ci = 2.108 * 1000  # J/kgC Specific heat ice
        rho_w = 1000  # Density of water
        rho_i = 916  # Density of Ice rho_i
        rho_a = 1.29  # kg/m3 air density at mean sea level
        k = 0.4  # Van Karman constant
        bc = 5.670367 * math.pow(10, -8)  # Stefan Boltzman constant
        g = 9.8  # gravity

        """Miscellaneous"""
        time_steps = 5 * 60  # s Model time steps
        p0 = 1013  # Standard air pressure hPa

        """ Estimating Albedo """
        df["a"] = albedo(df, surface)

        df["T_s"] = 0

        for i in range(1, df.shape[0]):

            """ Energy Balance starts """

            # Vapor Pressure empirical relations
            if "vp_a" not in list(df.columns):
                df.loc[i, "vp_a"] = (
                        6.11
                        * math.pow(
                    10, 7.5 * df.loc[i - 1, "T_a"] / (df.loc[i - 1, "T_a"] + 237.3)
                )
                        * df.loc[i, "RH"]
                        / 100
                )

            df.loc[i, "vp_ice"] = 6.112 * np.exp(
                22.46 * (df.loc[i - 1, "T_s"]) / ((df.loc[i - 1, "T_s"]) + 243.12)
            )

            # Sublimation only
            df.loc[i, "Ql"] = (
                    0.623
                    * Ls
                    * rho_a
                    / p0
                    * math.pow(k, 2)
                    * df.loc[i, "v_a"]
                    * (df.loc[i, "vp_a"] - df.loc[i, "vp_ice"])
                    / (
                            np.log(surface["h_aws"] / surface["z0mi"])
                            * np.log(surface["h_aws"] / surface["z0hi"])
                    )
            )

            # Short Wave Radiation SW
            df.loc[i, "SW"] = (1 - df.loc[i, "a"]) * (
                    df.loc[i, "Rad"] + df.loc[i, "DRad"]
            )

            # Cloudiness from diffuse fraction
            if df.loc[i, "Rad"] + df.loc[i, "DRad"] > 1:
                df.loc[i, "cld"] = df.loc[i, "DRad"] / (
                        df.loc[i, "Rad"] + df.loc[i, "DRad"]
                )
            else:
                df.loc[i, "cld"] = 0

            # atmospheric emissivity
            df.loc[i, "e_a"] = (
                                           1.24
                                           * math.pow(abs(df.loc[i, "vp_a"] / (df.loc[i, "T_a"] + 273.15)),
                                                      1 / 7)
                                   ) * (1 + 0.22 * math.pow(df.loc[i, "cld"], 2))

            # Long Wave Radiation LW
            if "oli000z0" not in list(df.columns):

                df.loc[i, "LW"] = df.loc[i, "e_a"] * bc * math.pow(
                    df.loc[i, "T_a"] + 273.15, 4
                ) - surface["ie"] * bc * math.pow(df.loc[i - 1, "T_s"] + 273.15, 4)
            else:
                df.loc[i, "LW"] = df.loc[i, "oli000z0"] - surface["ie"] * bc * math.pow(
                    df.loc[i - 1, "T_s"] + 273.15,
                    4)

            # Sensible Heat Qs
            df.loc[i, "Qs"] = (
                    ci
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

            # Total Energy W/m2
            df.loc[i, "TotalE"] = df.loc[i, "SW"] + df.loc[i, "LW"] + df.loc[i, "Qs"] + df.loc[i, "Ql"]

        x = df["When"]
        mask = df["TotalE"] < 0
        mask_index = df[mask].index
        df.loc[mask_index, "Fountain"] = 1
        mask = df["When"] >= dates["fountain_off_date"]
        mask_index = df[mask].index
        df.loc[mask_index, "Fountain"] = 0

    if option == 'schwarzsee':
        df.Discharge = df.Discharge * df.Fountain
    else:
        df.Discharge = fountain["discharge"] * df.Fountain

    return df["Fountain"], df["Discharge"]

def emissivity(df):
    S_0 = 1360  # W m-2

    for i in range(1, df.shape[0]) :

        # Vapor Pressure empirical relations
        if "vp_a" not in list(df.columns):
            df.loc[i, "vp_a"] = (
                    6.11
                    * math.pow(
                10, 7.5 * df.loc[i - 1, "T_a"] / (df.loc[i - 1, "T_a"] + 237.3)
            )
                    * df.loc[i, "RH"]
                    / 100
            )

        # Estimating Solar Area fraction
        df.loc[i, "theta_s"] = getSEA(
            df.loc[i, "When"],
            fountain["latitude"],
            fountain["longitude"],
            fountain["utc_offset"],
        )

        optical_air_mass = 35 * math.sin(df.loc[i, "theta_s"] * (math.pow((1224 * math.pow(math.sin(df.loc[i, "theta_s"]), 2) + 1), -0.5)))
        tau_R_tau_pg = 1.021 - 0.084* math.pow((optical_air_mass * (0.00949 * df.loc[i, "p_a"] + 0.051)), 0.5)
        precipitable_water = 4650 * df.loc[i, "vp_a"] / (df.loc[i, "T_a"] + 273.15)
        # print(optical_air_mass, precipitable_water)
        tau_w = 1 - 0.077 * math.pow((optical_air_mass * precipitable_water),0.3)
        tau_a = math.pow(0.935, optical_air_mass)
        df.loc[i, "S_clr"] = S_0 * math.sin(df.loc[i, "theta_s"]) * tau_R_tau_pg * tau_w * tau_a

        if df.loc[i, "S_clr"] != 0:
            df.loc[i, "s"] = (df.loc[i, "Rad"] + df.loc[i, "DRad"]) / df.loc[i, "S_clr"]
            if df.loc[i, "s"] > 1:
                df.loc[i, "s"] = 1
            df.loc[i, "e_a"] = (1 - df.loc[i, "s"]) + df.loc[i, "s"] * (0.83 - 0.18 * math.pow( 10, -0.067 * df.loc[i, "vp_a"]))
        else:
            if i - 96 > 0:
                for j in range(i - 96, i):
                    df.loc[i, "e_a"] += df.loc[j, "e_a"]
            else:
                for j in range(i, i + 96):
                    df.loc[i, "e_a"] += df.loc[j, "e_a"]

                    df.loc[i, "e_a"] = df.loc[i, "e_a"] / 96

        # # Cloudiness from diffuse fraction
        # if df.loc[i, "Rad"] + df.loc[i, "DRad"] > 1:
        #     df.loc[i, "cld"] = df.loc[i, "DRad"] / (
        #             df.loc[i, "Rad"] + df.loc[i, "DRad"]
        #     )
        # else:
        # # Night Cloudiness average of last 8 hours
        #     if i - 96 > 0:
        #         for j in range(i - 96, i):
        #             df.loc[i, "cld"] += df.loc[j, "cld"]
        #     else:
        #         for j in range(i, i + 96):
        #             df.loc[i, "cld"] += df.loc[j, "cld"]
        #
        #             df.loc[i, "cld"] = df.loc[i, "cld"] / 96

    return df["vp_a"], df["theta_s"], df["cld"]


if __name__ == '__main__':

    if site == "schwarzsee":

        # read files
        df_in = pd.read_csv(
            folders["data_file"],
            header=None,
            encoding="latin-1",
            skiprows=7,
            sep="\s+",
            names=[
                "Date",
                "Time",
                "Discharge",
                "Wind Direction",
                "Wind Speed",
                "Maximum Wind Speed",
                "Temperature",
                "Humidity",
                "Pressure",
                "Pluviometer",
            ],
        )

        # Drop
        df_in = df_in.drop(["Pluviometer"], axis=1)

        # Add Radiation data
        df_in2 = pd.read_csv(
            os.path.join(folders["dirname"], "data/raw/plaffeien_rad.txt"),
            sep="\s+",
            skiprows=2,
        )
        df_in2["When"] = pd.to_datetime(df_in2["time"], format="%Y%m%d%H%M")
        df_in2["ods000z0"] = pd.to_numeric(df_in2["ods000z0"], errors="coerce")
        df_in2["gre000z0"] = pd.to_numeric(df_in2["gre000z0"], errors="coerce")
        df_in2["Rad"] = df_in2["gre000z0"] - df_in2["ods000z0"]
        df_in2["DRad"] = df_in2["ods000z0"]

        # Add Precipitation data
        df_in2["Prec"] = pd.to_numeric(df_in2["rre150z0"], errors="coerce")
        df_in2["Prec"] = df_in2["Prec"] / 2  # 5 minute sum
        df_in2 = df_in2.set_index("When").resample("5T").ffill().reset_index()

        # Datetime
        df_in["When"] = pd.to_datetime(df_in["Date"] + " " + df_in["Time"])
        df_in["When"] = pd.to_datetime(df_in["When"], format="%Y.%m.%d %H:%M:%S")

        # Correct data errors
        i = 1
        while df_in.loc[i, "When"] != datetime(2019, 2, 6, 16, 15):
            if str(df_in.loc[i, "When"].year) != "2019":
                df_in.loc[i, "When"] = df_in.loc[i - 1, "When"] + pd.Timedelta(minutes=5)
            i = i + 1

        while df_in.loc[i, "When"] != datetime(2019, 3, 2, 15):
            if str(df_in.loc[i, "When"].year) != "2019":
                df_in.loc[i, "When"] = df_in.loc[i - 1, "When"] + pd.Timedelta(minutes=5)
            i = i + 1

        while df_in.loc[i, "When"] != datetime(2019, 3, 6, 16, 25):
            if str(df_in.loc[i, "When"].year) != "2019":
                df_in.loc[i, "When"] = df_in.loc[i - 1, "When"] + pd.Timedelta(minutes=5)
            i = i + 1

        df_in = df_in.resample("5Min", on="When").first().drop("When", 1).reset_index()

        # Fill missing data
        for i in range(1, df_in.shape[0]):
            if np.isnan(df_in.loc[i, "Temperature"]):
                df_in.loc[i, "Temperature"] = df_in.loc[i - 288, "Temperature"]
                df_in.loc[i, "Humidity"] = df_in.loc[i - 288, "Humidity"]
                df_in.loc[i, "Wind Speed"] = df_in.loc[i - 288, "Wind Speed"]
                df_in.loc[i, "Maximum Wind Speed"] = df_in.loc[
                    i - 288, "Maximum Wind Speed"
                ]
                df_in.loc[i, "Wind Direction"] = df_in.loc[i - 288, "Wind Direction"]
                df_in.loc[i, "Pressure"] = df_in.loc[i - 288, "Pressure"]
                df_in.loc[i, "Discharge"] = df_in.loc[i - 288, "Discharge"]

        mask = (df_in["When"] >= dates["start_date"]) & (df_in["When"] <= dates["end_date"])
        df_in = df_in.loc[mask]
        df_in = df_in.reset_index()

        mask = (df_in2["When"] >= dates["start_date"]) & (
            df_in2["When"] <= dates["end_date"]
        )
        df_in2 = df_in2.loc[mask]
        df_in2 = df_in2.reset_index()

        days = pd.date_range(start=dates["start_date"], end=dates["end_date"], freq="5T")
        days = pd.DataFrame({"When": days})

        df = pd.merge(
            days,
            df_in[
                [
                    "When",
                    "Discharge",
                    "Wind Speed",
                    "Maximum Wind Speed",
                    "Wind Direction",
                    "Temperature",
                    "Humidity",
                    "Pressure",
                ]
            ],
            on="When",
        )

        # Add Radiation DataFrame
        df["Rad"] = df_in2["Rad"]
        df["DRad"] = df_in2["DRad"]
        df["Prec"] = df_in2["Prec"] / 1000

        mask = (df["When"] >= dates["start_date"]) & (df["When"] <= dates["end_date"])
        df = df.loc[mask]
        df = df.reset_index()

        df["cld"] = 0
        df["SEA"] = 0
        df["e_a"] = 0

        # CSV output
        df.rename(
            columns={
                "Wind Speed": "v_a",
                "Temperature": "T_a",
                "Humidity": "RH",
                "Volume": "V",
                "Pressure": "p_a",
            },
            inplace=True,
        )


        for i in range(1, df.shape[0]):
            if np.isnan(df.loc[i, "Rad"]):
                df.loc[i, "Rad"] = df.loc[i - 1, "Rad"]
            if np.isnan(df.loc[i, "DRad"]):
                df.loc[i, "DRad"] = df.loc[i - 1, "DRad"]
            if np.isnan(df.loc[i, "Prec"]):
                df.loc[i, "Prec"] = df.loc[i - 1, "Prec"]

            """Solar Elevation Angle"""
            df.loc[i, "SEA"] = getSEA(
                df.loc[i, "When"],
                fountain["latitude"],
                fountain["longitude"],
                fountain["utc_offset"],
            )

            """Cloudiness"""
            # Cloudiness from diffuse fraction
            if df.loc[i, "Rad"] + df.loc[i, "DRad"] > 1:
                df.loc[i, "cld"] = df.loc[i, "DRad"] / (
                        df.loc[i, "Rad"] + df.loc[i, "DRad"]
                )
            else:
                # Night Cloudiness average of last 8 hours
                if i - 96 > 0:
                    for j in range(i - 96, i):
                        df.loc[i, "cld"] += df.loc[j, "cld"]
                    df.loc[i, "cld"] = df.loc[i, "cld"] / 96
                else:
                    for j in range(0, i):
                        df.loc[i, "cld"] += df.loc[j, "cld"]
                    df.loc[i, "cld"] = df.loc[i, "cld"] / i

            """ Vapour Pressure"""
            if "vpa" not in list(df.columns):
                df.loc[i, "vp_a"] = (6.11 * math.pow(10, 7.5 * df.loc[i - 1, "T_a"] / (df.loc[i - 1, "T_a"] + 237.3)) * df.loc[i, "RH"] / 100)
            else:
                df.loc[i, "vp_a"] = df.loc[i, "vpa"]

            # atmospheric emissivity
            df.loc[i, "e_a"] = ( 1.24 * math.pow(abs(df.loc[i, "vp_a"] / (df.loc[i, "T_a"] + 273.15)), 1 / 7)
                               ) * (1 + 0.22 * math.pow(df.loc[i, "cld"], 2))

        # v_a mean
        v_a = df["v_a"].replace(0, np.NaN).mean()  # m/s Average Humidity
        df["v_a"] = df["v_a"].replace(0, v_a)

        df_out = df[
            ["When", "T_a", "RH", "v_a", "Discharge", "Rad", "DRad", "Prec", "p_a", "SEA", "cld", "vp_a", "e_a"]
        ]

        df_out = df_out.round(5)

    if site == "plaffeien":

        """
        Parameter
        ---------
                  Unit                                 Description
        pva200s0  hPa                                  Vapour pressure 2 m above ground; current value
        prestas0  hPa                                  Pressure at station level (QFE); current value
        ods000z0  W/m²                                 diffuse radiation, average 10 minutes
        gre000z0  W/m²                                 Global radiation; ten minutes mean
        tre200s0  °C                                   Air temperature 2 m above ground; current value
        rre150z0  mm                                   Precipitation; ten minutes total
        ure200s0  %                                    Relative air humidity 2 m above ground; current value
        fkl010z0  m/s                                  Wind speed scalar; ten minutes mean
        """

        # read files
        df_in = pd.read_csv(folders["data_file"], encoding="latin-1", skiprows=2, sep=";")
        df_in["When"] = pd.to_datetime(df_in["time"], format="%Y%m%d%H%M")  # Datetime

        time_steps = 5 * 60  # s # Model time steps
        mask = (df_in["When"] >= dates["start_date"]) & (df_in["When"] <= dates["end_date"])
        df_in = df_in.loc[mask]
        df_in = df_in.reset_index()

        # Add Radiation data
        df_in["ods000z0"] = pd.to_numeric(df_in["ods000z0"], errors="coerce")
        df_in["gre000z0"] = pd.to_numeric(df_in["gre000z0"], errors="coerce")
        df_in["Rad"] = df_in["gre000z0"] - df_in["ods000z0"]
        df_in["DRad"] = df_in["ods000z0"]
        df_in["T_a"] = pd.to_numeric(
            df_in["tre200s0"], errors="coerce"
        )  # Add Temperature data
        df_in["Prec"] = pd.to_numeric(
            df_in["rre150z0"], errors="coerce"
        )  # Add Precipitation data
        df_in["RH"] = pd.to_numeric(df_in["ure200s0"], errors="coerce")  # Add Humidity data
        df_in["v_a"] = pd.to_numeric(
            df_in["fkl010z0"], errors="coerce"
        )  # Add wind speed data
        df_in["p_a"] = pd.to_numeric(df_in["prestas0"], errors="coerce")  # Air pressure
        df_in["vpa"] = pd.to_numeric(
            df_in["pva200s0"], errors="coerce"
        )  # Vapour pressure over air

        df_in["Prec"] = df_in["Prec"] / 1000

        # Fill nans
        df_in = df_in.fillna(method="ffill")

        df_out = df_in[["When", "T_a", "RH", "v_a", "Rad", "DRad", "Prec", "p_a", "vpa",]]

        # 5 minute sum
        cols = ["T_a", "RH", "v_a", "Rad", "DRad", "Prec", "p_a", "vpa"]
        df_out[cols] = df_out[cols] / 2
        df_out = df_out.set_index("When").resample("5T").ffill().reset_index()

        for i in range(1, df_out.shape[0]):

            """Solar Elevation Angle"""
            df_out.loc[i, "SEA"] = getSEA(
                df_out.loc[i, "When"],
                fountain["latitude"],
                fountain["longitude"],
                fountain["utc_offset"],
            )


        cols = [
            "When",
            "T_a",
            "RH",
            "v_a",
            "Rad",
            "DRad",
            "Prec",
            "p_a",
            "vpa",
            "SEA",
        ]
        df_out = df_out[cols]
        df_out = df_out.round(5)

    if site == "guttannen":

        """
        Parameter
        ---------
                  Unit                                 Description
        pva200s0  hPa                                  Vapour pressure 2 m above ground; current value
        prestas0  hPa                                  Pressure at station level (QFE); current value
        gre000z0  W/m²                                 Global radiation; ten minutes mean
        oli000z0  W/m²                                 Longwave incoming radiation; ten minute average
        tre200s0  °C                                   Air temperature 2 m above ground; current value
        rre150z0  mm                                   Precipitation; ten minutes total
        ure200s0  %                                    Relative air humidity 2 m above ground; current value
        fkl010z0  m/s                                  Wind speed scalar; ten minutes mean
        """

        # read files
        df_in = pd.read_csv(folders["data_file"], encoding="latin-1", skiprows=2, sep=";")
        df_in["When"] = pd.to_datetime(df_in["time"], format="%Y%m%d%H%M")  # Datetime

        mask = (df_in["When"] >= dates["start_date"]) & (df_in["When"] <= dates["end_date"])
        df_in = df_in.loc[mask]
        df_in = df_in.reset_index()

        # Convert to int
        df_in["oli000z0"] = pd.to_numeric(
            df_in["oli000z0"], errors="coerce"
        )  # Add Longwave Radiation data
        df_in["gre000z0"] = pd.to_numeric(
            df_in["gre000z0"], errors="coerce"
        )  # Add Radiation data
        df_in["T_a"] = pd.to_numeric(
            df_in["tre200s0"], errors="coerce"
        )  # Add Temperature data
        df_in["Prec"] = pd.to_numeric(
            df_in["rre150z0"], errors="coerce"
        )  # Add Precipitation data
        df_in["RH"] = pd.to_numeric(df_in["ure200s0"], errors="coerce")  # Add Humidity data
        df_in["v_a"] = pd.to_numeric(
            df_in["fkl010z0"], errors="coerce"
        )  # Add wind speed data
        df_in["p_a"] = pd.to_numeric(df_in["prestas0"], errors="coerce")  # Air pressure
        df_in["vpa"] = pd.to_numeric(
            df_in["pva200s0"], errors="coerce"
        )  # Vapour pressure over air

        df_in["Rad"] = df_in["gre000z0"] - df_in["gre000z0"] * 0.1
        df_in["DRad"] = df_in["gre000z0"] * 0.1
        df_in["LW"] = df_in["oli000z0"]
        df_in["Prec"] = df_in["Prec"] / 1000

        # Fill nans
        df_in = df_in.fillna(method="ffill")

        df_out = df_in[
            ["When", "T_a", "RH", "v_a", "Rad", "DRad", "oli000z0", "Prec", "p_a", "vpa",]
        ]
        df_out = df_out.round(5)

        # 5 minute sum
        cols = ["T_a", "RH", "v_a", "Rad", "DRad", "Prec", "p_a", "vpa", "oli000z0"]
        df_out[cols] = df_out[cols] / 2
        df_out = df_out.set_index("When").resample("5T").ffill().reset_index()

        cols = [
            "When",
            "T_a",
            "RH",
            "v_a",
            "Rad",
            "DRad",
            "Prec",
            "p_a",
            "vpa",
        ]
        df_out = df_out[cols]
        df_out = df_out.round(5)


    filename = folders["input_folder"] + site

    df_out.to_csv(filename + "_raw_input.csv")

    """ Derived Parameters"""

    l = [
        "a",
        "r_f",
        "Fountain",
    ]
    for col in l:
        df[col] = 0

    """Discharge Rate"""
    df["Fountain"], df["Discharge"] = discharge_rate(df,fountain)

    """Albedo Decay"""
    surface["decay_t"] = (
        surface["decay_t"] * 24 * 60 / 5
    )  # convert to 5 minute time steps
    s = 0
    f = 0

    """ Fountain Spray radius """
    Area = math.pi * math.pow(fountain["aperture_f"], 2) / 4


    for i in tqdm(range(1, df.shape[0])):

        if option == "schwarzsee":

            ti = surface["decay_t"]
            a_min = surface["a_i"]

            # Precipitation
            if (df.loc[i, "Fountain"] == 0) & (df.loc[i, "Prec"] > 0):
                if df.loc[i, "T_a"] < surface["rain_temp"]:  # Snow
                    s = 0
                    f = 0

            if df.loc[i, "Fountain"] > 0:
                f = 1
                s = 0

            if f == 0:  # last snowed
                df.loc[i, "a"] = a_min + (surface["a_s"] - a_min) * math.exp(-s / ti)
                s = s + 1
            else:  # last sprayed
                df.loc[i, "a"] = a_min
                s = s + 1
        else:
            df.loc[i, "a"] = surface["a_i"]

        """ Fountain Spray radius """
        v_f = df.loc[i, "Discharge"] / (60 * 1000 * Area)
        df.loc[i, "r_f"] = projectile_xy(
            v_f, fountain["h_f"]
        )

    df.to_csv(filename + "_input.csv")

    total = time.time() - start

    print("Total time : ", total / 60)

    # filename = folders["input_folder"] + site + "_input.csv"
    # df = pd.read_csv(filename, sep=",")
    # df["When"] = pd.to_datetime(df["When"], format="%Y.%m.%d %H:%M:%S")

    pp = PdfPages(folders["input_folder"] + site + "_derived_parameters" + ".pdf")

    x = df.When

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y1 = df.Discharge
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Discharge [$l\, min^{-1}$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()


    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y1 = df.r_f
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Spray Radius [$m$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y1 = df.e_a
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Atmospheric emissivity")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    y1 = df.cld
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Cloudiness")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    y1 = df.vp_a
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Vapour Pressure")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    y1 = df.a
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Albedo")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    pp.close()

    """Input Plots"""

    pp = PdfPages(folders["input_folder"] + site + "_data" + ".pdf")

    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(
        nrows=6, ncols=1, sharex="col", sharey="row", figsize=(15, 12)
    )

    # fig.suptitle("Field Data", fontsize=14)
    # Remove horizontal space between axes
    # fig.subplots_adjust(hspace=0)

    x = df.When

    if option == "schwarzsee":
        y1 = df.Discharge
    else:
        y1 = df.Fountain
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Discharge [$l\, min^{-1}$]")
    ax1.grid()

    ax1t = ax1.twinx()
    ax1t.plot(x, df.Prec * 1000, "b-", linewidth=0.5)
    ax1t.set_ylabel("Precipitation [$mm$]", color="b")
    for tl in ax1t.get_yticklabels():
        tl.set_color("b")

    y2 = df.T_a
    ax2.plot(x, y2, "k-", linewidth=0.5)
    ax2.set_ylabel("Temperature [$\degree C$]")
    ax2.grid()


    y3 = df.Rad + df.DRad
    ax3.plot(x, y3, "k-", linewidth=0.5)
    ax3.set_ylabel("Global [$W\,m^{-2}$]")
    ax3.grid()

    ax3t = ax3.twinx()
    ax3t.plot(x, df.DRad, "b-", linewidth=0.5)
    ax3t.set_ylim(ax3.get_ylim())
    ax3t.set_ylabel("Diffuse [$W\,m^{-2}$]", color="b")
    for tl in ax3t.get_yticklabels():
        tl.set_color("b")

    y4 = df.RH
    ax4.plot(x, y4, "k-", linewidth=0.5)
    ax4.set_ylabel("Humidity [$\%$]")
    ax4.grid()

    y5 = df.p_a
    ax5.plot(x, y5, "k-", linewidth=0.5)
    ax5.set_ylabel("Pressure [$hPa$]")
    ax5.grid()

    y6 = df.v_a
    ax6.plot(x, y6, "k-", linewidth=0.5)
    ax6.set_ylabel("Wind [$m\,s^{-1}$]")
    ax6.grid()

    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())

    # rotates and right aligns the x labels, and moves the bottom of the axes up to make room for them
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")

    plt.savefig(
        os.path.join(folders["input_folder"], site + "_data.jpg"),
        bbox_inches="tight",
        dpi=300,
    )

    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y1 = df.T_a
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Temperature [$\degree C$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)


    y2 = df.Discharge
    ax1.plot(x, y2, "k-", linewidth=0.5)
    ax1.set_ylabel("Discharge Rate ")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y3 = df.Rad
    ax1.plot(x, y3, "k-", linewidth=0.5)
    ax1.set_ylabel("Direct SWR [$W\,m^{-2}$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y31 = df.DRad
    ax1.plot(x, y31, "k-", linewidth=0.5)
    ax1.set_ylabel("Diffuse SWR [$W\,m^{-2}$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y4 = df.Prec * 1000
    ax1.plot(x, y4, "k-", linewidth=0.5)
    ax1.set_ylabel("Ppt [$mm$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y5 = df.p_a
    ax1.plot(x, y5, "k-", linewidth=0.5)
    ax1.set_ylabel("Pressure [$hPa$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y6 = df.v_a
    ax1.plot(x, y6, "k-", linewidth=0.5)
    ax1.set_ylabel("Wind [$m\,s^{-1}$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    pp.close()
