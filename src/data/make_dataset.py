import sys
sys.path.append('/home/surya/Programs/PycharmProjects/air_model')
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
from tqdm import tqdm
import os
import logging
from src.data.config import site, option, dates, folders, fountain
from scipy import stats

start = time.time()

def discharge_rate(df, fountain):

    if option == 'schwarzsee':
        df["Fountain"] = 0  # Fountain run time

        df_nights = pd.read_csv(
            os.path.join(folders["raw_folder"], "schwarzsee_fountain_time.txt"),
            sep="\\s+",
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
                    df.loc[i, "SW_direct"] + df.loc[i, "DRad"]
            )

            # Cloudiness from diffuse fraction
            if df.loc[i, "SW_direct"] + df.loc[i, "DRad"] > 1:
                df.loc[i, "cld"] = df.loc[i, "DRad"] / (
                        df.loc[i, "SW_direct"] + df.loc[i, "DRad"]
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


if __name__ == '__main__':

    if site == "schwarzsee":

        # read files
        df_in = pd.read_csv(
            folders["raw_folder"]+ site + "_aws.txt",
            header=None,
            encoding="latin-1",
            skiprows=7,
            sep="\\s+",
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

        # Datetime
        df_in["When"] = pd.to_datetime(df_in["Date"] + " " + df_in["Time"])
        df_in["When"] = pd.to_datetime(df_in["When"], format="%Y.%m.%d %H:%M:%S")

        # Correct datetime errors
        for i in tqdm(range(1, df_in.shape[0])):
            if str(df_in.loc[i, "When"].year) != "2019":
                df_in.loc[i, "When"] = df_in.loc[i - 1, "When"] + pd.Timedelta(minutes=5)

        df_in = df_in.set_index("When").resample("5T").last().reset_index()

        mask = (df_in["When"] >= dates["start_date"]) & (df_in["When"] <= dates["end_date"])
        df_in = df_in.loc[mask]
        df_in = df_in.reset_index()

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

        # CSV output
        df.rename(
            columns={
                "Wind Speed": "v_a",
                "Temperature": "T_a",
                "Humidity": "RH",
                "Pressure": "p_a",
            },
            inplace=True,
        )

        # ERA5 begins
        df_in3 = pd.read_csv(folders["raw_folder"]+ site + "_ERA5.csv", sep=",", header=0, parse_dates=["dataDate"])

        df_in3 = df_in3.drop(["Latitude", "Longitude"], axis=1)
        df_in3["time"] = df_in3["validityTime"].replace([0, 100, 200, 300, 400, 500, 600, 700, 800, 900],
                                                ["0000", "0100", "0200", "0300", "0400", "0500", "0600", "0700", "0800",
                                                 "0900"])
        df_in3["time"] = df_in3["time"].astype(str)
        df_in3["time"] = df_in3["time"].str[0:2] + ":" + df_in3["time"].str[2:4]
        df_in3["dataDate"] = df_in3["dataDate"].astype(str)
        df_in3["When"] = df_in3["dataDate"] + " " + df_in3["time"]
        df_in3["When"] = pd.to_datetime(df_in3["When"])

        days = pd.date_range(start=dates["start_date"], end=df_in3["When"].iloc[-1], freq="1H")
        df_out1 = pd.DataFrame({"When": days})

        df_out1 = df_out1.set_index("When")
        df_in3 = df_in3.set_index("When")

        time_steps = 60 * 60
        df_out1["10u"] = df_in3.loc[df_in3.shortName == '10u', "Value"]
        df_out1["10v"] = df_in3.loc[df_in3.shortName == '10v', "Value"]
        df_out1["2d"] = df_in3.loc[df_in3.shortName == '2d', "Value"]
        df_out1["2t"] = df_in3.loc[df_in3.shortName == '2t', "Value"]
        df_out1["sp"] = df_in3.loc[df_in3.shortName == 'sp', "Value"]
        df_out1["tcc"] = df_in3.loc[df_in3.shortName == 'tcc', "Value"]
        df_out1["tp"] = df_in3.loc[df_in3.shortName == 'tp', "Value"]
        df_out1["ssrd"] = df_in3.loc[df_in3.shortName == 'ssrd', "Value"] / time_steps
        df_out1["strd"] = df_in3.loc[df_in3.shortName == 'strd', "Value"] / time_steps
        df_out1["fdir"] = df_in3.loc[df_in3.shortName == 'fdir', "Value"] / time_steps

        df_out1["v_a"] = np.sqrt(df_out1["10u"] ** 2 + df_out1["10v"] ** 2)
        df_out1["RH"] = 100 * (np.exp((17.625 * df_out1["2d"]) / (243.04 + df_out1["2d"])) / np.exp(
            (17.625 * df_out1["2t"]) / (243.04 + df_out1["2t"])))
        df_out1["sp"] = df_out1["sp"] / 100
        df_out1["tp"] = df_out1["tp"] # mm/s
        df_out1["SW_diffuse"] = df_out1["ssrd"] - df_out1["fdir"]
        df_out1["2t"] = df_out1["2t"] - 273.15

        # CSV output
        df_out1.rename(
            columns={
                "2t": "T_a",
                "sp": "p_a",
                "tcc": "cld",
                "tp": "Prec",
                "fdir": 'SW_direct',
                "strd": 'LW_in',

            },
            inplace=True,
        )

        df_in3 = df_out1[
            ["T_a", "RH","Prec", "v_a", "SW_direct", "SW_diffuse", "LW_in", "cld", "p_a"]
        ]

        df_in3 = df_in3.round(5)

        upsampled = df_in3.resample("5T")
        interpolated = upsampled.interpolate(method='linear')
        interpolated = interpolated.reset_index()

        interpolated["Discharge"] = 0
        mask = (interpolated["T_a"] < fountain["crit_temp"]) & (interpolated["SW_direct"] < 100)
        mask_index = interpolated[mask].index
        interpolated.loc[mask_index, "Discharge"] = 2 * 60
        mask = interpolated["When"] >= dates["fountain_off_date"]
        mask_index = interpolated[mask].index
        interpolated.loc[mask_index, "Discharge"] = 0
        interpolated = interpolated.reset_index()

        df_in3 = interpolated[
            ["When", "T_a", "RH", "v_a", "SW_direct", "SW_diffuse", "LW_in", "cld", "p_a"]
        ]

        df_in3 = df_in3.reset_index()
        mask = (df_in3["When"] >= dates["start_date"]) & (df_in3["When"] <= dates["end_date"])
        df_in3 = df_in3.loc[mask]
        df_in3 = df_in3.reset_index()

        df_in3.to_csv(folders["input_folder"] + "raw_input_ERA5.csv")

        df_ERA5 = interpolated[["When", "T_a", "RH", "v_a", "SW_direct", "SW_diffuse", "LW_in", "cld", "p_a", "Prec", "Discharge"]]
        df_ERA5.loc[:,"Discharge"] = 0


        # Fill from ERA5
        df_ERA5 = df_ERA5.set_index("When")
        df = df.set_index("When")
        df.loc[df["T_a"].isnull(), [ 'T_a', 'RH', 'v_a', 'p_a', 'Discharge']] = df_ERA5[[ 'T_a', 'RH', 'v_a', 'p_a', 'Discharge']]

        df["v_a"] = df["v_a"].replace(0, np.NaN)
        df.loc[df["v_a"].isnull(), "v_a"] = df_ERA5["v_a"]

        df[['SW_direct', "SW_diffuse", 'cld']] = df_ERA5[['SW_direct', "SW_diffuse", 'cld']]

        # RMSE

        print("RMSE Temp", ((df.T_a - df_ERA5.T_a) ** 2).mean() ** .5)
        print("RMSE wind", ((df.v_a - df_ERA5.v_a) ** 2).mean() ** .5)

        slope, intercept, r_value1, p_value, std_err = stats.linregress(df.T_a.values , df_in3.T_a.values)
        slope, intercept, r_value2, p_value, std_err = stats.linregress(df.v_a.values , df_in3.v_a.values)

        print("R2 temp", r_value1**2)
        print("R2 wind", r_value2**2)

        # Add Precipitation data
        df_in2 = pd.read_csv(
            os.path.join(folders["raw_folder"], "plf_ppt.txt"),
            sep=";",
            skiprows=2,
        )
        df_in2["When"] = pd.to_datetime(df_in2["time"], format="%Y%m%d%H%M")

        df_in2["Prec"] = pd.to_numeric(df_in2["rre150z0"], errors="coerce")
        df_in2["Prec"] = df_in2["Prec"] / (10*60)  # ppt rate mm/s
        df_in2 = df_in2.set_index("When").resample("5T").interpolate(method='linear').reset_index()

        mask = (df_in2["When"] >= dates["start_date"]) & (
                df_in2["When"] <= dates["end_date"]
        )
        df_in3 = df_in2.loc[mask]
        df_in3 = df_in3.set_index("When")
        df["Prec"] = df_in3["Prec"]

        df = df.reset_index()


        """Discharge Rate"""
        df["Fountain"], df["Discharge"] = discharge_rate(df,fountain)

        df["Discharge"] = df["Discharge"] * df["Fountain"]

        df_out = df[
            ["When", "T_a", "RH", "v_a", "Discharge", "SW_direct", "SW_diffuse", "Prec", "p_a", 'cld']
        ]

        if df_out.isnull().values.any() :
            print( "Warning: Null values present")
            print(df[["When", "T_a", "RH", "v_a", "Discharge", "SW_direct", "SW_diffuse", "Prec", "p_a", 'cld']].isnull().sum())

        df_out = df_out.round(5)


        df_out.to_csv(folders["input_folder"] + "raw_input.csv")


        # Extend data
        df_ERA5["Prec"] = 0
        df_ERA5 = df_ERA5.reset_index()
        mask = (df_ERA5["When"] >= df_out["When"].iloc[-1])& (
                df_ERA5["When"] <= datetime(2019,5,30)
        )
        df_ERA5 = df_ERA5.loc[mask]

        mask = (df_in2["When"] >= dates["start_date"]) & (
                df_in2["When"] <= df_ERA5["When"].iloc[-1]
        )
        df_in2 = df_in2.loc[mask]
        df_in2 = df_in2.set_index("When")

        df_out = df_out.set_index("When")
        df_ERA5 = df_ERA5.set_index("When")

        df_ERA5 = df_ERA5.drop(["LW_in"], axis=1)
        concat = pd.concat([df_out, df_ERA5])
        concat["Prec"] = df_in2["Prec"]
        concat = concat.reset_index()

        if concat.isnull().values.any() :
            print( "Warning: Null values present")
            print(concat[["When", "T_a", "RH", "v_a", "Discharge", "SW_direct", "SW_diffuse", "Prec", "p_a", 'cld']].isnull().sum())

        concat.to_csv(folders["input_folder"] + "raw_input_extended.csv")
        concat.to_hdf(
            folders["input_folder"] + "raw_input_extended.h5", key="df", mode="w"
        )


    if site == "guttannen":

        # read files
        df_in = pd.read_csv(
            folders["data_file"],
            header=None,
            encoding="latin-1",
            skiprows=7,
            sep="\\s+",
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

        # Datetime
        df_in["When"] = pd.to_datetime(df_in["Date"] + " " + df_in["Time"])
        df_in["When"] = pd.to_datetime(df_in["When"], format="%Y.%m.%d %H:%M:%S")

        # Drop
        df_in = df_in.drop(["Pluviometer"], axis=1)
        df_in = df_in.drop(["Wind Direction"], axis=1)

        mask = (df_in["When"] >= dates["start_date"]) & (df_in["When"] <= dates["end_date"])
        df_in = df_in.loc[mask]
        df_in = df_in.reset_index()

        # Error Correction

        # v_a mean
        v_a = df_in["Wind Speed"].replace(5643.2, np.NaN).mean()  # m/s Average Humidity
        df_in["Wind Speed"] = df_in["Wind Speed"].replace(5643.2, v_a)

        """
                Parameter
                ---------
                          Unit                                 Description
                oli000z0
                prestas0  hPa                                  Pressure at station level (QFE); current value
                gre000z0  W/m²                                 Global radiation; ten minutes mean
                pva200s0  hPa                                  Vapour pressure 2 m above ground; current value
                rre150z0  mm                                   Precipitation; ten minutes total
                ure200s0  %                                    Relative air humidity 2 m above ground;
                fkl010z0  m/s                                  Wind speed scalar; ten minutes mean
                tre200s0  °C                                   Air temperature 2 m above ground; current
        """

        # Add Radiation data
        df_in2 = pd.read_csv(os.path.join(folders["dirname"], "data/raw/guttannen_1.txt"), encoding="latin-1", skiprows=2, sep='\\s+')
        df_in2["When"] = pd.to_datetime(df_in2["time"], format="%Y%m%d%H%M")  # Datetime

        # Convert to int
        df_in2["oli000z0"] = pd.to_numeric(
            df_in2["oli000z0"], errors="coerce"
        )  # Add Longwave Radiation data

        df_in2["gre000z0"] = pd.to_numeric(
            df_in2["gre000z0"], errors="coerce"
        )  # Add Radiation data

        # Add rest of data
        df_in3 = pd.read_csv(os.path.join(folders["dirname"], "data/raw/guttannen_2.txt"), encoding="latin-1",
                             skiprows=2, sep='\\s+')
        df_in3["When"] = pd.to_datetime(df_in3["time"], format="%Y%m%d%H%M")  # Datetime

        df_in3["Prec"] = pd.to_numeric(
            df_in3["rre150z0"], errors="coerce"
        )  # Add Precipitation data

        df_in3["vp_a"] = pd.to_numeric(
            df_in3["pva200s0"], errors="coerce"
        )  # Vapour pressure over air

        df_in3["Prec"] = df_in3["Prec"] / 1000



        df_in3["Temperature"] = pd.to_numeric(
            df_in3["tre200s0"], errors="coerce"
        )  # Air temperature

        df_in3["Wind Speed"] = pd.to_numeric(
            df_in3["fkl010z0"], errors="coerce"
        )  # Wind speed

        df_in3["Humidity"] = pd.to_numeric(
            df_in3["ure200s0"], errors="coerce"
        )

        df_in2 = df_in2.set_index("When").resample("5T").ffill().reset_index()
        df_in3 = df_in3.set_index("When").resample("5T").ffill().reset_index()

        mask = (df_in["When"] >= dates["start_date"]) & (df_in["When"] <= dates["error_date"])
        df_in = df_in.loc[mask]
        df_in = df_in.reset_index()

        mask = (df_in2["When"] >= dates["start_date"]) & (
                df_in2["When"] <= dates["error_date"]
        )
        df_2 = df_in2.loc[mask]
        df_2 = df_2.reset_index()

        mask = (df_in3["When"] >= dates["start_date"]) & (
                df_in3["When"] <= dates["error_date"]
        )
        df_3 = df_in3.loc[mask]
        df_3 = df_3.reset_index()

        mask = (df_in3["When"] >= dates["error_date"]) & (
                df_in3["When"] <= dates["end_date"]
        )
        df_4 = df_in3.loc[mask]
        df_4 = df_4.reset_index()

        mask = (df_in2["When"] >= dates["error_date"]) & (
                df_in2["When"] <= dates["end_date"]
        )
        df_5 = df_in2.loc[mask]
        df_5 = df_5.reset_index()
        df_4["Pressure"] = df_in["Pressure"].mean()
        df_4["gre000z0"] = df_5["gre000z0"]
        df_4["oli000z0"] = df_5["oli000z0"]


        df_4["Discharge"] = 0
        mask = df_4["Temperature"] < fountain["crit_temp"]
        mask_index = df_4[mask].index
        df_4.loc[mask_index, "Discharge"] = 15
        mask = df_4["When"] >= dates["fountain_off_date"]
        mask_index = df_4[mask].index
        df_4.loc[mask_index, "Discharge"] = 0

        days = pd.date_range(start=dates["start_date"], end=dates["error_date"], freq="5T")
        days = pd.DataFrame({"When": days})

        df = pd.merge(
            days,
            df_in[
                [
                    "When",
                    "Discharge",
                    "Temperature",
                    "Humidity",
                    "Pressure",
                ]
            ],
            on="When",
        )

        # Fill with other station-
        df["gre000z0"] = df_2["gre000z0"]
        df["oli000z0"] = df_2["oli000z0"]
        df["Prec"] = df_3["Prec"]
        df["vp_a"] = df_3["vp_a"]
        df["Wind Speed"] = df_3["Wind Speed"]


        mask = (df["When"] >= dates["start_date"]) & (df["When"] <= dates["error_date"])
        df = df.loc[mask]
        df = df.reset_index()

        # Add rest of data
        df_4 = df_4.reindex(df_4.index.drop(0)).reset_index(drop=True)
        df = df.append(df_4, ignore_index=True)

        # Error Correction
        df = df.fillna(method='ffill')

        cld = 0.5
        df["Rad"] = df_in2["gre000z0"] - df_in2["gre000z0"] * cld
        df["DRad"] = df_in2["gre000z0"] * cld
        df["cld"] = cld

        # CSV output
        df.rename(
            columns={
                "Wind Speed": "v_a",
                "Temperature": "T_a",
                "Humidity": "RH",
                "Pressure": "p_a",
                "Rad": 'SW_direct',
                "DRad": 'SW_diffuse',
                "oli000z0": 'LW_in',
            },
            inplace=True,
        )

        # Inactive spray
        mask = df["Discharge"] < 8
        mask_index = df[mask].index
        df.loc[mask_index, "Discharge"] = 0

        df_out = df[
            ["When", "T_a", "RH", "v_a", "Discharge", "SW_direct", "SW_diffuse", "LW_in", "cld", "Prec", "p_a", "vp_a"]
        ]

        df_out = df_out.round(5)

        print(df_out.tail())

        df_out.to_csv(folders["input_folder"] + "raw_input.csv")

        fig, ax = plt.subplots()
        ax.plot(df.When, df.Discharge)
        plt.show()
