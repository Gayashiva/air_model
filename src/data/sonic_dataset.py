import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
import math
import time
from tqdm import tqdm
import os
import glob
from src.data.config import site, dates, option, folders, fountain, surface
import fnmatch

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
        y = hs + (t * v * math.sin(theta_f)) - (g * t * t) / 2\
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


dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

start = time.time()

if __name__ == '__main__':

    path = folders["data"]  # use your path
    all_files = glob.glob(
        os.path.join(path, "TOA5__Flux*.dat"))  # advisable to use os.path.join as this makes concatenation OS independent
    pattern = "TOA5__FluxB*.dat"
    li = []
    li_B = []
    for filename in all_files:

        if 'B' in filename:
            df_inB = pd.read_csv(filename, header=1)
            df_inB = df_inB[2:].reset_index(drop=True)
            df_inB = df_inB.drop(["RECORD"], axis=1)
            li_B.append(df_inB)
        else:
            df_in = pd.read_csv(filename, header=1)
            df_in = df_in[2:].reset_index(drop=True)
            li.append(df_in)

    df_A = pd.concat(li, axis=0, ignore_index=True)
    df_B = pd.concat(li_B, axis=0, ignore_index=True)

    for col in df_B.columns:
        if 'B' not in col:
            if col != 'TIMESTAMP':
                df_B = df_B.drop([col], axis=1)


    df_A["TIMESTAMP"] = pd.to_datetime(df_A["TIMESTAMP"], format='%Y-%m-%d %H:%M:%S')
    df_B["TIMESTAMP"] = pd.to_datetime(df_B["TIMESTAMP"], format='%Y-%m-%d %H:%M:%S')

    # mask = (df_A["TIMESTAMP"] >= dates["start_date"]) & (df_A["TIMESTAMP"] <= dates["end_date"])
    # df_A = df_A.loc[mask]
    # df_A = df_A.reset_index()
    # mask = (df_B["TIMESTAMP"] >= dates["start_date"]) & (df_B["TIMESTAMP"] <= dates["end_date"])
    # df_B = df_B.loc[mask]
    # df_B = df_B.reset_index()

    df_A = df_A.sort_values(by='TIMESTAMP')
    df_B = df_B.sort_values(by='TIMESTAMP')

    df = pd.merge(df_A, df_B, how='inner', left_index=True, on='TIMESTAMP')
    print(df.info())

    df["H"] = pd.to_numeric(df["H"], errors="coerce")
    df["HB"] = pd.to_numeric(df["HB"], errors="coerce")
    df["SW_IN"] = pd.to_numeric(df["SW_IN"], errors="coerce")
    df["LW_IN"] = pd.to_numeric(df["LW_IN"], errors="coerce")
    df["amb_press_Avg"] = pd.to_numeric(df["amb_press_Avg"], errors="coerce")
    df["e_probe"] = pd.to_numeric(df["e_probe"], errors="coerce")
    df["NETRAD"] = pd.to_numeric(df["NETRAD"], errors="coerce")
    df["T_probe_Avg"] = pd.to_numeric(df["T_probe_Avg"], errors="coerce")
    df["T_SONIC"] = pd.to_numeric(df["T_SONIC"], errors="coerce")
    df["RH_probe_Avg"] = pd.to_numeric(df["RH_probe_Avg"], errors="coerce")
    df["Waterpressure"] = pd.to_numeric(df["Waterpressure"], errors="coerce")
    df["WaterFlow"] = pd.to_numeric(df["WaterFlow"], errors="coerce")
    df["WS"] = pd.to_numeric(df["WS"], errors="coerce")
    df["WS_MAX"] = pd.to_numeric(df["WS"], errors="coerce")
    df["WSB"] = pd.to_numeric(df["WSB"], errors="coerce")
    df["SnowHeight"] = pd.to_numeric(df["SnowHeight"], errors="coerce")

    # for i in range(1,9):
    #     col = 'Tice_Avg(' + str(i) + ')'
    #     df[col] = pd.to_numeric(df[col], errors="coerce")


    df.to_csv(folders["input_folder"] + "raw_output.csv")

    # Errors
    df['H'] = df['H'] / 1000
    df['HB'] = df['HB'] / 1000
    df['H'] = df['H'].apply(lambda x: x if abs(x) < 500 else np.NAN)
    df['HB'] = df['HB'].apply(lambda x: x if abs(x) < 500 else np.NAN)

    g = 9.81
    h_aws = 1.2
    z0 = 0.0017
    c_a = 1.01 * 1000
    rho_a = 1.29
    p0 = 1013
    k = 0.4
    df['Ri_b'] = 0



    for i in range(0,df.shape[0]):

        df.loc[i,'Ri_b'] = (
                g
                * (h_aws - z0)
                * (df.loc[i, "T_SONIC"])
                / ((df.loc[i, "T_SONIC"]+273) * df.loc[i, "WS"] ** 2))

        # Sensible Heat
        df.loc[i, "HC"] = (
                c_a
                * rho_a
                * df.loc[i, "amb_press_Avg"] * 10
                / p0
                * math.pow(k, 2)
                * df.loc[i, "WS"]
                * (df.loc[i, "T_SONIC"])
                / ((np.log(h_aws / z0)) ** 2)
        )

        if df.loc[i,'Ri_b'] < 0.2:

            if df.loc[i,'Ri_b'] > 0:
                df.loc[i, "HC"] = df.loc[i, "HC"] * (1- 5 * df.loc[i,'Ri_b']) ** 2
            else:
                df.loc[i, "HC"] = df.loc[i, "HC"] * math.pow((1- 16 * df.loc[i,'Ri_b']), 0.75)


    # print(df["amb_press_Avg"].head(), df["T_probe_Avg"].head(), df["WS"].head(), df["HC"].head() )

    dfd = df.set_index("TIMESTAMP").resample("D").mean().reset_index()

    """Input Plots"""

    pp = PdfPages(folders["input_folder"] + site + "_data" + ".pdf")

    x = df.TIMESTAMP

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y1 = df.T_probe_Avg
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Temperature [$\degree C$]")
    ax1.grid()

    ax1t = ax1.twinx()
    ax1t.plot(x, df.T_SONIC, "b-", linewidth=0.5)
    ax1t.set_ylabel("Temperature Sonic [$\degree C$]", color="b")
    for tl in ax1t.get_yticklabels():
        tl.set_color("b")

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)

    y2 = df.WaterFlow
    ax1.plot(x, y2, "k-", linewidth=0.5)
    ax1.set_ylabel("Discharge Rate ")
    ax1.grid()

    ax1t = ax1.twinx()
    ax1t.plot(x, df.Waterpressure, "b-", linewidth=0.5)
    ax1t.set_ylabel("Water Pressure [$Bar$]", color="b")
    for tl in ax1t.get_yticklabels():
        tl.set_color("b")

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)

    y3 = df.NETRAD
    ax1.plot(x, y3, "k-", linewidth=0.5)
    ax1.set_ylabel("Net Radiation [$W\,m^{-2}$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)

    x1 = dfd.TIMESTAMP
    y31 = -dfd.H
    ax1.plot(x1, y31, "k-", linewidth=0.5)
    ax1.set_ylabel("Sensible Heat A [$W\,m^{-2}$]")
    ax1.grid()

    ax1t = ax1.twinx()
    ax1t.plot(x1, dfd.HC, "b-", linewidth=0.5)
    ax1t.set_ylabel("Sensible Heat C [$W\,m^{-2}$]", color="b")
    for tl in ax1t.get_yticklabels():
        tl.set_color("b")

    ax1.set_ylim([-100,100])
    ax1t.set_ylim([-100,100])
    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)

    y31 = -dfd.HB
    ax1.plot(x1, y31, "k-", linewidth=0.5)
    ax1.set_ylabel("Sensible Heat B [$W\,m^{-2}$]")
    ax1.grid()

    ax1.set_ylim([-100,100])

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)
    y31 = dfd.HC
    ax1.plot(x1, y31, "k-", linewidth=0.5)
    ax1.set_ylabel("Sensible Heat C [$W\,m^{-2}$]")
    ax1.grid()

    # ax1.set_ylim([-100,100])

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)

    y4 = df.SnowHeight
    ax1.plot(x, y4, "k-", linewidth=0.5)
    ax1.set_ylabel("Snow Height [$cm$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)

    y5 = df.amb_press_Avg
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

    ax1 = fig.add_subplot(111)

    y6 = df.WS
    ax1.plot(x, y6, "k-", linewidth=0.5)
    ax1.set_ylabel("Wind A [$m\,s^{-1}$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)
    y7 = df.WSB
    ax1.plot(x, y7, "k-", linewidth=0.5)
    ax1.set_ylabel("Wind B [$m\,s^{-1}$]")
    ax1.grid()

    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    ax1 = fig.add_subplot(111)
    y7 = df.Ri_b
    ax1.plot(x, y7, "k-", linewidth=0.5)
    ax1.set_ylabel("Ri")
    ax1.grid()
    # ax1.set_ylim([0,0.2])


    # format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    pp.close()

    # # CSV output
    # df.rename(
    #     columns={
    #         "TIMESTAMP": "When",
    #         "LW_IN": "oli000z0",
    #         "SW_IN": "Rad",
    #         "amb_press_Avg": "p_a",
    #         "WSB": "v_a",
    #         "e_probe": "vp_a",
    #         "WaterFlow": "Discharge",
    #         "T_probe_Avg": "T_a",
    #         "RH_probe_Avg": "RH",
    #         "HB": "H_eddy",
    #     },
    #     inplace=True,
    # )
    #
    # for i in tqdm(range(1, df.shape[0])):
    #
    #     if np.isnan(df.loc[i, "Discharge"]) or np.isinf(df.loc[i, "Discharge"]):
    #         df.loc[i, "Discharge"] = df.loc[i - 1, "Discharge"]
    #     if np.isnan(df.loc[i, "H_eddy"]) or np.isinf(df.loc[i, "H_eddy"]):
    #         df.loc[i, "H_eddy"] = df.loc[i - 1, "H_eddy"]
    #
    #     """Solar Elevation Angle"""
    #     df.loc[i, "SEA"] = getSEA(
    #         df.loc[i, "When"],
    #         fountain["latitude"],
    #         fountain["longitude"],
    #         fountain["utc_offset"],
    #     )
    #
    # df = df.set_index("When").resample("30T").ffill().reset_index()
    #
    # df_out = df[
    #     ["When", "T_a", "RH", "v_a", "Discharge", "Rad", "DRad", "Prec", "p_a", "SEA", "vp_a", "H_eddy"]
    # ]
    #
    # print(df_out.head())
    #
    # df_out = df_out.round(5)
    #
    #
    # filename = folders["input_folder"] + site
    #
    # df_out.to_csv(filename + "_raw_input.csv")
    #
    # """ Derived Parameters"""
    #
    # l = [
    #     "a",
    #     "r_f",
    # ]
    # for col in l:
    #     df[col] = 0
    #
    # """Albedo Decay"""
    # surface["decay_t"] = (
    #         surface["decay_t"] * 24 * 60 / 30
    # )  # convert to 30 minute time steps
    # s = 0
    # f = 0
    #
    # for i in tqdm(range(1, df.shape[0])):
    #
    #     if option == "schwarzsee":
    #
    #         ti = surface["decay_t"]
    #         a_min = surface["a_i"]
    #
    #         # Precipitation
    #         if (df.loc[i, "Discharge"] == 0) & (df.loc[i, "Prec"] > 0):
    #             if df.loc[i, "T_a"] < surface["rain_temp"]:  # Snow
    #                 s = 0
    #                 f = 0
    #
    #         if df.loc[i, "Discharge"] > 0:
    #             f = 1
    #             s = 0
    #
    #         if f == 0:  # last snowed
    #             df.loc[i, "a"] = a_min + (surface["a_s"] - a_min) * math.exp(-s / ti)
    #             s = s + 1
    #         else:  # last sprayed
    #             df.loc[i, "a"] = a_min
    #             s = s + 1
    #     else:
    #         df.loc[i, "a"] = surface["a_i"]
    #
    #     """ Fountain Spray radius """
    #     if df.loc[i, 'When'] > datetime(2020, 1, 9):
    #         fountain["aperture_f"] = 0.003
    #     else:
    #         fountain["aperture_f"] = 0.005
    #
    #     Area = math.pi * math.pow(fountain["aperture_f"], 2) / 4
    #     v_f = df.loc[i, "Discharge"] / (60 * 1000 * Area)
    #     df.loc[i, "r_f"] = projectile_xy(
    #         v_f, fountain["h_f"]
    #     )
    #
    # df.to_csv(filename + "_input.csv")

    # pp = PdfPages(folders["input_folder"] + site + "_derived_parameters" + ".pdf")
    #
    # x = df.When
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    #
    # y1 = df.Discharge
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Discharge [$l\, min^{-1}$]")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    #
    # y1 = df.r_f
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Spray Radius [$m$]")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # y1 = df.vp_a
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Vapour Pressure")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # y1 = df.a
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Albedo")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # pp.close()

    # pp = PdfPages(folders["input_folder"] + "schwarzsee_2020.pdf")
    #
    # x = df.When
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    #
    # y1 = df.WaterFlow
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Discharge [$l\, min^{-1}$]")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    #
    # y1 = df.T_probe_Avg
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Temperature")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    #
    # y1 = df.RH_probe_Avg
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Humidity")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # y1 = df.WS
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Wind speed")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # y1 = df.SnowHeight
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("SnowHeight")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # # y1 = df_in["Tice"]
    # for i in range(1,3):
    #     col = 'Tice_Avg(' + str(i) + ')'
    #     plt.plot(x, df[col], label='id %s' % i)
    # plt.legend()
    #
    # ax1.set_ylabel("Ice Temperatures")
    # ax1.set_ylim([-1,3])
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # y1 = df.NETRAD
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("NETRAD")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # y1 = df.H
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Sensible heat")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # y1 = df['H'] + df['NETRAD']
    # ax1.plot(x, y1, "k-", linewidth=0.5)
    # ax1.set_ylabel("Sensible heat")
    # ax1.grid()
    #
    # # format the ticks
    # ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    # ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.grid()
    # fig.autofmt_xdate()
    # pp.savefig(bbox_inches="tight")
    # plt.clf()
    #
    # pp.close()