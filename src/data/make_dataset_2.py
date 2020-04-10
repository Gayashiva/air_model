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
import glob
import logging
from src.data.config import site, dates, option, folders, fountain, surface

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


dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

start = time.time()

if __name__ == '__main__':

    path = folders["data"]  # use your path
    all_files = glob.glob(
        os.path.join(path, "TOA5__Flux_CSFormat_*.dat"))  # advisable to use os.path.join as this makes concatenation OS independent

    li = []
    for filename in all_files:
        df_in = pd.read_csv(filename, header=1)
        df_in = df_in[2:].reset_index(drop=True)

        li.append(df_in)

    df = pd.concat(li, axis=0, ignore_index=True)

    df["TIMESTAMP"] = pd.to_datetime(df["TIMESTAMP"], format='%Y-%m-%d %H:%M:%S')
    df["H"] = pd.to_numeric(df["H"], errors="coerce")
    df["SW_IN"] = pd.to_numeric(df["SW_IN"], errors="coerce")
    df["LW_IN"] = pd.to_numeric(df["LW_IN"], errors="coerce")
    df["amb_press_Avg"] = pd.to_numeric(df["amb_press_Avg"], errors="coerce")
    df["e_probe"] = pd.to_numeric(df["e_probe"], errors="coerce")
    df["NETRAD"] = pd.to_numeric(df["NETRAD"], errors="coerce")
    df["T_probe_Avg"] = pd.to_numeric(df["T_probe_Avg"], errors="coerce")
    df["RH_probe_Avg"] = pd.to_numeric(df["RH_probe_Avg"], errors="coerce")
    df["Waterpressure"] = pd.to_numeric(df["Waterpressure"], errors="coerce")
    df["WaterFlow"] = pd.to_numeric(df["WaterFlow"], errors="coerce")
    df["WS"] = pd.to_numeric(df["WS"], errors="coerce")
    df["SnowHeight"] = pd.to_numeric(df["SnowHeight"], errors="coerce")

    for i in range(1,9):
        col = 'Tice_Avg(' + str(i) + ')'
        df[col] = pd.to_numeric(df[col], errors="coerce")

    mask = (df["TIMESTAMP"] >= dates["start_date"]) & (df["TIMESTAMP"] <= dates["end_date"])
    df = df.loc[mask]
    df = df.reset_index()

    # Errors
    df['H'] = df['H'] / 1000
    df['DRad'] = 0
    df['Prec'] = 0
    df['e_probe'] = df['e_probe'] * 10
    df['H'] = df['H'].apply(lambda x: x if abs(x) < 500 else np.NAN)

    df = df.sort_values(by='TIMESTAMP')

    # CSV output
    df.rename(
        columns={
            "TIMESTAMP": "When",
            "LW_IN": "oli000z0",
            "SW_IN": "Rad",
            "amb_press_Avg": "p_a",
            "WS": "v_a",
            "e_probe": "vp_a",
            "WaterFlow": "Discharge",
            "T_probe_Avg": "T_a",
            "RH_probe_Avg": "RH",
        },
        inplace=True,
    )

    for i in tqdm(range(1, df.shape[0])):

        if np.isnan(df.loc[i, "Discharge"]) or np.isinf(df.loc[i, "Discharge"]):
            df.loc[i, "Discharge"] = df.loc[i - 1, "Discharge"]

        """Solar Elevation Angle"""
        df.loc[i, "SEA"] = getSEA(
            df.loc[i, "When"],
            fountain["latitude"],
            fountain["longitude"],
            fountain["utc_offset"],
        )

    df_out = df[
        ["When", "T_a", "RH", "v_a", "Discharge", "Rad", "DRad", "Prec", "p_a", "SEA", "vp_a"]
    ]

    df_out = df_out.round(5)


    filename = folders["input_folder"] + site

    df_out.to_csv(filename + "_raw_input.csv")

    """ Derived Parameters"""

    l = [
        "a",
        "r_f",
    ]
    for col in l:
        df[col] = 0

    """Albedo Decay"""
    surface["decay_t"] = (
            surface["decay_t"] * 24 * 60 / 5
    )  # convert to 5 minute time steps
    s = 0
    f = 0

    for i in tqdm(range(1, df.shape[0])):

        if option == "schwarzsee":

            ti = surface["decay_t"]
            a_min = surface["a_i"]

            # Precipitation
            if (df.loc[i, "Discharge"] == 0) & (df.loc[i, "Prec"] > 0):
                if df.loc[i, "T_a"] < surface["rain_temp"]:  # Snow
                    s = 0
                    f = 0

            if df.loc[i, "Discharge"] > 0:
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
        if df.loc[i, 'When'] > datetime(2020, 1, 9):
            fountain["aperture_f"] = 0.003
        else:
            fountain["aperture_f"] = 0.005

        Area = math.pi * math.pow(fountain["aperture_f"], 2) / 4
        v_f = df.loc[i, "Discharge"] / (60 * 1000 * Area)
        df.loc[i, "r_f"] = projectile_xy(
            v_f, fountain["h_f"]
        )

    df.to_csv(filename + "_input.csv")

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

    x = df.When

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