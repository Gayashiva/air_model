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

dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

start = time.time()

if __name__ == '__main__':

    path = folders["data"]  # use your path
    all_files = glob.glob(
        os.path.join(path, "TOA5__Flux_CSFormat_*.dat"))  # advisable to use os.path.join as this makes concatenation OS independent

    li = []
    for filename in all_files:
        df = pd.read_csv(filename, header=1)
        df = df[2:].reset_index(drop=True)
        # df["TIMESTAMP"] = pd.to_datetime(df["TIMESTAMP"], format='%Y-%m-%d %H:%M:%S')
        df = df.set_index('RECORD')

        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=False)
    print(frame.info())
    print(frame.tail())
    frame = frame.reset_index()

    df_in = frame

    df_in["TIMESTAMP"] = pd.to_datetime(df_in["TIMESTAMP"], format='%Y-%m-%d %H:%M:%S')
    df_in["H"] = pd.to_numeric(df_in["H"], errors="coerce")
    df_in["NETRAD"] = pd.to_numeric(df_in["NETRAD"], errors="coerce")
    df_in["T_probe_Avg"] = pd.to_numeric(df_in["T_probe_Avg"], errors="coerce")
    df_in["RH_probe_Avg"] = pd.to_numeric(df_in["RH_probe_Avg"], errors="coerce")
    df_in["Waterpressure"] = pd.to_numeric(df_in["Waterpressure"], errors="coerce")
    df_in["WaterFlow"] = pd.to_numeric(df_in["WaterFlow"], errors="coerce")
    df_in["Tice_Avg(1)"] = pd.to_numeric(df_in["Tice_Avg(1)"], errors="coerce")
    df_in["WS"] = pd.to_numeric(df_in["WS"], errors="coerce")
    df_in["SnowHeight"] = pd.to_numeric(df_in["SnowHeight"], errors="coerce")

    # df_in = df_in.set_index('TIMESTAMP')

    # Errors
    df_in['H'] = df_in['H']/1000

    print(df_in['RECORD'])

    pp = PdfPages(folders["input_folder"] + "schwarzsee_2020.pdf")

    x = df_in.TIMESTAMP

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y1 = df_in.WaterFlow
    y2 = df_in.Waterpressure
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Discharge [$l\, min^{-1}$]")
    ax1.grid()

    ax1t = ax1.twinx()
    ax1t.plot(x, y2, "b-", linewidth=0.5)
    ax1t.set_ylabel("Pressure", color="b")
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

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    y1 = df_in.T_probe_Avg
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Temperature")
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

    y1 = df_in.RH_probe_Avg
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Humidity")
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
    y1 = df_in.WS
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Wind speed")
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
    y1 = df_in.SnowHeight
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("SnowHeight")
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
    y1 = df_in["Tice_Avg(1)"]
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Ice Temperature")
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
    y1 = df_in.NETRAD
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("NETRAD")
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
    y1 = df_in.H
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Sensible heat")
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