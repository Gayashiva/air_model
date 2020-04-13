import os
import glob
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
import math
import time
from src.data.config import site, dates, option, folders, fountain, surface
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

start = time.time()

if __name__ == '__main__':

    path = folders["data"]  # use your path
    all_filesA = glob.glob(
        os.path.join(path, "TOA5__Flux_CSFormat_5.dat"))
    all_filesB = glob.glob(
        os.path.join(path, "TOA5__FluxB_CSFormat_5.dat"))

    li = []
    for filename in all_filesA:
        df_in = pd.read_csv(filename, header=1)
        df_in = df_in[2:].reset_index(drop=True)

        li.append(df_in)

    dfa = pd.concat(li, axis=0, ignore_index=True)

    li = []
    for filename in all_filesB:
        df_in = pd.read_csv(filename, header=1)
        df_in = df_in[2:].reset_index(drop=True)

        li.append(df_in)

    dfb = pd.concat(li, axis=0, ignore_index=True)

    dfa["TIMESTAMP"] = pd.to_datetime(dfa["TIMESTAMP"], format='%Y-%m-%d %H:%M:%S.%f')
    dfa["WS"] = pd.to_numeric(dfa["WS"], errors="coerce")

    dfb["TIMESTAMP"] = pd.to_datetime(dfb["TIMESTAMP"], format='%Y-%m-%d %H:%M:%S.%f')
    dfb["WSB"] = pd.to_numeric(dfb["WSB"], errors="coerce")

    dfa = dfa.sort_values(by='TIMESTAMP')
    dfb = dfb.sort_values(by='TIMESTAMP')

    x =dfa.TIMESTAMP
    y1 = dfa["WS"] - dfb["WSB"]
    # y2 = dfb["WS"] + 1

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(x, y1, "k-", linewidth=0.5)
    ax1.set_ylabel("Sonic A")
    ax1.set_xlabel("Days")

    # ax2 = ax1.twinx()
    # ax2.plot(x, y2, "b-", linewidth=0.5)
    # ax2.set_ylabel("Sonic B", color="b")
    # for tl in ax2.get_yticklabels():
    #     tl.set_color("b")
    #
    # ax2.set_ylim(ax1.get_ylim())

    #  format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    plt.show()

    print(dfa.head())


