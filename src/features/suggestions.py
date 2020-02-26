import logging
import os
import os.path
import time
from datetime import datetime
from logging import StreamHandler
from pandas.plotting import register_matplotlib_converters

register_matplotlib_converters()
import matplotlib
from matplotlib.offsetbox import AnchoredText
import math
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from src.data.config import site, dates, option, folders, fountain, surface
from src.models.air_forecast import icestupa
from src.data.make_dataset import projectile_xy, discharge_rate

plt.rcParams["figure.figsize"] = (10,7)

start = time.time()

# Create the Logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create the Handler for logging data to a file
logger_handler = logging.FileHandler(
    os.path.join(os.path.join(folders["dirname"], "data/logs/"), site + "_site.log"),
    mode="w",
)
logger_handler.setLevel(logging.DEBUG)

# Create the Handler for logging data to console.
console_handler = StreamHandler()
console_handler.setLevel(logging.CRITICAL)

# Create a Formatter for formatting the log messages
logger_formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

# Add the Formatter to the Handler
logger_handler.setFormatter(logger_formatter)
console_handler.setFormatter(logger_formatter)

# Add the Handler to the Logger
logger.addHandler(logger_handler)
logger.addHandler(console_handler)


fountain['discharge'] = 8

#  read files
filename0 = os.path.join(folders["input_folder"] + site + "_raw_input.csv")
df_in = pd.read_csv(filename0, sep=",")
df_in["When"] = pd.to_datetime(df_in["When"], format="%Y.%m.%d %H:%M:%S")

""" Derived Parameters"""

l = [
    "a",
    "r_f",
    "Fountain",
]
for col in l:
    df_in[col] = 0

"""Discharge Rate"""
df_in["Fountain"], df_in["Discharge"] = discharge_rate(df_in, fountain)
df_in["Discharge"] = fountain["discharge"] * df_in["Fountain"]

"""Albedo Decay"""
surface["decay_t"] = (
        surface["decay_t"] * 24 * 60 / 5
)  # convert to 5 minute time steps
s = 0
f = 0

""" Fountain Spray radius """
Area = math.pi * math.pow(fountain["aperture_f"], 2) / 4

for i in range(1, df_in.shape[0]):

    if option == "schwarzsee":

        ti = surface["decay_t"]
        a_min = surface["a_i"]

        # Precipitation
        if (df_in.loc[i, "Fountain"] == 0) & (df_in.loc[i, "Prec"] > 0):
            if df_in.loc[i, "T_a"] < surface["rain_temp"]:  # Snow
                s = 0
                f = 0

        if df_in.loc[i, "Fountain"] > 0:
            f = 1
            s = 0

        if f == 0:  # last snowed
            df_in.loc[i, "a"] = a_min + (surface["a_s"] - a_min) * math.exp(-s / ti)
            s = s + 1
        else:  # last sprayed
            df_in.loc[i, "a"] = a_min
            s = s + 1
    else:
        df_in.loc[i, "a"] = surface["a_i"]

    """ Fountain Spray radius """
    v_f = df_in.loc[i, "Discharge"] / (60 * 1000 * Area)
    df_in.loc[i, "r_f"] = projectile_xy(
        v_f, fountain["h_f"]
    )


df = icestupa(df_in, fountain, surface)

# Full Output
if option == "temperature":
    filename2 = (
        folders["output_folder"] + site + "_" + option + "_" + str(fountain["crit_temp"])
    )
else:
    filename2 = folders["output_folder"] + site + "_" + option
filename4 = os.path.join(filename2 + "_model_results.csv")
df.to_csv(filename4, sep=",")

df = df.rename({'SW': '$SW_{net}$', 'LW': '$LW_{net}$', 'Qs': '$Q_S$', 'Ql': '$Q_L$', 'Qc': '$Q_C$' }, axis=1)

# Plots

filename2 = filename2 + "_" + str(fountain['discharge'])

pp = PdfPages(filename2 + "_results.pdf")

x = df.When
y1 = df.iceV

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "k-")
ax1.set_ylabel("Ice Volume [$m^3$]")
ax1.set_xlabel("Days")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

y1 = df.SA

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "k-")
ax1.set_ylabel("Surface Area [$m^2$]")
ax1.set_xlabel("Days")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

y1 = df.h_ice
y2 = df.r_ice

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "k-")
ax1.set_ylabel("Ice Cone Height [$m$]")
ax1.set_xlabel("Days")

ax2 = ax1.twinx()
ax2.plot(x, y2, "b-", linewidth=0.5)
ax2.set_ylabel("Ice Radius", color="b")
for tl in ax2.get_yticklabels():
    tl.set_color("b")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

y1 = df.SRf

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "k-", linewidth=0.5)
ax1.set_ylabel("Solar Area fraction")
ax1.set_xlabel("Days")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

y1 = df.iceV
y2 = df['TotalE'] + df['$Q_L$']

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "k-")
ax1.set_ylabel("Ice Volume [$m^3$]")
ax1.set_xlabel("Days")

ax2 = ax1.twinx()
ax2.plot(x, y2, "b-", linewidth=0.5)
ax2.set_ylabel("Energy [$W\,m^{-2}$]", color="b")
for tl in ax2.get_yticklabels():
    tl.set_color("b")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

y1 = df.T_s

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "k-", linewidth=0.5)
ax1.set_ylabel("Surface Temperature [$\degree C$]")
ax1.set_xlabel("Days")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

y1 = df.solid / 5

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "b-", linewidth=0.5)
ax1.set_ylabel("Ice Production rate [$l\,min^{-1}$]")
ax1.set_xlabel("Days")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

y1 = df.gas / 5
y2 = df.deposition / 5

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "k-", linewidth=0.5)
ax1.set_ylabel("Gas Production rate [$l\,min^{-1}$]")
ax1.set_xlabel("Days")

ax2 = ax1.twinx()
ax2.plot(x, y2, "b-", linewidth=0.5)
ax2.set_ylabel("Deposition rate [$l\,min^{-1}$]", color="b")
for tl in ax2.get_yticklabels():
    tl.set_color("b")

ax2.set_ylim(ax1.get_ylim())

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()
plt.close('all')

pp.close()

# Model Suggestion
if fountain["discharge"] == 10.5:
    filename = os.path.join(filename2 + "_suggestion.pdf")
    pp = PdfPages(filename)

    x = df.When
    y1 = df.iceV

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(x, y1, "b-", lw=1)
    ax1.set_ylabel("Ice Volume [$m^3$]")
    ax1.set_xlabel("Days")

    # Include Validation line segment 1
    ax1.plot(
        [datetime(2019, 2, 14, 16), datetime(2019, 2, 14, 16)],
        [0.67115, 1.042],
        color="green",
        lw=1,
    )
    ax1.scatter(datetime(2019, 2, 14, 16), 0.856575, color="green", marker="o")

    # Include Validation line segment 2
    ax1.plot(
        [datetime(2019, 3, 10, 18), datetime(2019, 3, 10, 18)],
        [0.037, 0.222],
        color="green",
        lw=1,
    )
    ax1.scatter(datetime(2019, 3, 10, 18), 0.1295, color="green", marker="o")

    #  format the ticks
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.grid()
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")

    plt.savefig(
        os.path.join(folders["output_folder"], site + "_result.jpg"),
        bbox_inches="tight",
        dpi=300,
    )

    plt.clf()

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(10, 5))

    # fig.suptitle("Mass and Energy balance", fontsize=14)

    x = df.When
    y1 = df.ice

    ax1.plot(x, y1, "b-")
    ax1.set_ylabel("Ice [kg]")
    ax1.grid()
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())

    y2 = df.SA
    ax2.plot(x, y2, "b-")
    ax2.set_ylabel("Surface Area [$m^2$]")
    ax2.grid()

    y3 = df['TotalE'] + df['$Q_L$']
    ax3.plot(x, y3, "b-")
    ax3.set_ylabel("E [$W\,m^{-2}$]")
    ax3.grid()

    # rotates and right aligns the x labels, and moves the bottom of the axes up to make room for them
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    fig, (ax1, ax2, ax3) = plt.subplots(
        nrows=3, ncols=1, sharex=True, sharey=True, figsize=(10, 5)
    )

    # fig.suptitle("Mass balance", fontsize=14)

    x = df.When
    y1 = df.ice

    ax1.plot(x, y1, "b-")
    ax1.set_ylabel("Ice [kg]")
    ax1.grid()
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_minor_locator(mdates.DayLocator())

    y2 = df.meltwater
    ax2.plot(x, y2, "b-")
    ax2.set_ylabel("Meltwater [kg]")
    ax2.grid()

    y3 = df.vapour
    ax3.plot(x, y3, "b-")
    ax3.set_ylabel("Vapour [kg]")
    ax3.grid()
    # rotates and right aligns the x labels, and moves the bottom of the axes up to make room for them
    fig.autofmt_xdate()
    pp.savefig(bbox_inches="tight")
    plt.clf()

    pp.close()

"""Paper Output"""
filename = os.path.join(filename2 + "_paper.pdf")
pp = PdfPages(filename)

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


x = df.When
y1 = df.iceV

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "b-", lw=1)
ax1.set_ylabel("Ice Volume [$m^3$]")
ax1.set_xlabel("Days")

# Include Validation line segment 1
ax1.plot(
    [datetime(2019, 2, 14, 16), datetime(2019, 2, 14, 16)],
    [0.67115, 1.042],
    color="green",
    lw=1,
)
ax1.scatter(datetime(2019, 2, 14, 16), 0.856575, color="green", marker="o")

# Include Validation line segment 2
ax1.plot(
    [datetime(2019, 3, 10, 18), datetime(2019, 3, 10, 18)],
    [0.037, 0.222],
    color="green",
    lw=1,
)
ax1.scatter(datetime(2019, 3, 10, 18), 0.1295, color="green", marker="o")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")

plt.savefig(
    os.path.join(folders["output_folder"], site + "_result.jpg"),
    bbox_inches="tight",
    dpi=300,
)

plt.clf()

# Day melt and Night freeze Plots

# df = df.reset_index()
for i in df.index:
    if df.loc[i, 'solid'] < 0:
        df.loc[i, 'solid'] = 0

dfds = df.set_index("When").resample("D").sum().reset_index()

dfd = df.set_index("When").resample("D").mean().reset_index()
dfd['When'] = dfd['When'].dt.strftime("%b %d")
dfd["Discharge"] = dfd["Discharge"] == 0
dfd["Discharge"] = dfd["Discharge"].astype(int)
dfd["Discharge"] = dfd["Discharge"].astype(str)

dfds['melted'] = dfds['melted'] * -1 / (dfd['SA'] * 1000) * 100
dfds['solid'] = dfds['solid'] / (dfd['SA'] * 916) * 100
dfds["When"] = pd.to_datetime(dfds["When"], format="%Y.%m.%d %H:%M:%S")
dfds['When'] = dfds['When'].dt.strftime("%b %d")

dfds2 = dfds[['When','solid', 'melted']]
dfds2 = dfds2.rename({'solid': 'Ice frozen', 'melted': 'Meltwater discharged'}, axis=1)

dfds2["label"] = ' '
labels = ["Jan 29", "Feb 05", "Feb 12", "Feb 19", "Feb 26", "Mar 05", "Mar 12", "Mar 19"]
for i in range(0, dfds2.shape[0]):
    for item in labels:
        if dfds2.When[i] == item:
            dfds2.loc[i, 'label'] = dfds2.When[i]

dfds2 = dfds2.set_index("label")

dfd["label"] = ' '
labels = ["Jan 29", "Feb 05", "Feb 12", "Feb 19", "Feb 26", "Mar 05", "Mar 12", "Mar 19"]
for i in range(0, dfd.shape[0]):
    for item in labels:
        if dfd.When[i] == item:
            dfd.loc[i, 'label'] = dfd.When[i]

dfd = dfd.set_index("label")

fig, (ax1, ax2, ax3) = plt.subplots(
    nrows=3, ncols=1, sharex=True, figsize=(15, 12)
)
fig.subplots_adjust(hspace=0)

y1 = dfds2[['Ice frozen', 'Meltwater discharged']]
z = dfd[['$SW_{net}$', '$LW_{net}$', '$Q_S$', '$Q_L$', '$Q_C$']]
y3 = dfd['SA']


y1.plot(kind='bar', stacked=True, edgecolor='black', linewidth=0.5, color=['#D9E9FA', '#0C70DE'], ax=ax1)
ax1.set_ylabel('Thickness [$cm$]')
ax1.set_ylim(-2.5, 2.5)
ax1.legend(loc='upper right' , prop={'size': 6})

ax1.grid( axis="y",color="black", alpha=.3, linewidth=.5, which="major")
at = AnchoredText("(a)",
                  prop=dict(size=6), frameon=True,
                  loc='upper left',
                  )
at.patch.set_boxstyle("round,pad=0,rounding_size=0.2")
ax1.add_artist(at)


z.plot.bar(stacked=True, edgecolor=dfd['Discharge'], linewidth=0.5, ax=ax2)
ax2.set_ylabel('Energy [$W\,m^{-2}$]')
ax2.legend(loc='lower right', prop={'size': 6})
ax2.set_ylim(-120, 120)
ax2.grid(axis="y", color="black", alpha=.3, linewidth=.5, which="major")
at = AnchoredText("(b)",
                  prop=dict(size=6), frameon=True,
                  loc='upper left',
                  )
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax2.add_artist(at)

y3.plot.bar(edgecolor='#0C70DE',  color=['#D9E9FA'], linewidth=0.5, ax=ax3)
ax3.set_ylabel('Area [$m^2$]')
ax3.grid(axis="y", color="black", alpha=.3, linewidth=.5, which="major")
at = AnchoredText("(c)",
                  prop=dict(size=6), frameon=True,
                  loc='upper left',
                  )
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax3.add_artist(at)

plt.xlabel('Days')
plt.xticks(rotation=45)
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

z.plot.bar(stacked=True, edgecolor=dfd['Discharge'], linewidth=0.5)
plt.xlabel('Days')
plt.ylabel('Energy [$W\,m^{-2}$]')
plt.legend(loc='upper left')
plt.ylim(-150, 150)
plt.xticks(rotation=45)
pp.savefig(bbox_inches="tight")
plt.clf()
plt.close('all')

fig, (ax1, ax2, ax3) = plt.subplots(
    nrows=3, ncols=1, sharex=True, sharey=True, figsize=(10, 5)
)

# fig.suptitle("Daily ice and melt", fontsize=14)

x = df.When
y1 = df.ice

ax1.plot(x, y1, "b-")
ax1.set_ylabel("Ice [kg]")
ax1.grid()
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())

y2 = df.meltwater
ax2.plot(x, y2, "b-")
ax2.set_ylabel("Meltwater [kg]")
ax2.grid()

y3 = df.vapour
ax3.plot(x, y3, "b-")
ax3.set_ylabel("Vapour [kg]")
ax3.grid()
# rotates and right aligns the x labels, and moves the bottom of the axes up to make room for them
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()


pp.close()
