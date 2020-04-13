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

#  read files
if option == "temperature":
    filename1 = (
        folders["output_folder"] + site + "_" + option + "_" + str(fountain["crit_temp"])
    )
else:
    filename1 = folders["output_folder"] + site + "_" + option

filename1 = os.path.join(filename1 + "_model_results.csv")

same = False

if same:
    if os.path.isfile(filename1):
        print("Simulation Exists")
        df = pd.read_csv(filename1, sep=",")
        df["When"] = pd.to_datetime(df["When"], format="%Y.%m.%d %H:%M:%S")

    else:
        filename0 = os.path.join(folders["input_folder"] + site + "_input.csv")
        df_in = pd.read_csv(filename0, sep=",")
        df_in["When"] = pd.to_datetime(df_in["When"], format="%Y.%m.%d %H:%M:%S")

        df = icestupa(df_in, fountain, surface)

        total = time.time() - start

        print("Total time : ", total / 60)

else:
    filename0 = os.path.join(folders["input_folder"] + site + "_input.csv")
    df_in = pd.read_csv(filename0, sep=",")
    df_in["When"] = pd.to_datetime(df_in["When"], format="%Y.%m.%d %H:%M:%S")

    df = icestupa(df_in, fountain, surface)

    total = time.time() - start

    print("Total time : ", total / 60)


# # Output for manim
# filename2 = os.path.join(folders["output_folder"], site + "_model_gif.csv")
# cols = ["When", "h_ice", "h_f", "r_ice", "ice", "T_a", "Discharge"]
# df[cols].to_csv(filename2, sep=",")

# Output for energy balance
# filename3 = os.path.join(folders["output_folder"], site + "_model_energy.csv")
# cols = ["When", "SW", "LW", "Qs", "Ql", "SA", "iceV"]
# df[cols].to_csv(filename3, sep=",")

# Full Output
if option == "temperature":
    filename2 = (
        folders["output_folder"] + site + "_" + option + "_" + str(fountain["crit_temp"])
    )
else:
    filename2 = folders["output_folder"] + site + "_" + option
filename4 = os.path.join(filename2 + "_model_results.csv")
df.to_csv(filename4, sep=",")

df['melt_thick'] = df['melted']/ (df['SA'] * 1000)

df = df.rename({'SW': '$SW_{net}$', 'LW': '$LW_{net}$', 'Qs': '$Q_S$', 'Ql': '$Q_L$', 'Qc': '$Q_C$' }, axis=1)

print(df['T_s'].tail())
print(df['$Q_S$'].tail())
# Plots


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

y1 = df.thickness * 1000

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "b-", linewidth=0.5)
ax1.set_ylabel("Thickness melted [$mm$]")
ax1.set_xlabel("Days")

#  format the ticks
ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.xaxis.set_minor_locator(mdates.DayLocator())
ax1.grid()
fig.autofmt_xdate()
pp.savefig(bbox_inches="tight")
plt.clf()

y1 = df['$Q_S$']
y2 = df.H_eddy

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, "k-", linewidth=0.5)
ax1.set_ylabel("Sensible Heat Derived")
ax1.set_xlabel("Days")

ax2 = ax1.twinx()
ax2.plot(x, y2, "b-", linewidth=0.5)
ax2.set_ylabel("Sensible Heat Original", color="b")
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

# Day melt and Night freeze Plots

for i in df.index:
    if df.loc[i, 'solid'] < 0:
        df.loc[i, 'solid'] = 0

dfds = df.set_index("When").resample("D").sum().reset_index()

dfd = df.set_index("When").resample("D").mean().reset_index()
dfd['When'] = dfd['When'].dt.strftime("%b %d")
dfd["Discharge"] = dfd["Discharge"] == 0
dfd["Discharge"] = dfd["Discharge"].astype(int)
dfd["Discharge"] = dfd["Discharge"].astype(str)

dfds['melted'] = dfds['melted'] * -1 / (dfd['SA'] * 916) * 1000
dfds['solid'] = dfds['solid'] / (dfd['SA'] * 916) * 1000
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

dfds2 = dfds2.set_index("When")

dfd["label"] = ' '
labels = ["Jan 29", "Feb 05", "Feb 12", "Feb 19", "Feb 26", "Mar 05", "Mar 12", "Mar 19"]
for i in range(0, dfd.shape[0]):
    for item in labels:
        if dfd.When[i] == item:
            dfd.loc[i, 'label'] = dfd.When[i]

dfd = dfd.set_index("When")

fig, (ax1, ax3) = plt.subplots(
    nrows=2, ncols=1, sharex=True, figsize=(15, 12)
)
fig.subplots_adjust(hspace=0)

y1 = dfds2[['Ice frozen', 'Meltwater discharged']]
z = dfd[['$SW_{net}$', '$LW_{net}$', '$Q_S$', '$Q_L$', '$Q_C$']]
y3 = dfd['SA']


y1.plot(kind='bar', stacked=True, edgecolor='black', linewidth=0.5, color=['#D9E9FA', '#0C70DE'], ax=ax1)
ax1.set_ylabel('Thickness Change [$mm$]')

ax1.legend(loc='upper right' , prop={'size': 6})

ax1.grid( axis="y",color="black", alpha=.3, linewidth=.5, which="major")
at = AnchoredText("(a)",
                  prop=dict(size=6), frameon=True,
                  loc='upper left',
                  )
at.patch.set_boxstyle("round,pad=0,rounding_size=0.2")
ax1.add_artist(at)

y3.plot.bar(edgecolor='k',  color=['tab:gray'], linewidth=0.5, ax=ax3)
ax3.set_ylabel('Area [$m^2$]')
ax3.grid(axis="y", color="black", alpha=.3, linewidth=.5, which="major")
at = AnchoredText("(b)",
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

pp.close()
