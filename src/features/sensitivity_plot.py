import logging
import os
import time
from datetime import datetime
from logging import StreamHandler
from pandas.plotting import register_matplotlib_converters

register_matplotlib_converters()
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from src.data.config import site, option, folders, fountain, surface
from src.models.air_forecast import icestupa

from SALib.sample import saltelli
from SALib.analyze import sobol
import matplotlib.colors


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

problem = {"num_vars": 1, "names": ["dx"], "bounds": [[1e-02, 1e-01]]}
param_values = saltelli.sample(problem, 4, calc_second_order=False)

filename = os.path.join(
        folders['sim_folder'], site + "_simulations_" + str(problem["names"]) + str(param_values.shape[0]) + ".csv"
    )

filename2 = os.path.join(
        folders['sim_folder'], site + "_simulations_" + str(problem["names"]) + "_full" + ".csv"
    )
#  read files
dfo = pd.read_csv(filename, sep=",")
dfo2 = pd.read_csv(filename2, sep=",")

print(dfo)

# x1 = dfo.dx * 1000
y1 = dataframe["temperature"]
print(y1)

fig, ax = plt.subplots()
ax.plot(dataframe)
ax.set_ylabel("Max Ice Volume[$m3$]")
ax.set_xlabel("Ice Layer thickness[mm]")
ax.grid()
plt.show()

# x2 = dfo.dx * 1000
# y2 = dfo.max_melt_thickness * 1000

# x2 = dfo2.dx
# y2 = dfo2.max_melt_thickness * -1

# ax.scatter(x2, y2)
# lims = [
#     np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
#     np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
# ]
#
# # now plot both limits against eachother
# ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
# ax.set_aspect('equal')
# ax.set_xlim(lims)
# ax.set_ylim(lims)

# # Plots
# fig = plt.figure()
# cmap = plt.cm.rainbow
# norm = matplotlib.colors.Normalize(
#     vmin=problem["bounds"][0][0], vmax=problem["bounds"][0][1]
# )
#
# for i in range(dfo.shape[0]):
#     x = dfo.Efficiency[i]
#     y1 = dfo.Max_IceV[i]
#
#     ax1 = fig.add_subplot(111)
#     ax1.scatter(x, y1, color=cmap(norm(dfo.dx[i])))
#     ax1.set_ylabel("Max_IceV[$m^3$]")
#     ax1.set_xlabel("Efficiency")
#
# # format the ticks
# ax1.grid()
#
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = fig.colorbar(sm)
# cbar.set_label("Layer Thickness[$m$]")
# # plt.show()
#
# plt.savefig(
#     os.path.join(
#         folders['sim_folder'], site + "_simulations_" + str(problem["names"]) + str(problem["bounds"])  + ".jpg"
#     ),
#     bbox_inches="tight",
#     dpi=300,
# )
