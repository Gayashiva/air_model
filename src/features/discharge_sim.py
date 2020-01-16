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

param_values = np.arange(1, 15, 0.2).tolist()


dfx = pd.DataFrame({'Discharge': []})
for i, X in enumerate(param_values):

    #  read files
    filename0 = os.path.join(folders['input_folder'], site + "_" + option + "_input.csv")
    df_in = pd.read_csv(filename0, sep=",")
    df_in["When"] = pd.to_datetime(df_in["When"], format="%Y.%m.%d %H:%M:%S")

    print(X)
    fountain['discharge'] = X
    df = pd.DataFrame({'A': []})
    df = icestupa(df_in, fountain, surface)
    print("Model runtime", df["When"].iloc[-1] - df["When"].iloc[0])

    # dfd = df.set_index("When").resample("D").mean().reset_index()
    # filename1 = os.path.join(
    #     folders['output_folder'], site + "_simulations_" + str(X) + '_' + str(param_values) + ".csv"
    # )
    # dfd.to_csv(filename1, sep =',')

    dfx = dfx.append({'Discharge': X, 'Max Growthrate' : df["solid"].max(), 'Max SA' : df["SA"].max(), 'MaxV': df["iceV"].max(), 'h/r': df["h_r"].iloc[-1], 'r': df["r_ice"].max(), 'Endice': df["iceV"].iloc[-1], 'Runtime': df["When"].iloc[-1] - df["When"].iloc[0]}, ignore_index=True)
    print(dfx)


filename2 = os.path.join(
    folders['output_folder'], site + "_simulations_discharge.csv"
)
print(dfx)
dfx.to_csv(filename2, sep=',')