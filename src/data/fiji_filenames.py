import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
import math
import time
from tqdm import tqdm
import shutil, os
import glob
import fnmatch
from src.data.config import site, dates, folders, fountain, surface
from os import listdir
from os.path import isfile, join
import shutil


dir = "/home/surya/Programs/PycharmProjects/air_model/data/raw/"

oldpath = "/home/surya/Pictures/Guttannen_optical/"
newpath = "/home/surya/Pictures/Guttannen_optical/converted/"
onlyfiles = [f for f in listdir(oldpath) if isfile(join(oldpath, f))]
df_names = pd.DataFrame({"col": onlyfiles})

df_names["Label"] = df_names["col"].str.split("m").str[-1]

df_names["Label"] = (
    "2020-"
    + df_names["Label"].str[2:4]
    + "-"
    + df_names["Label"].str[4:6]
    + " "
    + df_names["Label"].str[6:8]
)

df_names["When"] = pd.to_datetime(df_names["Label"], format="%Y-%m-%d %H")

# df_names["When"] = df_names["When"].dt.strftime('%b %d %H')

print(df_names.head())

for i in range(0, df_names.shape[0]):

    if 6 < df_names.loc[i, "When"].hour < 19 :
        if df_names.loc[i, "When"].hour % 2 == 0:
            shutil.copy(oldpath + df_names.loc[i, "col"], newpath)
            os.rename(newpath + df_names.loc[i, "col"], newpath + str(df_names.loc[i,"When"].strftime('%m-%d %H') + ":00"))
