import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors
import uncertainpy as un
import statistics as st


input = "/home/surya/Programs/PycharmProjects/air_model/data/processed/schwarzsee/simulations/data/"
figures = "/home/surya/Programs/PycharmProjects/air_model/data/processed/schwarzsee/simulations/figures/"

names = "full"
variance = []
mean = []
evaluations = []

data = un.Data()
filename1 = input + names + ".h5"
data.load(filename1)

eval = data["max_volume"].evaluations
print(f"95 percent confidence interval caused by {names} is {round(2 * st.stdev(eval),2)}")

print(data["max_volume"].mean)

print(data["UQ_Icestupa"])
data = data["UQ_Icestupa"]
fig, ax = plt.subplots()


ax.plot(data.time, data.mean, color = 'black')
ax.set_xlabel("Time [Days]")
ax.set_ylabel("Ice Volume[$m^3$]")


ax.fill_between(data["time"],
                         data.percentile_5,
                         data.percentile_95, color='gray', alpha=0.2)

# ax.set_xlim([min(time), max(time)])
ax.set_ylim([0, 1])
plt.legend(["Mean", "90% prediction interval"], loc="best")

plt.tight_layout()


plt.savefig(figures + "uncertainty.jpg", bbox_inches="tight", dpi=300)