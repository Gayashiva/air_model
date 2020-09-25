import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from src.data.config import site, option, dates, folders, fountain

file1 = "/home/surya/Programs/PycharmProjects/air_model/data/interim/schwarzsee/raw_input.csv"
file2 = "/home/surya/Programs/PycharmProjects/ERA5/Leh.csv"
# file3 = "/home/surya/Programs/PycharmProjects/ERA5/Others.csv"


df1 = pd.read_csv(file1, sep=",", header=0, parse_dates=["When"])

# Datetime
df = pd.read_csv(file2, sep=",", header=0, parse_dates=["dataDate"])
# df2 = pd.read_csv(file3, sep=",", header=0, parse_dates=["dataDate"])

df = df.drop(["Latitude", "Longitude"], axis=1)
df["time"] = df["validityTime"].replace([0,100,200,300,400, 500, 600,700,800,900], ["0000","0100","0200","0300","0400", "0500", "0600","0700","0800","0900"])
df["time"] = df["time"].astype(str)
df["time"] = df["time"].str[0:2] + ":" + df["time"].str[2:4]
df["dataDate"] = df["dataDate"].astype(str)
df["When"] = df["dataDate"] + " " + df["time"]
df["When"] = pd.to_datetime(df["When"])

start_date=df["When"].iloc[0]
end_date=df["When"].iloc[-1]
days = pd.date_range(start=start_date, end=end_date, freq="1H")
df_out = pd.DataFrame({"When": days})

mask = (df["When"] >= start_date) & (df["When"] <= end_date)
df = df.loc[mask]
df = df.reset_index()

df_out = df_out.set_index("When")
df = df.set_index("When")
df1 = df1.set_index("When")
# df2 = df2.set_index("When")

time_steps = 60*60
df_out["10u"] = df.loc[df.shortName == '10u',"Value"]
df_out["10v"] = df.loc[df.shortName == '10v',"Value"]
df_out["2d"] = df.loc[df.shortName == '2d',"Value"]
df_out["2t"] = df.loc[df.shortName == '2t',"Value"]
df_out["sp"] = df.loc[df.shortName == 'sp',"Value"]
df_out["tcc"] = df.loc[df.shortName == 'tcc',"Value"]
df_out["tp"] = df.loc[df.shortName == 'tp',"Value"]
df_out["ssrd"] = df.loc[df.shortName == 'ssrd',"Value"]/time_steps
df_out["strd"] = df.loc[df.shortName == 'strd',"Value"]/time_steps
df_out["fdir"] = df.loc[df.shortName == 'fdir',"Value"]/time_steps

df_out["v_a"] = np.sqrt(df_out["10u"]**2 + df_out["10v"]**2)
df_out["RH"] = 100*(np.exp((17.625*df_out["2d"])/(243.04+df_out["2d"]))/np.exp((17.625*df_out["2t"])/(243.04+df_out["2t"])))
df_out["sp"] = df_out["sp"]/100
df_out["tp"] = df_out["tp"]/12
df_out["SW_diffuse"] = df_out["ssrd"] - df_out["fdir"]
df_out["2t"] = df_out["2t"]-273.15



# CSV output
df_out.rename(
    columns={
        "2t": "T_a",
        "sp": "p_a",
        "tcc": "cld",
        "tp": "Prec",
        "fdir": 'SW_direct',
        "strd": 'LW_in',

    },
    inplace=True,
)



df = df_out[
            [ "T_a", "RH", "v_a", "SW_direct", "SW_diffuse", "LW_in", "cld", "p_a"]
        ]

df = df.round(5)

upsampled = df.resample("5T")
# interpolated = upsampled.interpolate(method='spline', order=2)
interpolated = upsampled.interpolate(method='linear')
interpolated["Prec"] = df_out["Prec"].resample("5T").bfill()

interpolated = interpolated.reset_index()

interpolated["Discharge"] = 0
mask = (interpolated["T_a"] < fountain["crit_temp"]) & (interpolated["SW_direct"] < 100)
mask_index = interpolated[mask].index
interpolated.loc[mask_index, "Discharge"] = 2*60
mask = interpolated["When"] >= dates["fountain_off_date"]
mask_index = interpolated[mask].index
interpolated.loc[mask_index, "Discharge"] = 0



# For Leh
interpolated["Prec"] = 0
interpolated["RH"] = 20
interpolated["cld"] = 0.1


interpolated.to_csv(folders["input_folder"] + "raw_input.csv")

print(interpolated.tail(10))
# print(df_out.T_a.corr(df1.T_a))
# df1 = df1.reset_index()
# fig, ax = plt.subplots()
# ax.plot(interpolated.When, interpolated.Discharge)
# plt.show()