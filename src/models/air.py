import pandas as pd
from datetime import datetime
from tqdm import tqdm
import os
import math
import time
from pandas.plotting import register_matplotlib_converters

register_matplotlib_converters()
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import seaborn as sns
from src.data.config import site, dates, option


class Icestupa:
    """Physical Constants"""

    L_s = 2848 * 1000  # J/kg Sublimation
    L_e = 2514 * 1000  # J/kg Evaporation
    L_f = 334 * 1000  # J/kg Fusion
    c_a = 1.01 * 1000  # J/kgC Specific heat air
    c_w = 4.186 * 1000  # J/kgC Specific heat water
    c_i = 2.097 * 1000  # J/kgC Specific heat ice
    rho_w = 1000  # Density of water
    rho_i = 917  # Density of Ice rho_i
    rho_a = 1.29  # kg/m3 air density at mean sea level
    k = 0.4  # Van Karman constant
    k_i = 2.123  # Thermal Conductivity Waite et al. 2006
    bc = 5.670367 * math.pow(10, -8)  # Stefan Boltzman constant
    p0 = 1013  # Standard air pressure hPa

    """Model constants"""
    time_steps = 5*60  # s Model time steps
    dx = 5e-03  # Ice layer thickness
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

    """Surface"""
    ie = 0.95  # Ice Emissivity ie
    a_i = 0.35  # Albedo of Ice a_i
    a_s = 0.85  # Albedo of Fresh Snow a_s
    t_decay = 10  # Albedo decay rate decay_t_d
    z_i = 0.0017  # Ice Momentum and Scalar roughness length
    T_rain = 1  # Temperature condition for liquid precipitation

    """Fountain"""
    dia_f = 0.005  # Fountain aperture diameter
    h_f = 1.35  # Fountain steps h_f
    theta_f = 45  # Fountain aperture diameter
    ftl = 0  # Fountain flight time loss ftl
    T_w = 5  # Fountain Water temperature

    """Miscellaneous"""
    h_aws = 3  # m height of AWS
    latitude = 46.693723
    longitude = 7.297543

    # site = "schwarzsee"
    state = 0
    utc_offset = 1

    def __init__(self, site="schwarzsee"):

        self.site = site

        self.folders = dict(
            input_folder=os.path.join(self.dirname, "data/interim/" + site + "/"),
            output_folder=os.path.join(self.dirname, "data/processed/" + site + "/"),
            sim_folder=os.path.join(
                self.dirname, "data/processed/" + site + "/simulations"
            ),
        )
        if option != "ERA5":
            input_file = self.folders["input_folder"] + "raw_input.csv"
        else:
            input_file = "/home/surya/Programs/PycharmProjects/ERA5/Eispalast_raw_input_ERA5.csv"
            print("Using ERA5")

        self.df = pd.read_csv(input_file, sep=",", header=0, parse_dates=["When"])

        if self.site == "guttannen":
            crit_temp = 0  # Fountain runtime temperature
            self.latitude = 46.649999
            self.longitude = 8.283333
            self.tree_height = 1.93
            self.tree_radius = 4.13 / 2
            self.dia_f = 0.005  # Fountain aperture diameter
            self.h_f = 3.93  # Fountain steps h_f
            self.theta_f = 0
            self.df_cam = pd.read_csv(self.folders["input_folder"] + "cam.csv", sep=",", header=0, parse_dates=["When"])

    def SEA(self, date):

        latitude = self.latitude
        longitude = self.longitude
        utc_offset = self.utc_offset
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

    def projectile_xy(self, v, h=0):
        if h == 0:
            hs = self.h_f
        else:
            hs = h
        g = 9.81
        data_xy = []
        t = 0.0
        theta_f = math.radians(self.theta_f)
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

    def albedo(self, row, s=0, f=0):

        i = row.Index

        """Albedo"""
        # Precipitation
        if (row.Discharge == 0) & (row.Prec > 0):
            if row.T_a < self.T_rain:  # Snow
                s = 0
                f = 0

        if row.Discharge > 0:
            f = 1
            s = 0

        if f == 0:  # last snowed
            self.df.loc[i, "a"] = self.a_i + (self.a_s - self.a_i) * math.exp(
                -s / self.t_decay
            )
            s = s + 1
        else:  # last sprayed

            self.df.loc[i, "a"] = self.a_i

        return s, f

    def spray_radius(self, r_mean=0, dia_f_new=0):

        Area_old = math.pi * math.pow(self.dia_f, 2) / 4
        v_old = self.df["Discharge"].replace(0, np.NaN).mean() / (60 * 1000 * Area_old)

        if r_mean != 0:
            self.r_mean = r_mean
        else:
            if dia_f_new != 0:
                """Keeping Discharge constant"""
                v_new = (
                        math.pi * self.dia_f ** 2 * v_old / (dia_f_new ** 2 * math.pi)
                )
                h_new = h_old - (v_new ** 2 - v_old ** 2) / (2 * 9.81)
                self.r_mean = self.projectile_xy(
                    v=v_new, h=h_new
                )  # todo implement other options
            else:
                self.r_mean = self.projectile_xy(v=v_old)

        return self.r_mean

    def derive_parameters(self):

        missing = ["a", "cld", "SEA", "e_a", "vp_a", "LW_in", "RH"]
        for col in missing:
            if col in list(self.df.columns):
                missing.remove(col)
            else:
                self.df[col] = 0
                print("Missing variables", col)

        """Albedo Decay"""
        self.t_decay = (
                self.t_decay * 24 * 60 * 60 / self.time_steps
        )  # convert to 5 minute time steps
        s = 0
        f = 0

        for row in tqdm(self.df[1:].itertuples(), total=self.df.shape[0]):

            """Solar Elevation Angle"""
            self.df.loc[row.Index, "SEA"] = self.SEA(row.When)

            """ Vapour Pressure"""
            if "vp_a" in missing:
                self.df.loc[row.Index, "vp_a"] = (
                        6.11
                        * math.pow(
                    10,
                    7.5
                    * self.df.loc[row.Index - 1, "T_a"]
                    / (self.df.loc[row.Index - 1, "T_a"] + 237.3),
                )
                        * row.RH
                        / 100
                )

            """LW incoming"""
            if "LW_in" in missing:

                # Cloudiness from diffuse fraction
                if row.SW_direct + row.SW_diffuse > 1:
                    self.df.loc[row.Index, "cld"] = row.SW_diffuse / (
                            row.SW_direct + row.SW_diffuse
                    )
                else:
                    # Night Cloudiness average of last 8 hours
                    if row.Index - 96 > 0:
                        for j in range(row.Index - 96, row.Index):
                            self.df.loc[row.Index, "cld"] += self.df.loc[j, "cld"]
                        self.df.loc[row.Index, "cld"] = (
                                self.df.loc[row.Index, "cld"] / 96
                        )
                    else:
                        for j in range(0, row.Index):
                            self.df.loc[row.Index, "cld"] += self.df.loc[j, "cld"]
                        self.df.loc[row.Index, "cld"] = (
                                self.df.loc[row.Index, "cld"] / row.Index
                        )

                self.df.loc[row.Index, "e_a"] = (
                                                        1.24
                                                        * math.pow(
                                                    abs(self.df.loc[row.Index, "vp_a"] / (row.T_a + 273.15)), 1 / 7
                                                )
                                                ) * (1 + 0.22 * math.pow(self.df.loc[row.Index, "cld"], 2))

                self.df.loc[row.Index, "LW_in"] = (
                        self.df.loc[row.Index, "e_a"]
                        * self.bc
                        * math.pow(row.T_a + 273.15, 4)
                )

            s, f = self.albedo(row, s, f)

        self.df = self.df.round(5)

        if option == "ERA5":
            data_store = pd.HDFStore(self.folders["input_folder"] + "model_input_ERA5.h5")
        else:
            data_store = pd.HDFStore(self.folders["input_folder"] + "model_input.h5")

        data_store["df"] = self.df
        data_store.close()

    def surface_area(self, i):

        if (self.df.solid[i - 1] > 0) & (self.df.loc[i - 1, "r_ice"] > self.r_mean):
            # Ice Radius
            self.df.loc[i, "r_ice"] = self.df.loc[i - 1, "r_ice"]

            # Ice Height
            self.df.loc[i, "h_ice"] = (
                    3
                    * self.df.loc[i - 1, "iceV"]
                    / (math.pi * self.df.loc[i, "r_ice"] ** 2)
            )

            # Height by Radius ratio
            self.df.loc[i, "h_r"] = (
                    self.df.loc[i - 1, "h_ice"] / self.df.loc[i - 1, "r_ice"]
            )

        else:

            # Height to radius ratio
            self.df.loc[i, "h_r"] = self.df.loc[i - 1, "h_r"]

            # Ice Radius
            self.df.loc[i, "r_ice"] = math.pow(
                self.df.loc[i - 1, "iceV"] / math.pi * (3 / self.df.loc[i, "h_r"]),
                1 / 3,
            )

            # Ice Height
            self.df.loc[i, "h_ice"] = self.df.loc[i, "h_r"] * self.df.loc[i, "r_ice"]

        # Area of Conical Ice Surface
        self.df.loc[i, "SA"] = (
                math.pi
                * self.df.loc[i, "r_ice"]
                * math.pow(
            (
                    math.pow(self.df.loc[i, "r_ice"], 2)
                    + math.pow((self.df.loc[i, "h_ice"]), 2)
            ),
            1 / 2,
        )
        )

        self.df.loc[i, "SRf"] = (
                                        0.5
                                        * self.df.loc[i, "h_ice"]
                                        * self.df.loc[i, "r_ice"]
                                        * math.cos(self.df.loc[i, "SEA"])
                                        + math.pi
                                        * math.pow(self.df.loc[i, "r_ice"], 2)
                                        * 0.5
                                        * math.sin(self.df.loc[i, "SEA"])
                                ) / self.df.loc[i, "SA"]

        # # Add Spray radius
        # if self.df.Discharge[i:].sum() != 0:
        #     self.df.loc[i, "SA"] += math.pi * ( self.r_mean **2 - self.df.loc[i, "r_ice"] ** 2)

    def energy_balance(self, row):
        i = row.Index

        self.df.loc[i, "vp_ice"] = 6.112 * np.exp(
            22.46 * (self.df.loc[i - 1, "T_s"]) / ((self.df.loc[i - 1, "T_s"]) + 272.62)
        )

        if self.liquid == 0:
            self.df.loc[i, "Ql"] = (
                    0.623
                    * self.L_s
                    * self.rho_a
                    / self.p0
                    * math.pow(self.k, 2)
                    * self.df.loc[i, "v_a"]
                    * (row.vp_a - self.df.loc[i, "vp_ice"])
                    / ((np.log(self.h_aws / self.z_i)) ** 2)
            )

            if self.df.loc[i, "Ql"] < 0:  # Sublimation
                self.gas -= (
                                    self.df.loc[i, "Ql"] * self.df.loc[i, "SA"] * self.time_steps
                            ) / self.L_s

                # Removing gas quantity generated from previous ice
                self.df.loc[i, "solid"] += (
                                                   self.df.loc[i, "Ql"] * (self.df.loc[i, "SA"]) * self.time_steps
                                           ) / self.L_s

                # Ice Temperature
                self.df.loc[i, "delta_T_s"] += (
                                                       self.df.loc[i, "Ql"] * self.time_steps
                                               ) / (self.rho_i * self.dx * self.c_i)

            else:  # Deposition

                self.df.loc[i, "dpt"] += (
                                                 self.df.loc[i, "Ql"] * self.df.loc[i, "SA"] * self.time_steps
                                         ) / self.L_s

                # Ice Temperature
                self.df.loc[i, "delta_T_s"] += (
                                                       self.df.loc[i, "Ql"] * self.time_steps
                                               ) / (self.rho_i * self.dx * self.c_i)

                self.df.loc[i, "solid"] += self.df.loc[i, "dpt"]

        # Sensible Heat Qs
        self.df.loc[i, "Qs"] = (
                self.c_a
                * self.rho_a
                * row.p_a
                / self.p0
                * math.pow(self.k, 2)
                * self.df.loc[i, "v_a"]
                * (self.df.loc[i, "T_a"] - self.df.loc[i - 1, "T_s"])
                / ((np.log(self.h_aws / self.z_i)) ** 2)
        )

        # Short Wave Radiation SW
        self.df.loc[i, "SW"] = (1 - row.a) * (
                row.SW_direct * self.df.loc[i, "SRf"] + row.SW_diffuse
        )

        # Long Wave Radiation LW
        self.df.loc[i, "LW"] = row.LW_in - self.ie * self.bc * math.pow(
            self.df.loc[i - 1, "T_s"] + 273.15, 4
        )

        if np.isnan(self.df.loc[i, "LW"]):
            print(f"LW {self.df.LW[i]}, LW_in {self.df.LW_in[i]}, T_s {self.df.T_s[i - 1]}")

        if self.liquid > 0:
            self.df.loc[i, "Qf"] = (
                    (
                            self.df.loc[i - 1, "solid"]
                            - self.df.loc[i - 1, "ppt"]
                            - self.df.loc[i - 1, "dpt"]
                    )
                    * self.c_w
                    * self.T_w
                    / (self.time_steps * self.df.loc[i, "SA"])
            )

            self.df.loc[i, "Qf"] += (
                    self.rho_i
                    * self.dx
                    * self.c_i
                    * (self.df.loc[i - 1, "T_s"])
                    / self.time_steps
            )

            self.df.loc[i, "delta_T_s"] = -self.df.loc[i - 1, "T_s"]

        self.df.loc[i, "Qg"] = (
                self.k_i
                * (self.df.loc[i - 1, "T_bulk"] - self.df.loc[i - 1, "T_s"])
                / (
                    self.df.loc[i, "r_ice"]

                )
        )

        # Bulk Temperature
        self.df.loc[i, "T_bulk"] = self.df.loc[i - 1, "T_bulk"] - self.df.loc[
            i, "Qg"
        ] * self.time_steps * self.df.loc[i, "SA"] / (
                                           self.df.loc[i - 1, "ice"] * self.c_i
                                   )

        # Total Energy W/m2
        self.df.loc[i, "TotalE"] = (
                self.df.loc[i, "SW"]
                + self.df.loc[i, "LW"]
                + self.df.loc[i, "Qs"]
                + self.df.loc[i, "Qf"]
                + self.df.loc[i, "Qg"]
        )

        # if np.isnan(self.df.loc[i, "TotalE"]) :
        #     print(f"When {self.df.When[i]}, SW {self.df.SW[i]}, LW {self.df.LW[i]}, Qs {self.df.Qs[i]}, Qf {self.df.Qf[i]}, Qg {self.df.Qg[i]}, SA {self.df.SA[i]}")

        # Total Energy Joules
        self.EJoules = self.df.loc[i, "TotalE"] * self.time_steps * self.df.loc[i, "SA"]

    def summary(self):
        Efficiency = (
                (self.df["meltwater"].iloc[-1] + self.df["ice"].iloc[-1])
                / self.df["input"].iloc[-1]
                * 100
        )

        print("\nIce Volume Max", float(self.df["iceV"].max()))
        print("Fountain efficiency", Efficiency)
        print("Ice Mass Remaining", self.df["ice"].iloc[-1])
        print("Meltwater", self.df["meltwater"].iloc[-1])
        print("Ppt", self.df["ppt"].sum())
        print("Deposition", self.df["dpt"].sum())

        # Full Output
        if option == "ERA5":
            filename4 = self.folders["output_folder"] + "model_results_ERA5.csv"
        else:
            filename4 = self.folders["output_folder"] + "model_results.csv"
        # filename4 = self.folders["output_folder"] + "model_results.csv"
        self.df.to_csv(filename4, sep=",")

        data_store = pd.HDFStore(
            self.folders["output_folder"] + "model_output.h5"
        )
        data_store["df_out"] = self.df
        data_store.close()

    def read_input(self):

        data_store = pd.HDFStore(
            self.folders["input_folder"] + "model_input.h5"
        )
        self.df = data_store["df"]
        data_store.close()

        mask = (self.df["When"] >= dates["start_date"]) & (self.df["When"] <= dates["end_date"])
        self.df = self.df.loc[mask]
        self.df = self.df.reset_index()

    def read_output(self):
        if option == "ERA5":
            self.df = pd.read_csv(
                self.folders["output_folder"] + "model_results_ERA5.csv",
                sep=",",
            )
        else:

            self.df = pd.read_csv(
                self.folders["output_folder"] + "model_results.csv",
                sep=",",
            )
        self.df["When"] = pd.to_datetime(self.df["When"], format="%Y.%m.%d %H:%M:%S")

        print(
            f"Mean of SW {self.df.SW.mean()}, LW {self.df.LW.mean()}, Qs {self.df.Qs.mean()}, Ql {self.df.Ql.mean()}, Qf {self.df.Qf.mean()}, Qg {self.df.Qg.mean()}"
        )
        print(
            f"Range of SW {self.df.SW.min()}-{self.df.SW.max()}, LW {self.df.LW.min()}-{self.df.LW.max()}, Qs {self.df.Qs.min()}-{self.df.Qs.max()}, Ql {self.df.Ql.min()}-{self.df.Ql.max()}, Qf {self.df.Qf.min()}-{self.df.Qf.max()}, Qg {self.df.Qg.min()}-{self.df.Qg.max()}"
        )
        print(
            f"Mean of emissivity {self.df.e_a.mean()}, Range of SRf {self.df.e_a.min()}-{self.df.e_a.max()}"
        )
        print(f"Max SA {self.df.SA.max()}")
        print(
            f"M_input {self.df.input.iloc[-1]}, M_R {self.df.ppt.sum()}, M_D {self.df.dpt.sum()}, M_F {self.df.Discharge.sum() * 5 + self.df.iceV.iloc[0] * self.rho_i}"
        )
        print(
            f"M_U {self.df.unfrozen_water.iloc[-1]}, M_solid {self.df.ice.iloc[-1]}, M_gas {self.df.vapour.iloc[-1]}, M_liquid {self.df.meltwater.iloc[-1]}"
        )
        print(
            f"Max_growth {self.df.solid.max() / 5}, average_discharge {self.df.Discharge.replace(0, np.NaN).mean()}"
        )

        # Output for manim
        filename2 = os.path.join(self.folders["output_folder"], self.site + "_model_gif.csv")
        self.df["h_f"] = self.h_f
        cols = ["When", "h_ice", "h_f", "r_ice", "ice", "T_a", "Discharge"]
        self.df[cols].to_csv(filename2, sep=",")
        print(self.df["When"].iloc[-1])

        # self.corr_plot()

    def melt_freeze(self):

        l = [
            "T_s",  # Surface Temperature
            "T_bulk",  # Bulk Temperature
            "SRf",
            "ice",
            "iceV",
            "solid",
            "vapour",
            "melted",
            "delta_T_s",
            "unfrozen_water",
            "TotalE",
            "SW",
            "LW",
            "Qs",
            "Ql",
            "Qf",
            "Qg",
            "meltwater",
            "SA",
            "h_ice",
            "r_ice",
            "ppt",
            "dpt",
        ]
        for col in l:
            self.df[col] = 0

        self.liquid, self.gas, self.EJoules = [0] * 3

        self.sum_T_s = 0  # weighted_sums
        self.sum_SA = 0  # weighted_sums

        for row in self.df[1:].itertuples():
            i = row.Index

            # Initialize
            if self.df.Discharge[i] > 0 and self.state == 0:
                self.state = 1

                if self.site == 'guttannen':
                    self.spray_radius()
                    self.df.loc[i - 1, "r_ice"] = 9.8655
                    self.df.loc[i - 1, "h_ice"] = self.tree_height

                if self.site == 'schwarzsee':
                    self.df.loc[i - 1, "r_ice"] = self.spray_radius()
                    self.df.loc[i - 1, "h_ice"] = self.dx

                if self.site == 'Leh':
                    self.df.loc[i - 1, "r_ice"] = self.spray_radius()
                    self.df.loc[i - 1, "h_ice"] = self.dx

                self.df.loc[i - 1, "h_r"] = self.df.loc[i - 1, "h_ice"] / self.df.loc[i - 1, "r_ice"]
                self.df.loc[i - 1, "iceV"] = math.pi / 3 * self.df.loc[i - 1, "r_ice"] ** 2 * self.df.loc[
                    i - 1, "h_ice"]
                self.df.loc[i - 1, "ice"] = math.pi / 3 * self.df.loc[i - 1, "r_ice"] ** 2 * self.dx
                self.df.loc[i - 1, "input"] = self.df.loc[i - 1, "ice"]

                self.start = i - 1

            # Ice Melted
            if self.df.loc[i - 1, "ice"] < 0:
                self.df.loc[i - 2, "meltwater"] += self.df.loc[i - 2, "ice"]
                self.df.loc[i - 2, "ice"] = 0
                if self.df.Discharge[i:].sum() == 0:  # If ice melted after fountain run
                    self.df = self.df[self.start:i - 1]
                    self.df = self.df.reset_index(drop=True)
                    break
                else:  # If ice melted in between fountain run
                    self.state = 0

            if self.state == 1:

                self.surface_area(i)

                # Precipitation to ice quantity
                if row.T_a < self.T_rain and row.Prec > 0:
                    self.df.loc[i, "ppt"] = (
                            self.rho_w
                            * row.Prec
                            * math.pi
                            * math.pow(self.df.loc[i, "r_ice"], 2)
                    )

                    self.df.loc[i, "solid"] += self.df.loc[i, "ppt"]

                # Fountain water output
                self.liquid = row.Discharge * (1 - self.ftl) * self.time_steps / 60

                self.energy_balance(row)

                if self.EJoules < 0 and self.liquid > 0:

                    """Freezing water"""

                    self.liquid -= (self.EJoules) / (-self.L_f)

                    if self.liquid < 0:
                        self.liquid += (self.EJoules) / (-self.L_f)
                        self.df.loc[i, "solid"] += self.liquid

                        self.liquid = 0
                    else:
                        self.df.loc[i, "solid"] += (-self.EJoules) / (self.L_f)

                else:
                    # Heating Ice
                    self.df.loc[i, "delta_T_s"] += (
                                                           self.df.loc[i, "TotalE"] * self.time_steps
                                                   ) / (self.rho_i * self.dx * self.c_i)

                    """Hot Ice"""
                    if (self.df.loc[i - 1, "T_s"] + self.df.loc[i, "delta_T_s"]) > 0:
                        # Melting Ice by Temperature
                        self.df.loc[i, "solid"] -= (
                                (self.rho_i * self.dx * self.c_i * self.df.loc[i, "SA"])
                                * (-(self.df.loc[i - 1, "T_s"] + self.df.loc[i, "delta_T_s"]))
                                / (-self.L_f)
                        )

                        self.df.loc[i, "melted"] += (
                                (self.rho_i * self.dx * self.c_i * self.df.loc[i, "SA"])
                                * (-(self.df.loc[i - 1, "T_s"] + self.df.loc[i, "delta_T_s"]))
                                / (-self.L_f)
                        )

                        self.df.loc[i - 1, "T_s"] = 0
                        self.df.loc[i, "delta_T_s"] = 0

                """ Quantities of all phases """
                self.df.loc[i, "T_s"] = (
                        self.df.loc[i - 1, "T_s"] + self.df.loc[i, "delta_T_s"]
                )
                self.df.loc[i, "meltwater"] = (
                        self.df.loc[i - 1, "meltwater"] + self.df.loc[i, "melted"]
                )
                self.df.loc[i, "ice"] = self.df.loc[i - 1, "ice"] + self.df.loc[i, "solid"]
                self.df.loc[i, "vapour"] = self.df.loc[i - 1, "vapour"] + self.gas
                self.df.loc[i, "unfrozen_water"] = (
                        self.df.loc[i - 1, "unfrozen_water"] + self.liquid
                )
                self.df.loc[i, "iceV"] = self.df.loc[i, "ice"] / self.rho_i + self.df.loc[self.start, "iceV"]
                self.df.loc[i, "input"] = (
                        self.df.loc[i - 1, "input"]
                        + self.df.loc[i, "ppt"]
                        + self.df.loc[i, "dpt"]
                        + self.df.loc[i, "Discharge"] * 5
                )
                self.df.loc[i, "thickness"] = self.df.loc[i, "solid"] / (
                        self.df.loc[i, "SA"] * self.rho_i
                )

                self.liquid, self.gas, self.EJoules = [0] * 3

    def corr_plot(self):

        data = self.df

        data = data[data.columns.drop(list(data.filter(regex="Unnamed")))]

        data["$q_{net}$"] = data["TotalE"] + data["Ql"]

        data["$\\Delta M_{input}$"] = data["Discharge"] * 5 + data["dpt"] + data["ppt"]

        data["$SW_{in}$"] = data["SW_direct"] + data["SW_diffuse"]

        # data = data.drop(["When", "input", "ppt", "ice", "T_s", "vapour", "Discharge", "TotalE", "T_a", "SEA", "SW_direct", "a", "cld", "SEA", "e_a", "vp_a", "LW_in", "vp_ice", "SRf", "SW_diffuse", "h_r", "RH", "iceV", "melted", "Qf", "SW", "LW", "Qs", "Ql", "dpt", "p_a", "thickness", "h_ice", "r_ice", "Prec", "v_a", "unfrozen_water", "meltwater"], axis=1)

        data = data.rename(
            {
                "solid": "$\\Delta M_{ice}$",
                "delta_T_s": "$\\Delta T_{ice}$",
                "SA": "A",
                "T_a": "$T_a$",
                "v_a": "$v_a$",
                "p_a": "$p_a$",
            },
            axis=1,
        )

        data = data[["$q_{net}$", "$T_a$", "$v_a$", "$p_a$", "RH", "$SW_{in}$", "$\\Delta M_{ice}$", "A"]]

        print(data.drop("$q_{net}$", axis=1).apply(lambda x: x.corr(data["$q_{net}$"])))

        corr = data.corr()
        ax = sns.heatmap(
            corr,
            vmin=-1,
            vmax=1,
            center=0,
            cmap=sns.diverging_palette(20, 220, n=200),
            square=True,
        )
        ax.set_xticklabels(
            ax.get_xticklabels(), rotation=45, horizontalalignment="right"
        )
        plt.show()


class PDF(Icestupa):
    def print_input(self, filename="derived_parameters.pdf"):
        if filename == "derived_parameters.pdf":
            filename = self.folders["input_folder"]

        """Input Plots"""
        if option == "ERA5":
            filename = self.folders["input_folder"] + "data_ERA5.pdf"
        else:
            filename = self.folders["input_folder"] + "data.pdf"

        pp = PdfPages(filename)

        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(
            nrows=6, ncols=1, sharex="col", sharey="row", figsize=(15, 12)
        )

        x = self.df.When

        y1 = self.df.Discharge
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Fountain Spray [$l\\, min^{-1}$]")
        ax1.grid()

        ax1t = ax1.twinx()
        ax1t.plot(x, self.df.Prec * 1000, "b-", linewidth=0.5)
        ax1t.set_ylabel("Precipitation [$mm$]", color="b")
        for tl in ax1t.get_yticklabels():
            tl.set_color("b")

        y2 = self.df.T_a
        ax2.plot(x, y2, "k-", linewidth=0.5)
        ax2.set_ylabel("Temperature [$\\degree C$]")
        ax2.grid()

        y3 = self.df.SW_direct + self.df.SW_diffuse
        ax3.plot(x, y3, "k-", linewidth=0.5)
        ax3.set_ylabel("Global Rad.[$W\\,m^{-2}$]")
        ax3.grid()

        ax3t = ax3.twinx()
        ax3t.plot(x, self.df.SW_diffuse, "b-", linewidth=0.5)
        ax3t.set_ylim(ax3.get_ylim())
        ax3t.set_ylabel("Diffuse Rad.[$W\\,m^{-2}$]", color="b")
        for tl in ax3t.get_yticklabels():
            tl.set_color("b")

        y4 = self.df.RH
        ax4.plot(x, y4, "k-", linewidth=0.5)
        ax4.set_ylabel("Humidity [$\\%$]")
        ax4.grid()

        y5 = self.df.p_a
        ax5.plot(x, y5, "k-", linewidth=0.5)
        ax5.set_ylabel("Pressure [$hPa$]")
        ax5.grid()

        y6 = self.df.v_a
        ax6.plot(x, y6, "k-", linewidth=0.5)
        ax6.set_ylabel("Wind [$m\\,s^{-1}$]")
        ax6.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        fig.autofmt_xdate()
        plt.savefig(self.folders["output_folder"] + "data.jpg", bbox_inches="tight")
        pp.savefig(bbox_inches="tight")

        plt.clf()

        ax1 = fig.add_subplot(111)
        y1 = self.df.T_a
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Temperature [$\\degree C$]")
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        ax1 = fig.add_subplot(111)

        y2 = self.df.Discharge
        ax1.plot(x, y2, "k-", linewidth=0.5)
        ax1.set_ylabel("Discharge Rate ")
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        ax1 = fig.add_subplot(111)
        y3 = self.df.SW_direct
        ax1.plot(x, y3, "k-", linewidth=0.5)
        ax1.set_ylabel("Direct SWR [$W\\,m^{-2}$]")
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        ax1 = fig.add_subplot(111)
        y31 = self.df.SW_diffuse
        ax1.plot(x, y31, "k-", linewidth=0.5)
        ax1.set_ylabel("Diffuse SWR [$W\\,m^{-2}$]")
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        ax1 = fig.add_subplot(111)
        y4 = self.df.Prec * 1000
        ax1.plot(x, y4, "k-", linewidth=0.5)
        ax1.set_ylabel("Ppt [$mm$]")
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        ax1 = fig.add_subplot(111)
        y5 = self.df.p_a
        ax1.plot(x, y5, "k-", linewidth=0.5)
        ax1.set_ylabel("Pressure [$hPa$]")
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        ax1 = fig.add_subplot(111)
        y6 = self.df.v_a
        ax1.plot(x, y6, "k-", linewidth=0.5)
        ax1.set_ylabel("Wind [$m\\,s^{-1}$]")

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()
        plt.close("all")
        pp.close()

    def print_output(self, filename="model_results.pdf"):

        if filename == "model_results.pdf":
            if option == "ERA5":
                filename = self.folders["output_folder"] + "model_results_ERA5.pdf"
            else:
                filename = self.folders["output_folder"] + "model_results.pdf"

        self.df = self.df.rename(
            {
                "SW": "$q_{SW}$",
                "LW": "$q_{LW}$",
                "Qs": "$q_S$",
                "Ql": "$q_L$",
                "Qf": "$q_{F}$",
                "Qg": "$q_{G}$",
            },
            axis=1,
        )

        # Plots
        pp = PdfPages(filename)

        fig = plt.figure()
        x = self.df.When
        y1 = self.df.iceV

        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Ice Volume [$m^3$]")
        # ax1.set_xlabel("Days")

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
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.SA
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Surface Area [$m^2$]")
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.h_ice
        y2 = self.df.r_ice
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Ice Cone Height [$m$]")

        ax2 = ax1.twinx()
        ax2.plot(x, y2, "b-", linewidth=0.5)
        ax2.set_ylabel("Ice Radius[$m$]", color="b")
        for tl in ax2.get_yticklabels():
            tl.set_color("b")

        # Include Validation line segment 1

        ax1.scatter(datetime(2019, 2, 14, 16), 0.7, color="black", marker="o")
        ax2.scatter(datetime(2019, 2, 14, 16), 1.15, color="blue", marker="o")
        ax1.grid()
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.a
        y2 = self.df.SRf
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Albedo")
        # ax1.set_xlabel("Days")
        ax1t = ax1.twinx()
        ax1t.plot(x, y2, "b-", linewidth=0.5)
        ax1t.set_ylabel("$f_{cone}$", color="b")
        for tl in ax1t.get_yticklabels():
            tl.set_color("b")
        # ax1.grid()

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.T_s
        y2 = self.df.T_bulk
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-", linewidth=0.5, alpha=0.5)
        ax1.set_ylabel("Surface Temperature [$\\degree C$]")
        # ax1.grid()
        ax2 = ax1.twinx()
        ax2.plot(x, y2, "b-", linewidth=0.5)
        ax2.set_ylabel("Bulk Temperature [$\\degree C$]", color="b")
        for tl in ax2.get_yticklabels():
            tl.set_color("b")
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        ax1.set_ylim([-20, 1])
        ax2.set_ylim([-20, 1])
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.close("all")

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(
            nrows=4, ncols=1, sharex="col", sharey="row", figsize=(15, 12)
        )

        x = self.df.When

        y1 = self.df.a
        y2 = self.df.SRf
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Albedo")
        # ax1.set_xlabel("Days")
        ax1t = ax1.twinx()
        ax1t.plot(x, y2, "b-", linewidth=0.5)
        ax1t.set_ylabel("$f_{cone}$", color="b")
        for tl in ax1t.get_yticklabels():
            tl.set_color("b")
        ax1.grid()

        y1 = self.df.e_a
        y2 = self.df.cld
        ax2.plot(x, y1, "k-")
        ax2.set_ylabel("Atmospheric Emissivity")
        # ax1.set_xlabel("Days")
        ax2t = ax2.twinx()
        ax2t.plot(x, y2, "b-", linewidth=0.5)
        ax2t.set_ylabel("Cloudiness", color="b")
        for tl in ax2t.get_yticklabels():
            tl.set_color("b")
        ax2.grid()

        y3 = self.df.vp_a - self.df.vp_ice
        ax3.plot(x, y3, "k-", linewidth=0.5)
        ax3.set_ylabel("Vapour gradient [$hPa$]")
        ax3.grid()

        y4 = self.df.T_a - self.df.T_s
        ax4.plot(x, y4, "k-", linewidth=0.5)
        ax4.set_ylabel("Temperature gradient [$\\degree C$]")
        ax4.grid()

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.close("all")

        dfd = self.df.set_index("When").resample("D").mean().reset_index()
        dfd["Discharge"] = dfd["Discharge"] == 0
        dfd["Discharge"] = dfd["Discharge"].astype(int)
        dfd["Discharge"] = dfd["Discharge"].astype(str)
        dfd["When"] = dfd["When"].dt.strftime("%b %d")

        dfd["label"] = " "
        labels = [
            "Jan 30",
            "Feb 05",
            "Feb 12",
            "Feb 19",
            "Feb 26",
            "Mar 05",
            "Mar 12",
            "Mar 19",
        ]
        for i in range(0, dfd.shape[0]):
            for item in labels:
                if dfd.When[i] == item:
                    dfd.loc[i, "label"] = dfd.When[i]

        dfd = dfd.set_index("label")

        z = dfd[["$q_{SW}$", "$q_{LW}$", "$q_S$", "$q_L$", "$q_{F}$", "$q_{G}$"]]
        ax = z.plot.bar(stacked=True, edgecolor=dfd["Discharge"], linewidth=0.5)
        ax.xaxis.set_label_text("")
        plt.grid(axis="y", color="black", alpha=0.3, linewidth=0.5, which="major")
        # plt.xlabel('Date')
        plt.ylabel("Energy Flux [$W\\,m^{-2}$]")
        plt.legend(loc="lower right")
        # plt.ylim(-150, 150)
        plt.xticks(rotation=45)
        pp.savefig(bbox_inches="tight")
        plt.close("all")

        dfds = self.df[["When", "thickness", "SA"]]

        with pd.option_context("mode.chained_assignment", None):
            dfds["negative"] = dfds.loc[dfds.thickness < 0, "thickness"]
            dfds["positive"] = dfds.loc[dfds.thickness >= 0, "thickness"]

        dfds1 = dfds.set_index("When").resample("D").sum().reset_index()
        dfds2 = dfds.set_index("When").resample("D").mean().reset_index()

        dfds2["When"] = dfds2["When"].dt.strftime("%b %d")
        dfds2["label"] = " "
        labels = [
            "Jan 30",
            "Feb 05",
            "Feb 12",
            "Feb 19",
            "Feb 26",
            "Mar 05",
            "Mar 12",
            "Mar 19",
        ]
        for i in range(0, dfds2.shape[0]):
            for item in labels:
                if dfds2.When[i] == item:
                    dfds2.loc[i, "label"] = dfds2.When[i]

        dfds2 = dfds2.set_index("label")

        dfds1 = dfds1.rename(
            columns={"positive": "Ice thickness", "negative": "Meltwater thickness"}
        )
        y1 = dfds1[["Ice thickness", "Meltwater thickness"]]
        y3 = dfds2["SA"]

        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        y1.plot(
            kind="bar",
            stacked=True,
            edgecolor="black",
            linewidth=0.5,
            color=["#D9E9FA", "#0C70DE"],
            ax=ax1,
        )
        plt.ylabel("Thickness ($m$ w. e.)")
        plt.xticks(rotation=45)
        # plt.legend(loc='upper right')
        ax1.set_ylim(-0.025, 0.025)
        ax1.yaxis.set_minor_locator(AutoMinorLocator())
        ax1.grid(axis="y", color="black", alpha=0.3, linewidth=0.5, which="major")
        x_axis = ax1.axes.get_xaxis()
        x_axis.set_visible(False)

        ax3 = fig.add_subplot(2, 1, 2)
        ax = y3.plot.bar(
            y="SA", linewidth=0.5, edgecolor="black", color="xkcd:grey", ax=ax3
        )
        ax.xaxis.set_label_text("")
        ax3.set_ylabel("Surface Area ($m^2$)")
        ax3.grid(axis="y", color="black", alpha=0.3, linewidth=0.5, which="major")
        plt.xticks(rotation=45)
        # plt.savefig(self.folders["output_folder"] + "thickness.pdf")
        pp.savefig(bbox_inches="tight")
        plt.clf()

        plt.close("all")

        pp.close()

    def print_output_guttannen(self, filename="model_results.pdf"):

        if filename == "model_results.pdf":
            filename = self.folders["output_folder"] + "model_results.pdf"

        self.df = self.df.rename(
            {
                "SW": "$q_{SW}$",
                "LW": "$q_{LW}$",
                "Qs": "$q_S$",
                "Ql": "$q_L$",
                "Qf": "$q_{F}$",
                "Qg": "$q_{G}$",
            },
            axis=1,
        )

        df_temp = pd.read_csv(self.folders["input_folder"] + "lumtemp.csv")

        df_temp["When"] = pd.to_datetime(df_temp["When"])

        # Plots
        pp = PdfPages(filename)

        fig = plt.figure()
        x = self.df.When
        y1 = self.df.iceV

        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Ice Volume [$m^3$]")

        ax1t = ax1.twinx()
        ax1t.plot(self.df_cam.When, self.df_cam.Volume, "o-", color="b", alpha=0.5, linewidth=0.5)
        ax1t.set_ylabel("Cam Volume [$m^3$]", color="b")
        for tl in ax1t.get_yticklabels():
            tl.set_color("b")

        ax1.set_ylim([0, 500])
        ax1t.set_ylim([0, 500])

        ax1.scatter(datetime(2020, 1, 3), (54.15), color="green", marker="o", label="Drone Estimate")
        ax1.scatter(datetime(2020, 1, 24), (120.61), color="green", marker="o")
        ax1.scatter(datetime(2020, 2, 15), (128.32), color="green", marker="o")
        ax1.scatter(datetime(2020, 4, 14, 18), 0, color="green", marker="o")

        ax1.grid()
        ax1.legend()

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.SA
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Surface Area [$m^2$]")
        ax1.grid()

        ax1t = ax1.twinx()
        ax1t.plot(self.df_cam.When, self.df_cam.SA, "o-", color="b", alpha=0.5, linewidth=0.5)
        ax1t.set_ylabel("Cam SA [$m^2$]", color="b")
        for tl in ax1t.get_yticklabels():
            tl.set_color("b")

        ax1.set_ylim([0, 600])
        ax1t.set_ylim([0, 600])

        ax1.scatter(datetime(2020, 1, 3), (334.78), color="green", marker="o", label="Drone Estimate")
        ax1.scatter(datetime(2020, 1, 24), (374.61), color="green", marker="o")
        ax1.scatter(datetime(2020, 2, 15), (564.12), color="green", marker="o")
        ax1.legend()

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.h_ice
        y2 = self.df.r_ice
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Ice Cone Height [$m$]")
        ax1.grid()

        ax2 = ax1.twinx()
        ax2.plot(x, y2, "b-", linewidth=0.5)
        ax2.set_ylabel("Ice Radius[$m$]", color="b")
        for tl in ax2.get_yticklabels():
            tl.set_color("b")

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.h_r
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Ice Cone Slope [$m$]")
        ax1.grid()

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())

        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.T_s
        y2 = self.df.T_bulk
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-", linewidth=0.5, alpha=0.5)
        ax1.set_ylabel("Surface Temperature[$\\degree C$]")
        # ax1.grid()
        ax2 = ax1.twinx()
        ax2.plot(x, y2, "b-", linewidth=0.5)
        ax2.set_ylabel("Bulk Temperature[$\\degree C$]", color="b")
        for tl in ax2.get_yticklabels():
            tl.set_color("b")
        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        ax1.set_ylim([-20, 1])
        ax2.set_ylim([-20, 1])
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.clf()

        y1 = self.df.T_s
        ax1 = fig.add_subplot(111)
        ax1.plot(x, y1, "k-", linewidth=0.5, alpha=0.5)
        ax1.set_ylabel("Surface Temperature[$\\degree C$]")

        ax1t = ax1.twinx()
        ax1t.scatter(df_temp.When, df_temp.Temp, color="b", alpha=0.5, s=1)
        ax1t.set_ylabel("Cam Temp [$\\degree C$]", color="b")
        for tl in ax1t.get_yticklabels():
            tl.set_color("b")

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        ax1.set_ylim([-20, 5])
        ax1t.set_ylim([-20, 5])
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.close("all")

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(
            nrows=4, ncols=1, sharex="col", sharey="row", figsize=(15, 12)
        )

        x = self.df.When

        y1 = self.df.a
        y2 = self.df.SRf
        ax1.plot(x, y1, "k-")
        ax1.set_ylabel("Albedo")
        # ax1.set_xlabel("Days")
        ax1t = ax1.twinx()
        ax1t.plot(x, y2, "b-", linewidth=0.5)
        ax1t.set_ylabel("$f_{cone}$", color="b")
        for tl in ax1t.get_yticklabels():
            tl.set_color("b")
        ax1.grid()

        y1 = self.df.e_a
        y2 = self.df.cld
        ax2.plot(x, y1, "k-")
        ax2.set_ylabel("Atmospheric Emissivity")
        # ax1.set_xlabel("Days")
        ax2t = ax2.twinx()
        ax2t.plot(x, y2, "b-", linewidth=0.5)
        ax2t.set_ylabel("Cloudiness", color="b")
        for tl in ax2t.get_yticklabels():
            tl.set_color("b")
        ax2.grid()

        y3 = self.df.vp_a - self.df.vp_ice
        ax3.plot(x, y3, "k-", linewidth=0.5)
        ax3.set_ylabel("Vapour gradient [$hPa$]")
        ax3.grid()

        y4 = self.df.T_a - self.df.T_s
        ax4.plot(x, y4, "k-", linewidth=0.5)
        ax4.set_ylabel("Temperature gradient [$\\degree C$]")
        ax4.grid()

        ax1.xaxis.set_major_locator(mdates.WeekdayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax1.xaxis.set_minor_locator(mdates.DayLocator())
        fig.autofmt_xdate()
        pp.savefig(bbox_inches="tight")
        plt.close("all")

        dfd = self.df.set_index("When").resample("D").mean().reset_index()
        dfd["Discharge"] = dfd["Discharge"] == 0
        dfd["Discharge"] = dfd["Discharge"].astype(int)
        dfd["Discharge"] = dfd["Discharge"].astype(str)
        dfd["When"] = dfd["When"].dt.strftime("%b %d")

        dfd["label"] = " "

        for i in range(0, dfd.shape[0]):
            if i % 7 == 0:
                dfd.loc[i, "label"] = dfd.When[i]

        dfd = dfd.set_index("label")

        z = dfd[["$q_{SW}$", "$q_{LW}$", "$q_S$", "$q_L$", "$q_{F}$", "$q_{G}$"]]
        ax = z.plot.bar(stacked=True, edgecolor=dfd["Discharge"], linewidth=0.5)
        ax.xaxis.set_label_text("")
        plt.grid(axis="y", color="black", alpha=0.3, linewidth=0.5, which="major")
        # plt.xlabel('Date')
        plt.ylabel("Energy Flux Density [$W\\,m^{-2}$]")
        plt.legend(loc="best")
        plt.xticks(rotation=45)
        pp.savefig(bbox_inches="tight")
        plt.close("all")

        dfds = self.df[["When", "thickness", "SA"]]

        with pd.option_context("mode.chained_assignment", None):
            dfds["negative"] = dfds.loc[dfds.thickness < 0, "thickness"]
            dfds["positive"] = dfds.loc[dfds.thickness >= 0, "thickness"]

        dfds1 = dfds.set_index("When").resample("D").sum().reset_index()
        dfds2 = dfds.set_index("When").resample("D").mean().reset_index()

        dfds2["When"] = dfds2["When"].dt.strftime("%b %d")
        dfds2["label"] = " "

        for i in range(0, dfd.shape[0]):
            if i % 7 == 0:
                dfds2.loc[i, "label"] = dfd.When[i]

        dfds2 = dfds2.set_index("label")

        dfds1 = dfds1.rename(
            columns={"positive": "Ice thickness", "negative": "Meltwater thickness"}
        )
        y1 = dfds1[["Ice thickness", "Meltwater thickness"]]
        y3 = dfds2["SA"]

        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        y1.plot(
            kind="bar",
            stacked=True,
            edgecolor="black",
            linewidth=0.5,
            color=["#D9E9FA", "#0C70DE"],
            ax=ax1,
        )
        plt.ylabel("Thickness ($m$ w. e.)")
        plt.xticks(rotation=45)
        # plt.legend(loc='upper right')
        ax1.set_ylim(-0.055, 0.055)
        ax1.yaxis.set_minor_locator(AutoMinorLocator())
        ax1.grid(axis="y", color="black", alpha=0.3, linewidth=0.5, which="major")
        x_axis = ax1.axes.get_xaxis()
        x_axis.set_visible(False)

        ax3 = fig.add_subplot(2, 1, 2)
        ax = y3.plot.bar(
            y="SA", linewidth=0.5, edgecolor="black", color="xkcd:grey", ax=ax3
        )
        ax.xaxis.set_label_text("")
        ax3.set_ylabel("Surface Area ($m^2$)")
        ax3.grid(axis="y", color="black", alpha=0.3, linewidth=0.5, which="major")
        plt.xticks(rotation=45)
        # plt.savefig(self.folders["output_folder"] + "thickness.pdf")
        pp.savefig(bbox_inches="tight")
        plt.clf()

        plt.close("all")

        pp.close()


if __name__ == "__main__":
    start = time.time()

    schwarzsee = PDF(site=site)

    schwarzsee.derive_parameters()

    schwarzsee.print_input()

    schwarzsee.read_input()

    schwarzsee.melt_freeze()

    # schwarzsee.read_output()

    # schwarzsee.corr_plot()

    schwarzsee.summary()

    schwarzsee.print_output()

    total = time.time() - start

    print("Total time : ", total / 60)
