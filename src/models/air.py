import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from tqdm import tqdm
import os
import math
import time
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib
from matplotlib.offsetbox import AnchoredText
import math
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class Icestupa:


    constants = dict(
        L_s = 2848 * 1000,  # J/kg Sublimation
        L_e = 2514 * 1000, # J/kg Evaporation
        L_f = 334 * 1000,  # J/kg Fusion
        c_w = 4.186 * 1000,  # J/kgC Specific heat water
        c_i = 2.108 * 1000,  # J/kgC Specific heat ice
        rho_w = 1000,  # Density of water
        rho_i = 916,  # Density of Ice rho_i
        rho_a = 1.29,  # kg/m3 air density at mean sea level
        k = 0.4,  # Van Karman constant
        bc = 5.670367 * math.pow(10, -8),  # Stefan Boltzman constant

        time_steps=5 * 60,  # s Model time steps
        vp_w = 6.112, #Saturation Water vapour pressure
        p0 = 1013,  # Standard air pressure hPa
        
    )

    surface = dict(
        ie=0.95,  # Ice Emissivity ie
        a_i=0.35,  # Albedo of Ice a_i
        a_s=0.85,  # Albedo of Fresh Snow a_s
        decay_t=10,  # Albedo decay rate decay_t_d
        dx=1e-02,  # Ice layer thickness
        z0mi=0.0017,  # Ice Momentum roughness length
        z0hi=0.0017,  # Ice Scalar roughness length
        snow_fall_density=250,  # Snowfall density
        rain_temp=1,  # Temperature condition for liquid precipitation
        h_aws=3,  # m height of AWS
    )

    fountain = dict(
        aperture_f=0.005,  # Fountain aperture diameter
        theta_f=45,  # Fountain aperture diameter
        h_f=1.35,  # Fountain steps h_f
        latitude=46.693723,
        longitude=7.297543,
        utc_offset=1,
        ftl=0  # Fountain flight time loss ftl,
    )


    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


    def __init__(self, site ='schwarzsee'):

        self.site = site
        self.folders = dict(
            input_folder=os.path.join(self.dirname, "data/interim/" + site + "/"),
            output_folder=os.path.join(self.dirname, "data/processed/" + site + "/"),
            sim_folder=os.path.join(self.dirname, "data/processed/" + site + "/simulations"),
        )

        input_file = self.folders["input_folder"] + "raw_input.csv"

        self.state = 0
        self.df = pd.read_csv(input_file, sep=",", header=0, parse_dates= ["When"])

        i = 0
        start = 0
        while self.df.loc[i, "Discharge"] == 0:
            start = i
            i += 1

        i = start
        fountain_off = 0
        while self.df.Discharge[i:].any() > 0:
            fountain_off = i
            i += 1

        self.dates = dict(
            start_date=self.df.loc[start, "When"],
            fountain_off_date=self.df.loc[fountain_off, "When"],
        )

        self.df = self.df[start:fountain_off]
        self.df = self.df.reset_index(drop=True)

    def SEA(self, date):

        latitude = self.fountain['latitude']
        longitude = self.fountain['longitude']
        utc_offset = self.fountain['utc_offset']
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

    def projectile_xy(self, v, theta_f, g =9.81 ):
        hs = self.fountain["h_f"]
        data_xy = []
        t = 0.0
        theta_f = math.radians(theta_f)
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

    def albedo(self, i, s = 0, f = 0):

        a_min = self.surface["a_i"]

        """Albedo"""
        # Precipitation
        if (self.df.loc[i, "Discharge"] == 0) & (self.df.loc[i, "Prec"] > 0):
            if self.df.loc[i, "T_a"] < self.surface["rain_temp"]:  # Snow
                s = 0
                f = 0

        if self.df.loc[i, "Discharge"] > 0:
            f = 1
            s = 0

        if f == 0:  # last snowed
            self.df.loc[i, "a"] = self.surface["a_i"] + (self.surface["a_s"] - self.surface["a_i"]) * math.exp(-s / self.surface["decay_t"])
            s = s + 1
        else:  # last sprayed
            self.df.loc[i, "a"] = self.surface["a_i"]

        return s , f

    def derive_parameters_old(self):

        missing = [
            "a",
            "cld",
            "SEA",
            "vp_a",
            "r_f",
            "LW_in",
        ]
        for col in missing:
            if col in list(self.df.columns):
                missing = missing.remove(col)
            else:
                self.df[col] = 0

        """ Fountain Spray radius """
        Area = math.pi * math.pow(self.fountain["aperture_f"], 2) / 4


        """Albedo Decay"""
        self.surface["decay_t"] = (
                self.surface["decay_t"] * 24 * 60 * 60 / self.constants['time_steps']
        )  # convert to 5 minute time steps
        s = 0
        f = 0

        for i in tqdm(range(1, self.df.shape[0])):

            """Solar Elevation Angle"""
            self.df.loc[i, "SEA"] = self.SEA(self.df.loc[i, "When"])

            """ Vapour Pressure"""
            if "vp_a" in missing:
                self.df.loc[i, "vp_a"] = (6.11 * math.pow(10, 7.5 * self.df.loc[i - 1, "T_a"] / (self.df.loc[i - 1, "T_a"] + 237.3)) * self.df.loc[i, "RH"] / 100)

            """LW incoming"""
            if "LW_in" in missing:

                # Cloudiness from diffuse fraction
                if self.df.loc[i, "SW_direct"] + self.df.loc[i, "SW_diffuse"] > 1:
                    self.df.loc[i, "cld"] = self.df.loc[i, "SW_diffuse"] / (
                            self.df.loc[i, "SW_direct"] + self.df.loc[i, "SW_diffuse"]
                    )
                else:
                    # Night Cloudiness average of last 8 hours
                    if i - 96 > 0:
                        for j in range(i - 96, i):
                            self.df.loc[i, "cld"] += self.df.loc[j, "cld"]
                        self.df.loc[i, "cld"] = self.df.loc[i, "cld"] / 96
                    else:
                        for j in range(0, i):
                            self.df.loc[i, "cld"] += self.df.loc[j, "cld"]
                        self.df.loc[i, "cld"] = self.df.loc[i, "cld"] / i

                self.df.loc[i, "e_a"] = ( 1.24 * math.pow(abs(self.df.loc[i, "vp_a"] / (self.df.loc[i, "T_a"] + 273.15)), 1 / 7)
                                   ) * (1 + 0.22 * math.pow(self.df.loc[i, "cld"], 2))

                self.df.loc[i, "LW_in"] = self.df.loc[i, "e_a"] * self.constants['bc'] * math.pow(
                        self.df.loc[i, "T_a"] + 273.15, 4
                    )

            """ Fountain Spray radius """
            v_f = self.df.loc[i, "Discharge"] / (60 * 1000 * Area)
            self.df.loc[i, "r_f"] = self.projectile_xy(v_f, self.fountain['theta_f'])

            s,f = self.albedo(i, s, f)

        self.df = self.df.round(5)

        self.df.to_csv(self.folders["input_folder"] + "model_input.csv")

        self.print_input()

    def print_input(self):

        pp = PdfPages(self.folders["input_folder"] + "derived_parameters.pdf")

        x = self.df.When

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        y1 = self.df.Discharge
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Discharge [$l\, min^{-1}$]")
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

        y1 = self.df.r_f
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Spray Radius [$m$]")
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

        y1 = self.df.e_a
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Atmospheric emissivity")
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
        y1 = self.df.cld
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Cloudiness")
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
        y1 = self.df.vp_a
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Vapour Pressure")
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
        y1 = self.df.a
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Albedo")
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
        y1 = self.df.LW_in
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("LW_in")
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

        """Input Plots"""

        pp = PdfPages(self.folders["input_folder"] +"data" + ".pdf")

        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(
            nrows=6, ncols=1, sharex="col", sharey="row", figsize=(15, 12)
        )

        # fig.suptitle("Field Data", fontsize=14)
        # Remove horizontal space between axes
        # fig.subplots_adjust(hspace=0)

        x = self.df.When

        y1 = self.df.Discharge
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Discharge [$l\, min^{-1}$]")
        ax1.grid()

        ax1t = ax1.twinx()
        ax1t.plot(x, self.df.Prec * 1000, "b-", linewidth=0.5)
        ax1t.set_ylabel("Precipitation [$mm$]", color="b")
        for tl in ax1t.get_yticklabels():
            tl.set_color("b")

        y2 = self.df.T_a
        ax2.plot(x, y2, "k-", linewidth=0.5)
        ax2.set_ylabel("Temperature [$\degree C$]")
        ax2.grid()

        y3 = self.df.SW_direct + self.df.SW_diffuse
        ax3.plot(x, y3, "k-", linewidth=0.5)
        ax3.set_ylabel("Global [$W\,m^{-2}$]")
        ax3.grid()

        ax3t = ax3.twinx()
        ax3t.plot(x, self.df.SW_diffuse, "b-", linewidth=0.5)
        ax3t.set_ylim(ax3.get_ylim())
        ax3t.set_ylabel("Diffuse [$W\,m^{-2}$]", color="b")
        for tl in ax3t.get_yticklabels():
            tl.set_color("b")

        y4 = self.df.RH
        ax4.plot(x, y4, "k-", linewidth=0.5)
        ax4.set_ylabel("Humidity [$\%$]")
        ax4.grid()

        y5 = self.df.p_a
        ax5.plot(x, y5, "k-", linewidth=0.5)
        ax5.set_ylabel("Pressure [$hPa$]")
        ax5.grid()

        y6 = self.df.v_a
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
            os.path.join(self.folders["input_folder"], "data.jpg"),
            bbox_inches="tight",
            dpi=300,
        )

        plt.clf()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        y1 = self.df.T_a
        ax1.plot(x, y1, "k-", linewidth=0.5)
        ax1.set_ylabel("Temperature [$\degree C$]")
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

        y2 = self.df.Discharge
        ax1.plot(x, y2, "k-", linewidth=0.5)
        ax1.set_ylabel("Discharge Rate ")
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

        y3 = self.df.SW_direct
        ax1.plot(x, y3, "k-", linewidth=0.5)
        ax1.set_ylabel("Direct SWR [$W\,m^{-2}$]")
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

        y31 = self.df.SW_diffuse
        ax1.plot(x, y31, "k-", linewidth=0.5)
        ax1.set_ylabel("Diffuse SWR [$W\,m^{-2}$]")
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

        y4 = self.df.Prec * 1000
        ax1.plot(x, y4, "k-", linewidth=0.5)
        ax1.set_ylabel("Ppt [$mm$]")
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

        y5 = self.df.p_a
        ax1.plot(x, y5, "k-", linewidth=0.5)
        ax1.set_ylabel("Pressure [$hPa$]")
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

        y6 = self.df.v_a
        ax1.plot(x, y6, "k-", linewidth=0.5)
        ax1.set_ylabel("Wind [$m\,s^{-1}$]")
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

    def surface_area(self, i):

        if (self.df.Discharge[i] > 0) & (self.df.loc[i - 1, "r_ice"] >= self.r_mean):
            # Ice Radius
            self.df.loc[i, "r_ice"] = self.df.loc[i - 1, "r_ice"]

            # Ice Height
            self.df.loc[i, "h_ice"] = (
                    3 * self.df.loc[i - 1, "iceV"] / (math.pi * self.df.loc[i, "r_ice"] ** 2)
            )

            # Height by Radius ratio
            self.df.loc[i, "h_r"] = self.df.loc[i - 1, "h_ice"] / self.df.loc[i - 1, "r_ice"]

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

        else:

            # Height to radius ratio
            self.df.loc[i, "h_r"] = self.df.loc[i - 1, "h_r"]

            # Ice Radius
            self.df.loc[i, "r_ice"] = math.pow(
                self.df.loc[i - 1, "iceV"] / math.pi * (3 / self.df.loc[i, "h_r"]), 1 / 3
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
                        + math.pow(self.df.loc[i, "r_ice"] * self.df.loc[i, "h_r"], 2)
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
                           ) / (
                                   math.pi
                                   * math.pow(
                               (math.pow(self.df.loc[i, "h_ice"], 2) + math.pow(self.df.loc[i, "r_ice"], 2)),
                               1 / 2,
                           )
                                   * self.df.loc[i, "r_ice"]
                           )

    def energy_balance(self, i):

        self.df.loc[i, "vp_ice"] = 6.112 * np.exp(
            22.46 * (self.df.loc[i - 1, "T_s"]) / ((self.df.loc[i - 1, "T_s"]) + 272.62)
        )

        # Water Boundary
        if self.df.Discharge[i] > 0:
            self.df.loc[i, "vp_s"] = self.constants['vp_w']
            self.L = self.constants['L_e']
            self.c_s = self.constants['c_w']

        else:
            self.df.loc[i, "vp_s"] = self.df.loc[i, "vp_ice"]
            self.L = self.constants['L_s']
            self.c_s = self.constants['c_i']

        self.df.loc[i, "Ql"] = (
                0.623
                * self.L
                * self.constants['rho_a']
                / self.constants['p0']
                * math.pow(self.constants['k'], 2)
                * self.df.loc[i, "v_a"]
                * (self.df.loc[i, "vp_a"] - self.df.loc[i, "vp_s"])
                / (
                        np.log(self.surface["h_aws"] / self.surface["z0mi"])
                        * np.log(self.surface["h_aws"] / self.surface["z0hi"])
                )
        )

        if self.df.loc[i, "Ql"] < 0:
            self.df.loc[i, "gas"] -= (self.df.loc[i, "Ql"] * self.df.loc[i, "SA"] * self.constants['time_steps']) / self.L

            # Removing gas quantity generated from previous ice
            self.df.loc[i, "solid"] += (
                                          self.df.loc[i, "Ql"] * (self.df.loc[i, "SA"]) * self.constants['time_steps']
                                  ) / self.L

            # Ice Temperature
            self.df.loc[i, "delta_T_s"] += (self.df.loc[i, "Ql"] * self.constants['time_steps']) / (
                    self.constants['rho_i'] * self.surface["dx"] * self.c_s
            )

        else:  # Deposition

            self.df.loc[i, "deposition"] += (
                                               self.df.loc[i, "Ql"] * self.df.loc[i, "SA"] * self.constants['time_steps']
                                       ) / self.L

        # Sensible Heat Qs
        self.df.loc[i, "Qs"] = (
                self.c_s
                * self.constants['rho_a']
                * self.df.loc[i, "p_a"]
                / self.constants['p0']
                * math.pow(self.constants['k'], 2)
                * self.df.loc[i, "v_a"]
                * (self.df.loc[i, "T_a"] - self.df.loc[i - 1, "T_s"])
                / (
                        np.log(self.surface["h_aws"] / self.surface["z0mi"])
                        * np.log(self.surface["h_aws"] / self.surface["z0hi"])
                )
        )

        # Short Wave Radiation SW
        self.df.loc[i, "SW"] = (1 - self.df.loc[i, "a"]) * (
                self.df.loc[i, "SW_direct"] * self.df.loc[i, "SRf"] + self.df.loc[i, "SW_diffuse"]
        )

        # Long Wave Radiation LW
        if "LW_in" not in list(self.df.columns):

            self.df.loc[i, "LW"] = self.df.loc[i, "e_a"] * self.constants['bc'] * math.pow(
                self.df.loc[i, "T_a"] + 273.15, 4
            ) - self.surface["ie"] * self.constants['bc'] * math.pow(self.df.loc[i - 1, "T_s"] + 273.15, 4)
        else:
            self.df.loc[i, "LW"] = self.df.loc[i, "LW_in"] - self.surface["ie"] * self.constants['bc'] * math.pow(
                self.df.loc[i - 1, "T_s"] + 273.15, 4
            )

        # Conduction Freezing
        if (self.df.loc[i, "liquid"] > 0) & (self.df.loc[i - 1, "T_s"] < 0):
            self.df.loc[i, "Qc"] = (
                    self.constants['rho_i'] * self.surface["dx"] * self.constants['c_i'] * (-self.df.loc[i - 1, "T_s"]) / self.constants['time_steps']
            )
            self.df.loc[i, "delta_T_s"] = -self.df.loc[i - 1, "T_s"]

        # Total Energy W/m2
        self.df.loc[i, "TotalE"] = (
                self.df.loc[i, "SW"] + self.df.loc[i, "LW"] + self.df.loc[i, "Qs"] + self.df.loc[i, "Qc"]
        )

        # Total Energy Joules
        self.df.loc[i, "EJoules"] = self.df.loc[i, "TotalE"] * self.constants['time_steps'] * self.df.loc[i, "SA"]

    def summary(self,i):

        self.df = self.df[:i]
        Efficiency = float(
            (self.df["meltwater"].tail(1) + self.df["ice"].tail(1))
            / (self.df["sprayed"].tail(1) + self.df["ppt"].sum() + self.df["deposition"].sum())
            * 100
        )

        print("\nIce Volume Max", float(self.df["iceV"].max()))
        print("Fountain efficiency", Efficiency)
        print("Ice Mass Remaining", float(self.df["ice"].tail(1)))
        print("Meltwater", float(self.df["meltwater"].tail(1)))
        print("Ppt", self.df["ppt"].sum())
        print("Model runtime", self.df.loc[i - 1, "When"] - self.df.loc[0, "When"])

        # Full Output
        filename4 = os.path.join(self.folders["output_folder"] + "model_results.csv")
        self.df.to_csv(filename4, sep=",")

        return float(self.df["iceV"].max()), Efficiency

    def melt_freeze(self):

        self.df = pd.read_csv(self.folders["input_folder"] + "model_input.csv", sep=",", header=0, parse_dates= ["When"])

        mask = self.df["When"] <= self.dates['fountain_off_date']
        self.df = self.df.loc[mask]
        self.df = self.df.reset_index(drop = True)

        l = [
            "T_s",  # Surface Temperature
            "delta_T_s",  # Temperature Change
            "ice",
            "iceV",
            "solid",
            "liquid",
            "vapour",
            "melted",
            "gas",
            "water",
            "sprayed",
            "TotalE",
            "SW",
            "LW",
            "Qs",
            "Ql",
            "Qc",
            "meltwater",
            "SA",
            "h_ice",
            "r_ice",
            "SRf",
            "vp_ice",
            "ppt",
            "deposition",
        ]
        for col in l:
            self.df[col] = 0

        """Initialize"""
        self.r_mean = self.df['r_f'].replace(0, np.NaN).mean()
        self.df.loc[0, "r_ice"] = self.r_mean
        self.df.loc[0, "h_ice"] = self.surface["dx"]
        self.df.loc[0, "iceV"] = self.surface["dx"] * math.pi * self.df.loc[0, "r_ice"] ** 2

        for i in tqdm(range(1, self.df.shape[0])):

            # Ice Melted
            if self.df.loc[i - 1, "iceV"] <= 0:
                self.df.loc[i - 1, "solid"] = 0
                self.df.loc[i - 1, "ice"] = 0
                self.df.loc[i - 1, "iceV"] = 0
                if self.df.Discharge[i:].sum() == 0:  # If ice melted after fountain run
                    break
                else:  # If ice melted in between fountain run
                    state = 0

            self.surface_area(i)

            # Precipitation to ice quantity
            if (self.df.loc[i, "T_a"] < self.surface["rain_temp"]) and self.df.loc[i, "Prec"] > 0:

                if self.df.loc[i, 'When'] <= self.dates['fountain_off_date']:
                    self.df.loc[i, "ppt"] = (
                            self.surface["snow_fall_density"]
                            * self.df.loc[i, "Prec"]
                            * math.pi
                            * self.r_mean ** 2)
                else:

                    self.df.loc[i, "ppt"] = (
                            self.surface["snow_fall_density"]
                            * self.df.loc[i, "Prec"]
                            * math.pi
                            * math.pow(self.df.loc[i, "r_ice"], 2)
                    )

            # Fountain water output
            self.df.loc[i, "liquid"] = self.df.loc[i, "Discharge"] * (1 - self.fountain['ftl']) * self.constants['time_steps'] / 60

            self.energy_balance(i)

            if self.df.loc[i, "EJoules"] < 0:

                """ And fountain on """
                if self.df.loc[i - 1, "liquid"] > 0:

                    """Freezing water"""

                    self.df.loc[i, "liquid"] -= (self.df.loc[i, "EJoules"]) / (-self.constants['L_f'])

                    if self.df.loc[i, "liquid"] < 0:
                        self.df.loc[i, "liquid"] += (self.df.loc[i, "EJoules"]) / (-self.constants['L_f'])
                        self.df.loc[i, "solid"] += self.df.loc[i, "liquid"]
                        self.df.loc[i, "liquid"] = 0
                    else:
                        self.df.loc[i, "solid"] += (self.df.loc[i, "EJoules"]) / (-self.constants['L_f'])

                else:
                    """ When fountain off and energy negative """
                    # Cooling Ice
                    self.df.loc[i, "delta_T_s"] += (self.df.loc[i, "TotalE"] * self.constants['time_steps']) / (
                        self.constants['rho_i'] * self.surface["dx"] * self.c_s
                    )

            else:
                # Heating Ice
                self.df.loc[i, "delta_T_s"] += (self.df.loc[i, "TotalE"] * self.constants['time_steps']) / (
                    self.constants['rho_i'] * self.surface["dx"] * self.c_s
                )

                """Hot Ice"""
                if (self.df.loc[i - 1, "T_s"] + self.df.loc[i, "delta_T_s"]) > 0:

                    # Melting Ice by Temperature
                    self.df.loc[i, "solid"] -= (
                        (self.constants['rho_i'] * self.surface["dx"] * self.c_s * self.df.loc[i, "SA"])
                        * (-(self.df.loc[i - 1, "T_s"] + self.df.loc[i, "delta_T_s"]))
                        / (-self.constants['L_f'])
                    )

                    self.df.loc[i, "melted"] += (
                        (self.constants['rho_i'] * self.surface["dx"] * self.c_s * self.df.loc[i, "SA"])
                        * (-(self.df.loc[i - 1, "T_s"] + self.df.loc[i, "delta_T_s"]))
                        / (-self.constants['L_f'])
                    )

                    self.df.loc[i,"thickness"] = self.df.loc[i, 'melted'] / (self.df.loc[i, 'SA'] * self.constants['rho_i'])


                    self.df.loc[i - 1, "T_s"] = 0
                    self.df.loc[i, "delta_T_s"] = 0

            """ Quantities of all phases """
            self.df.loc[i, "T_s"] = self.df.loc[i - 1, "T_s"] + self.df.loc[i, "delta_T_s"]
            self.df.loc[i, "meltwater"] = self.df.loc[i - 1, "meltwater"] + self.df.loc[i, "melted"]
            self.df.loc[i, "ice"] = (
                self.df.loc[i - 1, "ice"]
                + self.df.loc[i, "solid"]
                + self.df.loc[i, "ppt"]
                + self.df.loc[i, "deposition"]
            )
            self.df.loc[i, "vapour"] = self.df.loc[i - 1, "vapour"] + self.df.loc[i, "gas"]
            self.df.loc[i, "sprayed"] = (
                self.df.loc[i - 1, "sprayed"] + self.df.loc[i, "Discharge"] * self.constants['time_steps'] / 60
            )
            self.df.loc[i, "water"] = self.df.loc[i - 1, "water"] + self.df.loc[i, "liquid"]
            self.df.loc[i, "iceV"] = (self.df.loc[i, "ice"] - self.df.loc[i, "ppt"]) / self.constants['rho_i'] + self.df.loc[
                i, "ppt"
            ] / self.surface["snow_fall_density"]


        IceV, Efficiency = self.summary(i)

        self.print_output()

        return IceV, Efficiency

    def print_output(self):

        self.df['melt_thick'] = self.df['melted'] / (self.df['SA'] * 1000)

        self.df = self.df.rename({'SW': '$SW_{net}$', 'LW': '$LW_{net}$', 'Qs': '$Q_S$', 'Ql': '$Q_L$', 'Qc': '$Q_C$'}, axis=1)

        # Plots
        pp = PdfPages(self.folders["output_folder"] + "model_results.pdf")

        x = self.df.When
        y1 = self.df.iceV

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

        y1 = self.df.SA

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

        y1 = self.df.h_ice
        y2 = self.df.r_ice

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

        y1 = self.df.SRf

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

        y1 = self.df.iceV
        y2 = self.df['TotalE'] + self.df['$Q_L$']

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

        y1 = self.df.T_s

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

        y1 = self.df.solid / 5

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

        y1 = self.df.thickness * 1000

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

        y1 = self.df.gas / 5
        y2 = self.df.deposition / 5

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

if __name__ == '__main__':

    start = time.time()

    schwarzsee = Icestupa()

    schwarzsee.derive_parameters()

    # schwarzsee.melt_freeze()

    total = time.time() - start

    print("Total time : ", total / 60)

