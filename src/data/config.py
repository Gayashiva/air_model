import os
import time
from datetime import datetime

dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# site = input("Input the Field Site Name: ") or "guttannen"

site = "schwarzsee_2020"
option = "schwarzsee"

print("Site is", site)

max = False

surface = dict(
    ie=0.95,  # Ice Emissivity ie
    a_i=0.35,  # Albedo of Ice a_i
    a_s=0.85,  # Albedo of Fresh Snow a_s
    decay_t=10,  # Albedo dry decay rate decay_t_d
    dx=1e-02,   #Ice layer thickness
    z0mi=0.0017,  # Ice Momentum roughness length
    z0hi=0.0017,  # Ice Scalar roughness length
    snow_fall_density= 250, # Snowfall density
    rain_temp=1, # Temperature condition for liquid precipitation
    h_aws = 3,  # m height of AWS
)

# if max :
#     surface = dict(
#         ie=0.9975,  # Ice Emissivity ie
#         a_i=0.36175,  # Albedo of Ice a_i
#         a_s=0.8925,  # Albedo of Fresh Snow a_s
#         decay_t=10.5,  # Albedo dry decay rate decay_t_d
#         dx=0.0015,  # Ice layer thickness
#         z0mi=0.0017,  # Ice Momentum roughness length
#         z0hi=0.0017,  # Ice Scalar roughness length
#         snow_fall_density=250,  # Snowfall density
#         rain_temp=1,  # Temperature condition for liquid precipitation
#         h_aws=3,  # m height of AWS
#
#     )
#
# else :
#     surface = dict(
#         ie=0.9025,  # Ice Emissivity ie
#         a_i=0.3325,  # Albedo of Ice a_i
#         a_s=0.8075,  # Albedo of Fresh Snow a_s
#         decay_t=9.5,  # Albedo dry decay rate decay_t_d
#         dx=0.0005,  # Ice layer thickness
#         z0mi=0.0017,  # Ice Momentum roughness length
#         z0hi=0.0017,  # Ice Scalar roughness length
#         snow_fall_density=250,  # Snowfall density
#         rain_temp=1,  # Temperature condition for liquid precipitation
#         h_aws=3,  # m height of AWS
#     )

if site == "schwarzsee_2020":
    folders = dict(
        dirname=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
        input_folder=os.path.join(dir, "data/interim/schwarzsee_2020/"),
        output_folder=os.path.join(dir, "data/processed/schwarzsee_2020/"),
        sim_folder=os.path.join(dir, "data/processed/schwarzsee_2020/simulations/"),
        data=os.path.join(dir, "data/raw/CardConvert/"),
    )

    dates = dict(
        start_date=datetime(2019, 1, 30, 17),
        end_date=datetime(2019, 3, 20, 18),
        fountain_off_date=datetime(2019, 3, 10, 18),
    )
    fountain = dict(
        aperture_f=0.005,  # Fountain aperture diameter
        h_f=1.35,  # Fountain steps h_f
        discharge=3.58,  # Fountain on discharge
        crit_temp=-5,  # Fountain runtime temperature
        latitude = 46.693723,
        longitude = 7.297543,
        utc_offset = 1,
    )

if site == "schwarzsee":
    folders = dict(
        dirname=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
        input_folder=os.path.join(dir, "data/interim/schwarzsee/"),
        output_folder=os.path.join(dir, "data/processed/schwarzsee/"),
        sim_folder=os.path.join(dir, "data/processed/schwarzsee/simulations/"),
        data_file=os.path.join(dir, "data/raw/" + site + "_aws.txt"),
    )

    dates = dict(
        start_date=datetime(2019, 1, 30, 17),
        end_date=datetime(2019, 3, 20, 18),
        fountain_off_date=datetime(2019, 3, 10, 18),
    )
    fountain = dict(
        aperture_f=0.005,  # Fountain aperture diameter
        h_f=1.35,  # Fountain steps h_f
        discharge=3.58,  # Fountain on discharge
        crit_temp=-5,  # Fountain runtime temperature
        latitude = 46.693723,
        longitude = 7.297543,
        utc_offset = 1,
    )

if site == "plaffeien":
    folders = dict(
        dirname=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
        input_folder=os.path.join(dir, "data/interim/plaffeien/"),
        output_folder=os.path.join(dir, "data/processed/plaffeien/"),
        sim_folder=os.path.join(dir, "data/processed/plaffeien/simulations/"),
        data_file=os.path.join(dir, "data/raw/" + site + "_aws.txt"),
    )

    dates = dict(
        start_date=datetime(2018, 11, 15),
        end_date=datetime(2019, 7, 1),
        fountain_off_date=datetime(2019, 3, 1),
    )
    fountain = dict(
        aperture_f=0.005,  # Fountain hole diameter
        h_f=1,  # Fountain steps h_f
        discharge=10.5,  # Fountain on discharge in LPM
        crit_temp=-1,  # Fountain runtime temperature
    )

if site == "guttannen":
    folders = dict(
        dirname=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
        input_folder=os.path.join(dir, "data/interim/guttannen/"),
        output_folder=os.path.join(dir, "data/processed/guttannen/"),
        data_file=os.path.join(dir, "data/raw/" + site + "_aws.txt"),
    )

    dates = dict(
        start_date=datetime(2018, 1, 10),
        end_date=datetime(2018, 3, 1),
        fountain_off_date=datetime(2018, 2, 1),
    )
    fountain = dict(
        aperture_f=0.005,  # Fountain hole diameter
        h_f=5,  # Fountain steps h_f
        discharge=6,  # Fountain on discharge
        crit_temp=-1,  # Fountain runtime temperature
    )
