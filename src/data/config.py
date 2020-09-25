import os
from datetime import datetime

dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

site = "Leh"
option = ''

# site = input("Input the Field Site Name: ") or site

print("Site is", site)

if site == "schwarzsee":
    folders = dict(
        dirname=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
        input_folder=os.path.join(dir, "data/interim/schwarzsee/"),
        output_folder=os.path.join(dir, "data/processed/schwarzsee/"),
        sim_folder=os.path.join(dir, "data/processed/schwarzsee/simulations/"),
        data_file=os.path.join(dir, "data/raw/schwarzsee/" + site + "_aws.txt"),
    )

    dates = dict(
        start_date=datetime(2019, 1, 30, 17),
        end_date=datetime(2019, 3, 17),
        fountain_off_date=datetime(2019, 3, 10, 18),
    )
    fountain = dict(
        aperture_f=0.005,  # Fountain aperture diameter
        theta_f=45,  # Fountain aperture diameter
        h_f=1.35,  # Fountain steps h_f
        discharge=3.58,  # Fountain on discharge
        crit_temp=0,  # Fountain runtime temperature
        latitude=46.693723,
        longitude=7.297543,
        utc_offset=1,
    )

if site == "Leh":
    folders = dict(
        dirname=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
        input_folder=os.path.join(dir, "data/interim/Leh/"),
        output_folder=os.path.join(dir, "data/processed/Leh/"),
        sim_folder=os.path.join(dir, "data/processed/Leh/simulations/"),
    )

    dates = dict(
        start_date=datetime(2019, 1, 1),
        end_date=datetime(2019, 6, 30),
        fountain_off_date=datetime(2019, 5, 1),
    )
    fountain = dict(
        aperture_f=0.005,  # Fountain aperture diameter
        theta_f=45,  # Fountain aperture diameter
        h_f=5,  # Fountain steps h_f
        discharge=3.58,  # Fountain on discharge
        crit_temp=0,  # Fountain runtime temperature
        latitude=46.693723,
        longitude=7.297543,
        utc_offset=1,
    )

if site == "schwarzsee_2020":
    folders = dict(
        dirname=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
        input_folder=os.path.join(dir, "data/interim/schwarzsee_2020/"),
        output_folder=os.path.join(dir, "data/processed/schwarzsee_2020/"),
        sim_folder=os.path.join(dir, "data/processed/schwarzsee_2020/simulations/"),
        data=os.path.join(dir, "data/raw/schwarzsee/CR6_DATA/CardConvert/"),
    )

    dates = dict(
        start_date=datetime(2020, 2, 15),
        end_date=datetime(2020, 2, 18),
        fountain_off_date=datetime(2020, 2, 10),
    )
    fountain = dict(
        aperture_f=0.005,  # Fountain aperture diameter
        h_f=4,  # Fountain steps h_f
        discharge=3.58,  # Fountain on discharge
        crit_temp=-5,  # Fountain runtime temperature
        latitude=46.693723,
        longitude=7.297543,
        utc_offset=1,
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
        data_file=os.path.join(dir, "data/raw/guttannen/" + site + "_2020.txt"),
    )

    dates = dict(
        start_date=datetime(2020, 1, 1, 18),
        end_date=datetime(2020, 5, 1),
        error_date=datetime(2020, 1, 19),
        fountain_off_date=datetime(2020, 3, 1),
    )
    fountain = dict(
        aperture_f=0.005,  # Fountain aperture diameter
        theta_f=0,
        h_f=3.93,  # Fountain steps h_f
        discharge=3.58,  # Fountain on discharge
        crit_temp=1,  # Fountain runtime temperature
        latitude=46.649999,
        longitude=8.283333,
        utc_offset=1,
        tree_height=1.93,
        tree_radius=4.13 / 2,
    )
