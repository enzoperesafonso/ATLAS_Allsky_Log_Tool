import os
import numpy as np
from astropy.time import Time
import astropy.units as u

_ATLAS_SITE_CODES = ["01k", "02k", "03k", "04k"]
_LOG_ARCHIVE_DIRECTORY_PATH = "archives/"


def seconds_to_days(time: float):
    return time / 86400


def gregorian_to_mjd(time: str):
    t = Time(time, format="iso", scale="ut1")
    t.format = "mjd"
    return t.value


def mjd_to_gregorian(time: float):
    t = Time(time, format="mjd")
    t.format = "iso"
    return t.value


# Check the log files for a specified MJD and Site
def check_log_files(MJD: float, atlas_site_code: str):
    if atlas_site_code not in _ATLAS_SITE_CODES:
        raise ValueError("invalid ATLAS site code")

    root_archive_directory_path = "{}{}".format(
        _LOG_ARCHIVE_DIRECTORY_PATH, atlas_site_code
    )

    mjd_file_names = [str(int(MJD)), str(int(MJD) - 1), str(int(MJD) + 1)]

    for mjd_file_name in mjd_file_names:
        log_file_path = os.path.join(
            root_archive_directory_path,
            mjd_file_name,
            atlas_site_code + mjd_file_name + ".log",
        )

        if os.path.isfile(log_file_path):
            with open(log_file_path, "r") as log_file:
                next(log_file)

                for line in log_file:
                    (
                        observation_number,
                        start_time_mjd,
                        exposure_time,
                    ) = line.strip().split()[:3]

                    if (
                        float(start_time_mjd)
                        < MJD
                        < float(start_time_mjd) + seconds_to_days(float(exposure_time))
                    ):
                        return [observation_number]

    return None

# Get images in a time window specified in UTC
def getLogs(window_start, delta_hours, ATLAS_site_code):  # inefficient but works...
    delta_time = np.linspace(0, delta_hours, delta_hours * 60 * 60) * u.hour
    observing_window = Time(window_start, format="iso", scale="ut1") + delta_time
    observing_window.format = "mjd"

    observing_window_mjd = observing_window.value

    image_filenames = []

    for mjd in observing_window_mjd:
        line = check_log_files(mjd, ATLAS_site_code)

        if line and (line not in image_filenames):
            image_filenames.append(line)

    with open("{}_{}.txt".format(window_start, ATLAS_site_code), "w") as file:
        for item in image_filenames:
            file.write(item[0] + "\n")
