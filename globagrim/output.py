import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import os

from .variables import global_const, global_str, global_int, global_array


# Later move to variables
def init_output():
    global out, PSG, SE, T, U, V, time

    print("Init: ", os.path.join(os.getcwd(), global_const.output_path))
    out = Dataset(global_const.output_path, "w")

    level = out.createDimension("level", global_const.NL)
    time = out.createDimension("time", None)  # unlimited
    lon = out.createDimension("lon", global_const.NJ)
    lat = out.createDimension("lat", global_const.NK)

    # create variables
    longitudes = out.createVariable("lon", 'd', "lon", fill_value=np.nan)
    latitudes = out.createVariable("lat", 'd', "lat", fill_value=np.nan)
    time = out.createVariable("time", 'd', "time", fill_value=np.nan)

    SE = out.createVariable("SE", 'f', ("lat", "lon"), fill_value=np.nan)

    PSG = out.createVariable(
        "PSG", 'f', ("time", "lat", "lon"), fill_value=np.nan
    )  # NetCDF has (level, time, lat, lon) as standard

    T = out.createVariable(
        "T", 'f', ("time", "level", "lat", "lon"), fill_value=np.nan
    )

    U = out.createVariable(
        "U", 'f', ("time", "level", "lat", "lon"), fill_value=np.nan
    )

    V = out.createVariable(
        "V", 'f', ("time", "level", "lat", "lon"), fill_value=np.nan
    )

    #    W = out.createVariable(
    #        "W", 'f', ("time", "level", "lat", "lon"), fill_value=np.nan
    #    )

    #    GP = out.createVariable(
    #        "GP", 'f', ("time", "level", "lat", "lon"), fill_value=np.nan
    #    )

    # assign units
    longitudes.units = "degrees_east"
    latitudes.units = "degrees_north"
    time.units = "seconds since 1992-10-8 15:15:42.5"  # Coordinated Universal Time

    SE.units = "m"
    PSG.units = "hPa"
    T.units = "K"
    U.units = "m/s"
    V.units = "m/s"
    #    W.units = "1/s"
    #    GP.units = "J/kg"

    # assign names
    SE.long_name = "Surface Elevation"
    PSG.long_name = "Sea Level Pressure"
    T.long_name = "Temperature"
    U.long_name = "Zonal Wind"
    V.long_name = "Meridional Wind"
    #    W.long_name = "Vertical Speed"
    #    GP.long_name = "Geopotential"

    # fill lon/lat
    longitudes[:] = global_array.flam_deg[1 : global_const.NJ + 1]
    latitudes[:] = global_array.phi_deg[1 : global_const.NK + 1]

    time[0] = (
        datetime.strptime(global_str.start_time, "%d.%m.%Y %H:%M:%S")
        - datetime.strptime("08.10.1992 15:15:42.5", "%d.%m.%Y %H:%M:%S.%f")
    ).total_seconds()


def fill_output():
    #    global PSG, U, V

    if global_int.ntout == 0:
        print("Write initial conditions to output.")
        SE[:, :] = (
            np.swapaxes(
                global_array.phis[1 : global_const.NJ + 1, 1 : global_const.NK + 1],
                0,
                1,
            )
            / global_const.G
        )  # NetCDF has (level, time, lat, lon) as standard
    else:
        print("Write to output")
    
    time[global_int.ntout] = global_int.ti + time[0]

    PSG[global_int.ntout, :, :] = (
        np.swapaxes(
            global_array.psg[1 : global_const.NJ + 1, 1 : global_const.NK + 1], 0, 1
        )
        / 100
    )

    T[global_int.ntout, :, :, :] = (
        np.moveaxis(
            global_array.tw[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
            [0, 1, 2],
            [2, 1, 0],
        )
        + global_const.T0
    )  # NetCDF has (level, time, lat, lon) as standard

    U[global_int.ntout, :, :, :] = np.moveaxis(
        global_array.uw[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
        [0, 1, 2],
        [2, 1, 0],
    )  # NetCDF has (level, time, lat, lon) as standard

    V[global_int.ntout, :, :, :] = np.moveaxis(
        global_array.vw[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
        [0, 1, 2],
        [2, 1, 0],
    )  # NetCDF has (level, time, lat, lon) as standard

    # W[global_int.ntout, :, :, :] = np.moveaxis(
    #     global_array.??[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
    #     [0, 1, 2],
    #     [2, 1, 0],
    # )  # NetCDF has (level, time, lat, lon) as standard

    # GP[global_int.ntout, :, :, :] = np.moveaxis(
    #     global_array.??[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
    #     [0, 1, 2],
    #     [2, 1, 0],
    # )  # NetCDF has (level, time, lat, lon) as standard


def close():
    global out
    # close output and write to file
    out.close()


if __name__ == "__main__":
    import grid

    grid.grid()

    import init

    init.init_case()

    import boundary_conditions

    boundary_conditions.boundary_conditions()

    init_output()
    fill_output()
    close()
