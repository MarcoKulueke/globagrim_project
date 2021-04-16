import numpy as np
from netCDF4 import Dataset
import os

from .variables import global_const, global_int, global_array


# Later move to variables
def init_output():
    global out, PSG, PHIS, T, U, V

    print("Init: ", os.path.join(os.getcwd(), "output.nc"))
    out = Dataset("output.nc", "w")

    level = out.createDimension("level", global_const.NL)
    time = out.createDimension("time", None)  # unlimited
    lon = out.createDimension("lon", global_const.NJ)
    lat = out.createDimension("lat", global_const.NK)

    # create variables
    longitudes = out.createVariable("lon", np.float, "lon", fill_value=np.nan)
    latitudes = out.createVariable("lat", np.float, "lat", fill_value=np.nan)

    PSG = out.createVariable(
        "PSG", np.float, ("time", "lat", "lon"), fill_value=np.nan
    )  # NetCDF has (level, time, lat, lon) as standard

    PHIS = out.createVariable(
        "PHIS", np.float, ("time", "lat", "lon"), fill_value=np.nan
    )

    T = out.createVariable(
        "T", np.float, ("level", "time", "lat", "lon"), fill_value=np.nan
    )

    U = out.createVariable(
        "U", np.float, ("level", "time", "lat", "lon"), fill_value=np.nan
    )

    V = out.createVariable(
        "V", np.float, ("level", "time", "lat", "lon"), fill_value=np.nan
    )

    #    W = out.createVariable(
    #        "W", np.float, ("level", "time", "lat", "lon"), fill_value=np.nan
    #    )

    #    GP = out.createVariable(
    #        "GP", np.float, ("level", "time", "lat", "lon"), fill_value=np.nan
    #    )

    # assign units
    longitudes.units = "degrees_east"
    latitudes.units = "degrees_north"

    PSG.units = "hPa"
    PHIS.units = "J/kg"
    T.units = "K"
    U.units = "m/s"
    V.units = "m/s"
    #    W.units = "1/s"
    #    GP.units = "J/kg"

    # assign names
    PSG.long_name = "Sea Level Pressure"
    PHIS.long_name = "Sea Level Geopotential"
    T.long_name = "Temperature"
    U.long_name = "Zonal Wind"
    V.long_name = "Meridional Wind"
    #    W.long_name = "Vertical Speed"
    #    GP.long_name = "Geopotential"

    # fill lon/lat
    longitudes[:] = global_array.flam_deg[1 : global_const.NJ + 1]
    latitudes[:] = global_array.phi_deg[1 : global_const.NK + 1]


def fill_output():
    global PSG, U, V

    if global_int.ntout == 0:
        print("Write initial conditions to output.")
    else:
        print("Write to output")

    PSG[global_int.ntout, :, :] = (
        np.swapaxes(
            global_array.psg[1 : global_const.NJ + 1, 1 : global_const.NK + 1], 0, 1
        )
        / 100
    )  # NetCDF has (level, time, lat, lon) as standard

    PHIS[global_int.ntout, :, :] = np.swapaxes(
        global_array.phis[1 : global_const.NJ + 1, 1 : global_const.NK + 1], 0, 1
    )  # NetCDF has (level, time, lat, lon) as standard

    T[:, global_int.ntout, :, :] = (
        np.moveaxis(
            global_array.tw[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
            [0, 1, 2],
            [2, 1, 0],
        )
        + global_const.T0
    )  # NetCDF has (level, time, lat, lon) as standard

    U[:, global_int.ntout, :, :] = np.moveaxis(
        global_array.uw[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
        [0, 1, 2],
        [2, 1, 0],
    )  # NetCDF has (level, time, lat, lon) as standard

    V[:, global_int.ntout, :, :] = np.moveaxis(
        global_array.vw[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
        [0, 1, 2],
        [2, 1, 0],
    )  # NetCDF has (level, time, lat, lon) as standard

    # W[:, global_int.ntout, :, :] = np.moveaxis(
    #     global_array.??[1 : global_const.NJ + 1, 1 : global_const.NK + 1, :],
    #     [0, 1, 2],
    #     [2, 1, 0],
    # )  # NetCDF has (level, time, lat, lon) as standard

    # GP[:, global_int.ntout, :, :] = np.moveaxis(
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
