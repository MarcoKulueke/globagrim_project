import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
import xarray as xr


def plot(var_name, time_step):

    colormap = {"PSG": "jet", "PHIS": "terrain", "T": "cool", "U": "PiYG", "V": "PiYG"}
    # open data set
    ds = xr.open_dataset("output.nc")

    # open var
    if len(ds[var_name].dims) == 3:
        plot_var = ds[var_name][time_step, :, :]
    elif len(ds[var_name].dims) == 4:
        plot_var = ds[var_name][time_step, 0, :, :]

    # add cyclic point
    plot_var, lon = add_cyclic_point(plot_var, ds.lon)
    # plot_u, dummy = add_cyclic_point(ds.U[time_step,0,::4,::4], ds.lon[::4])
    # plot_v, dummy = add_cyclic_point(ds.V[time_step,0,::4,::4], ds.lon[::4])

    # plot var
    fig = plt.figure()
    ax = fig.add_subplot(
        1,
        1,
        1,
        projection=ccrs.Orthographic(central_longitude=-90.0, central_latitude=0.0),
    )

    plt.contourf(
        lon, ds.lat, plot_var, transform=ccrs.PlateCarree(), cmap=colormap[var_name]
    )
    ax.coastlines()
    ax.set_global()

    # Add a color bar
    plt.colorbar(ax=ax, label=ds[var_name].units)

    # Add a title
    plt.title(ds[var_name].long_name + " (model time step: " + str(time_step) + ")")

    # ax.quiver(dummy, ds.lat[::4], plot_u, plot_v,
    #       transform=ccrs.PlateCarree(), scale=300, width=0.005)
    # Show
    plt.show()
