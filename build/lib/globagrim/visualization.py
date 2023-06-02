import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def plot(var_name, time_step, center, min_max, save):
    from .variables import global_const

    var_set = {"PSG": "jet", "SE": "terrain", "T": "cool", "U": "PiYG", "V": "PiYG"}
    
    if min_max == None:
        levels=10
    else:
        levels=np.linspace(min_max[0], min_max[1], 10)
        
    # open data set
    ds = xr.open_dataset(global_const.output_path)

    # open var
    if len(ds[var_name].dims) == 2:
        plot_var = ds[var_name][:, :]
    elif len(ds[var_name].dims) == 3:
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
        projection=ccrs.Orthographic(
            central_longitude=center[0], central_latitude=center[1]
        ),
    )
    cnplot = plt.contourf(
        lon, ds.lat, plot_var, transform=ccrs.PlateCarree(), levels=levels, cmap=var_set[var_name], extend='both'
    )
    ax.coastlines()
    ax.set_global()

    # Add a color bar
    plt.colorbar(ax=ax, label=ds[var_name].units)

    # Add a title
    plt.title(ds[var_name].long_name + " (output time step " + str(time_step) + ")")

    # ax.quiver(dummy, ds.lat[::4], plot_u, plot_v,
    #       transform=ccrs.PlateCarree(), scale=300, width=0.005)
    # Show
    if save == None:
        plt.show()
    else:
        fig.savefig(save, bbox_inches='tight', dpi=100)
        plt.close()
