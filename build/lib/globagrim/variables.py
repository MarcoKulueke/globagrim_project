import numpy as np
from .constants import global_const

#
# string variables
#
class GLOBAL_STR:
    start_time = "18.09.2010 01:55:19"
    # cmon = "dec"  # month
    # cday = "30"  # day
    # chour = "00"  # hour
    # cmin = "00"  # minute
    
global_str = GLOBAL_STR()

#
# integer variables
#
class GLOBAL_INT:
    
    nyear = 2014  # year
    nmonth = 12  # month
    nday = 30  # day
    nhour = 0  # hour
    nmin = 0  # minutes
    #
    # float variables
    #
    ti = 0.0  # model time in seconds
    #
    # scalars
    #
    dy = np.nan
    dsig = np.nan
    ntout = 0  # model output step

global_int = GLOBAL_INT()
#
# arrays
#
class GLOBAL_ARRAY:
    #
    # 1D arrays
    #
    flam_deg = np.full(global_const.NJ + 2, np.nan)  # longitude of grid points (deg)
    phi_deg = np.full(global_const.NK + 2, np.nan)  # latitude of grid points (deg)
    flam = np.full(global_const.NJ + 2, np.nan)  # longitude of grid points (rad)
    phi = np.full(global_const.NK + 2, np.nan)  # latitude of grid points (rad)
    dx = np.full(global_const.NK + 2, np.nan)  # zonal grid point distance (metric)
    cs = np.full(global_const.NK + 2, np.nan)  # cosinus latitude
    sn = np.full(global_const.NK + 2, np.nan)  # sinus latitude
    f = np.full(global_const.NK + 2, np.nan)  # coriolis parameter
    sigma = np.full(global_const.NL + 1, np.nan)  # SIGMA on full levels
    sigmah = np.full(global_const.NL, np.nan)  # SIGMA on half levels
    alp = np.full(global_const.NL, np.nan)  # height dependent coeffizient alpha
    gp0 = np.full(global_const.NL, np.nan)  # reference geopotential
    #
    # 2d arrays
    #
    ps = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2]
    )  # surfance pressure anomaly
    psn = np.zeros([global_const.NJ + 2, global_const.NK + 2])  # ps future
    pst = np.zeros([global_const.NJ + 2, global_const.NK + 2])  # Trend ps
    psg = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2]
    )  # absolute surface pressure
    phis = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2]
    )  # geopotential at surface (orography)
    #
    # 3d arrays with boundary conditions
    #
    u = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # mass-weighted zonal wind
    v = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # mass-weighted meridional wind
    t = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # mass-weighted temperature
    un = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # U future
    vn = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # V future
    tn = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # T future
    psa = np.zeros([global_const.NJ + 2, global_const.NK + 2])  # ps past
    ua = np.zeros([global_const.NJ + 2, global_const.NK + 2, global_const.NL])  # U past
    va = np.zeros([global_const.NJ + 2, global_const.NK + 2, global_const.NL])  # V past
    ta = np.zeros([global_const.NJ + 2, global_const.NK + 2, global_const.NL])  # T past
    ut = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # trend U
    vt = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # trend V
    tt = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # trend T
    uw = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # true zonal wind
    vw = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # true meridional wind
    tw = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # true temperature
    gp = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # geopotential
    ###########################
    # 3d arrays without boundary condiitons (index 0, NK +1 and NJ +1 are left empty)
    ###########################
    dsdt = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # SIGMA vertical velocity
    comp = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # adiabatic compression heat
    apx = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # zonal pressure gradient
    apy = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # meridional pressure gradient
    acx = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # zonal coriolis force
    acy = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # meridional coriolis force
    d = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # divergence of mass-weighted wind
    dm = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # divergence of mass flow
    du = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # zonal divergence of momentum
    dv = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # meridional divergence of momentum
    dvt = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # divergence of temperature flow
    vdivu = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # divergence of vertiacl U-flow
    vdivv = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # divergence of vertiacl V-flow
    vdivt = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # divergence of vertiacl T-flow
    diffu = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # diffusion of zonal momentum
    diffv = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # diffusion of meridional momentum
    difft = np.zeros(
        [global_const.NJ + 2, global_const.NK + 2, global_const.NL]
    )  # diffuson of temperature

global_array = GLOBAL_ARRAY()