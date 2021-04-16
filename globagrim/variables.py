import numpy as np

#
# define model constants
#
class global_const:
    pi = np.arctan(1.0) * 4.0
    NJ = 144  # west-east / number of circles of longitude
    NK = 72  # north-south / number of circles of latitude
    NL = 20  # number of SIGMA-surfaces
    DT = 15.0  # time step in SECONDS
    #    TF = 240.0  # total integration time in HOURS
    TF = 0.25  # total integration time in HOURS
    NOUT = 60  # output every NOUT time steps
    PS0 = 100000.0  # reference pressure experiment 1
    RE = 6.371e6  # earth radius
    RHOS = 1.3  # mean air density at surface
    OM = 7.292e-5  # angular velocity of earth
    RD = 287.0  # specific gas constant
    CP = 1005.0  # specific heat capacity
    FKAP = RD / CP  # heat capacity ratio
    T0 = 250.0  # reference temperature
    IEXP = 1  # experiment number
    NTFIL = 8640  # number of filter time steps if lfin=.true.
    FKD = 2.0e5  # diffusion coefficient if LDIFF=.true.
    #    LDIFF = True  # switch for horizontal diffusion
    #    lfil = False  # switch for dynamic filter
    #
    ndmon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


#
# string variables
#
class global_char:
    cmon = "dec"  # month
    cday = "30"  # day
    chour = "00"  # hour
    cmin = "00"  # minute


#
# integer variables
#
class global_int:
    nyear = 2014  # year
    nmonth = 12  # month
    nday = 30  # day
    nhour = 0  # hour
    nmin = 0  # minutes
    #
    # float variables
    #
    ti = 0.0  # model time in HOURS
    #
    # scalars
    #
    dy = np.nan
    dsig = np.nan
    ntout = 0  # model output step


#
# arrays
#
class global_array:
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
    )  # geopotential at surface
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
