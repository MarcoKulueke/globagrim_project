import numpy as np
import os


class global_const:
    pi = np.arctan(1.0) * 4.0
    NJ = 144  # west-east / number of circles of longitude
    NK = 72  # north-south / number of circles of latitude
    NL = 20  # number of SIGMA-surfaces
    DT = 15.0  # time step in SECONDS
    #    TF = 240.0  # total integration time in HOURS
    TF = 0.2  # total integration time in HOURS
    NOUT = 60  # output every NOUT time steps
    PS0 = 100000.0  # reference pressure experiment 1
    RE = 6.371e6  # earth radius
    RHOS = 1.3  # mean air density at surface
    OM = 7.292e-5  # angular velocity of earth
    G = 9.81  # m/s2
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
    output_path = os.path.join("output.nc")
    ndmon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
