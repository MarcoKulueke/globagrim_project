import numpy as np
from .variables import global_array, global_const

nt = int(global_const.TF * 3600 / global_const.DT + 0.5)  # number of time steps
dtout = global_const.DT * global_const.NOUT / 60  # output interval in min
ntout = nt / global_const.NOUT + 1  # output time steps

#
# define cases (curglobal_const.REntly only one)
#

# Case 1: Initial low pglobal_const.REssuglobal_const.RE system
def low_pressure_system():
    delp = 1000
    phic = global_const.pi / 4
    for k in range(1, global_const.NK +1):
        for j in range(1, global_const.NJ +1):
            global_array.ps[j, k] = -delp * np.exp(
                -(
                    (global_array.flam[j] - global_const.pi) ** 2
                    + (global_array.phi[k] - phic)
                    ** 2  # indexing of flam and phic seems to be wrong
                )
                / 0.05
            )

# Case 2: Stglobal_const.REam over an isolated montain (Williamson test similation)
def montain_flow():
    u0=20.
    phis0=9.81*2000.0
    flamc=3.*global_const.pi/2.
    phic=global_const.pi/6.
    distm=global_const.pi/9.
    for k in range (1,global_const.NK +1):
        for j in range (1,global_const.NJ +1):
            dist=np.sqrt((global_array.flam[j]-flamc)**2+(global_array.phi[k]-phic)**2)
            if(dist < distm):
                global_array.phis[j,k]=phis0*(1.-dist/distm)
            else:
                global_array.phis[j,k]=0.
            global_array.ps[j,k]=-global_const.RHOS*global_array.phis[j,k]-global_const.RHOS*u0*global_const.RE/2.*(2.*global_const.OM+u0/global_const.RE)*global_array.sn[k]**2
            global_array.u[j,k,:]=(global_array.ps[j,k]+global_const.PS0)*u0*global_array.cs[k]
    #
    #     Vorgabe einer Temperaturanomalie
    #
    for k in range(1,global_const.NK +1):
        for j in range(1,global_const.NJ +1):
            dist=np.sqrt((global_array.flam[j]-flamc)**2+(global_array.phi[k]-1.*phic)**2)
            if(dist < distm):
                global_array.t[j,k,:]=(1.-dist/distm)
            else:
                global_array.t[j,k,:]=0.
            global_array.t[j,k,:]=(global_array.ps[j,k]+global_const.PS0)*global_array.t[j,k,:]

#
# initialize case according to selection
#
def init_case():
    if global_const.IEXP == 1:
        low_pressure_system()
    elif global_const.IEXP == 2:
        montain_flow()
    else:
        print("Experiment number not defined.")


###############################################


if __name__ == "__main__":
    import grid

    grid.grid()
    init_case()
