import numpy as np
from variables import global_array, global_const

nt = int(global_const.TF * 3600 / global_const.DT + 0.5)  # number of time steps
print("Number of model time steps: ", nt)
dtout = global_const.DT * global_const.NOUT / 60  # output interval in min
ntout = nt / global_const.NOUT + 1  # output time steps

#
# define cases (currently only one)
#

# Case 1
def low_pressure_system():
    delp = 1000
    phic = global_const.pi / 4
    for k in range(0, global_const.NK):
        for j in range(0, global_const.NJ):
            global_array.ps[j, k] = -delp * np.exp(
                -(
                    (global_array.flam[j] - global_const.pi) ** 2
                    + (global_array.phi[k] - phic)
                    ** 2  # indexing of flam and phic seems to be wrong
                )
                / 0.05
            )
    global_array.phis = -global_array.ps / global_const.RHOS


#
# initialize case study according to selection
#
def init_case():
    if global_const.IEXP == 1:
        low_pressure_system()
    else:
        print("Experiment number not defined.")


###############################################


if __name__ == "__main__":
    import grid

    grid.grid()
    init_case()
