from variables import global_array, global_const, global_int
import numpy as np


def grid():
    #
    # determine longitude and latitude and add boundary conditions
    #
    global_array.dlam = 360 / global_const.NJ
    global_array.flam_deg[1 : global_const.NJ + 1] = np.arange(
        0, 360, global_array.dlam
    )
    global_array.flam_deg[0] = global_array.flam_deg[global_const.NJ]
    global_array.flam_deg[global_const.NJ + 1] = global_array.flam_deg[1]
    #
    dphi = 180 / global_const.NK
    global_array.phi_deg[1 : global_const.NK + 1] = np.arange(
        -90 + 0.5 * dphi, 90 + 0.5 * dphi, dphi
    )
    global_array.phi_deg[0] = global_array.phi_deg[1]
    global_array.phi_deg[global_const.NK + 1] = global_array.phi_deg[global_const.NK]
    #
    # convert to radian measuglobal_const.RE
    #
    global_array.dlam = global_const.pi / 180.0 * global_array.dlam
    global_array.flam[:] = global_const.pi / 180.0 * global_array.flam_deg
    dphi = global_const.pi / 180.0 * dphi
    global_array.phi[:] = global_const.pi / 180.0 * global_array.phi_deg
    #
    # calculation of metric grid width
    #
    global_array.dx[:] = (
        global_const.RE * np.cos(global_array.phi) * global_array.dlam
    )  # dependent on latitude
    global_array.dy = global_const.RE * dphi
    #
    # calculation latitide cosine and sine
    #
    global_array.cs[:] = np.cos(global_array.phi)
    global_array.sn[:] = np.sin(global_array.phi)
    #
    # calculation coriolis parameter
    #
    global_array.f[:] = 2 * global_const.OM * global_array.sn[0 : global_const.NK + 2]
    #
    # calculation vertical cooglobal_const.RDinate SIGMA
    #
    #
    # increment
    #
    global_int.dsig = 1.0 / global_const.NL
    #
    # full SIGMA surfaces
    #
    global_array.sigma[:] = (np.arange(global_const.NL + 1) + 1) * global_int.dsig
    #
    # half SIGMA surfaces
    #
    global_array.sigmah[:] = (np.arange(global_const.NL) + 0.5) * global_int.dsig
    #
    # hight dependent coefficient alpha
    #
    global_array.alp[:] = np.full(global_const.NL, np.nan)
    for l in range(0, global_const.NL - 1):
        global_array.alp[l] = 0.5 * np.log(
            global_array.sigmah[l + 1] / global_array.sigmah[l]
        )
    global_array.alp[global_const.NL - 1] = np.log(
        1.0 / global_array.sigmah[global_const.NL - 1]
    )

    #
    # calculation of geopotential
    #
    global_array.gp0[global_const.NL - 1] = (
        global_const.RD * global_const.T0 * global_array.alp[global_const.NL - 1]
    )
    for l in range(global_const.NL - 2, -1, -1):
        global_array.gp0[l] = (
            global_array.gp0[l + 1]
            + 2.0 * global_const.RD * global_const.T0 * global_array.alp[l]
        )


####################################################


if __name__ == "__main__":
    grid()
