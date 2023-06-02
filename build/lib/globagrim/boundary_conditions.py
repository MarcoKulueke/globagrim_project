from .variables import global_const, global_array


def boundary_conditions():
    #
    # poles
    #
    for j in range(0, global_const.NJ + 1):
        joppos = j + int(global_const.NJ / 2)
        if joppos > global_const.NJ:
            joppos = joppos - global_const.NJ
        global_array.ps[j, 0] = global_array.ps[joppos, 1]
        global_array.ps[j, global_const.NK + 1] = global_array.ps[
            joppos, global_const.NK
        ]
        global_array.phis[j, 0] = global_array.phis[joppos, 1]
        global_array.phis[j, global_const.NK + 1] = global_array.phis[
            joppos, global_const.NK
        ]
        global_array.u[j, 0, :] = -global_array.u[joppos, 1, :]
        global_array.u[j, global_const.NK + 1, :] = -global_array.u[
            joppos, global_const.NK, :
        ]
        global_array.v[j, 0, :] = -global_array.v[joppos, 1, :]
        global_array.v[j, global_const.NK + 1, :] = -global_array.v[
            joppos, global_const.NK, :
        ]
        global_array.t[j, 0, :] = global_array.t[joppos, 1, :]
        global_array.t[j, global_const.NK + 1, :] = global_array.t[
            joppos, global_const.NK, :
        ]

    #
    # east/west
    #
    for k in range(0, global_const.NK + 2):
        global_array.ps[0, k] = global_array.ps[global_const.NJ, k]
        global_array.ps[global_const.NJ + 1, k] = global_array.ps[1, k]
        global_array.phis[0, k] = global_array.phis[global_const.NJ, k]
        global_array.phis[global_const.NJ + 1, k] = global_array.phis[1, k]
        global_array.u[0, k, :] = global_array.u[global_const.NJ, k, :]
        global_array.u[global_const.NJ + 1, k, :] = global_array.u[1, k, :]
        global_array.v[0, k, :] = global_array.v[global_const.NJ, k, :]
        global_array.v[global_const.NJ + 1, k, :] = global_array.v[1, k, :]
        global_array.t[0, k, :] = global_array.t[global_const.NJ, k, :]
        global_array.t[global_const.NJ + 1, k, :] = global_array.t[1, k, :]


#########################


if __name__ == "__main__":
    import grid

    grid.grid()

    import init

    init.init_case()

    boundary_conditions()
