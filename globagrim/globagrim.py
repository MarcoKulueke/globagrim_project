def globagrim():
    from .variables import global_const, global_int, global_array
    from . import grid
    from . import init
    from . import boundary_conditions
    from . import trend
    from . import output

    #
    #     global atmospheric grid point model [GlobAGiM]
    #
    print("Number of longitudes: ", global_const.NJ)
    print("Number of latitudes: ", global_const.NK)
    print("Number of model time steps: ", init.nt)

    #
    #     init model grid
    #
    grid.grid()
    #
    #     init variabales
    #
    init.init_case()
    #
    #     apply boundary conditions
    #
    boundary_conditions.boundary_conditions()
    #
    #     calculate trend
    #
    trend.trend()
    #
    #     init output
    #
    output.init_output()
    #
    #     fill output
    #
    output.fill_output()
    print("---")
    #
    #     first time step with Euler method
    #
    global_array.psn[:, :] = (
        global_array.ps[:, :] + global_const.DT * global_array.pst[:, :]
    )
    global_array.un[:, :, :] = (
        global_array.u[:, :, :] + global_const.DT * global_array.ut[:, :, :]
    )
    global_array.vn[:, :, :] = (
        global_array.v[:, :, :] + global_const.DT * global_array.vt[:, :, :]
    )
    global_array.tn[:, :, :] = (
        global_array.t[:, :, :] + global_const.DT * global_array.tt[:, :, :]
    )
    global_array.ti = global_int.ti + global_const.DT / 3600.0
    n = 0
    global_int.ntout += 1
    #
    #     rewrite results
    #
    global_array.psa[:, :] = global_array.ps[:, :]
    global_array.ua[:, :, :] = global_array.u[:, :, :]
    global_array.va[:, :, :] = global_array.v[:, :, :]
    global_array.ta[:, :, :] = global_array.t[:, :, :]
    global_array.ps[:, :] = global_array.psn[:, :]
    global_array.u[:, :, :] = global_array.un[:, :, :]
    global_array.v[:, :, :] = global_array.vn[:, :, :]
    global_array.t[:, :, :] = global_array.tn[:, :, :]
    #
    print(
        "Model time step: ", n, ", Elapsed model time: ", (n + 1) * global_const.DT / 60, " minutes"
    )
    output.fill_output()
    print("---")
    #
    #     time loop
    #
    for n in range(1, init.nt):
        #
        #       calculate trend
        #
        trend.trend()
        #
        #       time step with Leap-Frog
        #
        global_array.psn[:, :] = (
            global_array.psa[:, :] + 2.0 * global_const.DT * global_array.pst[:, :]
        )
        global_array.un[:, :, :] = (
            global_array.ua[:, :, :] + 2.0 * global_const.DT * global_array.ut[:, :, :]
        )
        global_array.vn[:, :, :] = (
            global_array.va[:, :, :] + 2.0 * global_const.DT * global_array.vt[:, :, :]
        )
        global_array.tn[:, :, :] = (
            global_array.ta[:, :, :] + 2.0 * global_const.DT * global_array.tt[:, :, :]
        )
        global_array.ti = global_int.ti + global_const.DT / 3600.0
        #
        #       rewrite Results
        #
        global_array.psa[:, :] = global_array.ps[:, :]
        global_array.ua[:, :, :] = global_array.u[:, :, :]
        global_array.va[:, :, :] = global_array.v[:, :, :]
        global_array.ta[:, :, :] = global_array.t[:, :, :]
        global_array.ps[:, :] = global_array.psn[:, :]
        global_array.u[:, :, :] = global_array.un[:, :, :]
        global_array.v[:, :, :] = global_array.vn[:, :, :]
        global_array.t[:, :, :] = global_array.tn[:, :, :]
        #
        #       apply boundary conditions
        #
        boundary_conditions.boundary_conditions()
        #
        print(
            "Model time step: ",
            n,
            ", Time: ",
            (n + 1) * global_const.DT / 60,
            " minutes",
        )

        #        if n % global_const.NOUT == 0:
        if 0 == 0:  # output every model time step
            global_int.nmin = global_int.nmin + int(init.dtout + 0.5)
            global_int.ntout += 1
            output.fill_output()

        print("---")

    output.close()


if __name__ == "__main__":
    globagrim()
