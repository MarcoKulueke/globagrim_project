from .variables import global_const, global_int, global_array

#                             #
# Berechnen der Zeittendenzen #
#                             #


def true_wind_and_abs_ps():
    #
    # calculation of true wind component and absolute surface pressure
    #
    for k in range(0, global_const.NK + 2):
        for j in range(0, global_const.NJ + 2):
            global_array.psg[j, k] = global_array.ps[j, k] + global_const.PS0
            global_array.uw[j, k, :] = global_array.u[j, k, :] / global_array.psg[j, k]
            global_array.vw[j, k, :] = global_array.v[j, k, :] / global_array.psg[j, k]
            global_array.tw[j, k, :] = global_array.t[j, k, :] / global_array.psg[j, k]


def geopential():
    #
    # calculation of geopotential with hydrostatic equation
    #
    global_array.gp[:, :, global_const.NL - 1] = (
        global_array.phis[:, :]
        + global_const.RD
        * global_array.tw[:, :, global_const.NL - 1]
        * global_array.alp[global_const.NL - 1]
    )
    for l in range(global_const.NL - 2, -1, -1):
        global_array.gp[:, :, l] = (
            global_array.gp[:, :, l + 1]
            + global_const.RD
            * (global_array.tw[:, :, l] + global_array.tw[:, :, l + 1])
            * global_array.alp[l]
        )


def zonal_pressure_gradient_force():
    #
    # zonal pressure gradient force
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.apx[j, k, :] = (
                -global_const.RD
                * (global_array.tw[j, k, :] + global_const.T0)
                * (global_array.ps[j + 1, k] - global_array.ps[j - 1, k])
                / global_array.dx[k]
                / 2.0
                - global_array.psg[j, k]
                * (global_array.gp[j + 1, k, :] - global_array.gp[j - 1, k, :])
                / global_array.dx[k]
                / 2.0
            )


def meridional_pressure_gradient_force():
    #
    # meridional pressure gradient force
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.apy[j, k, :] = (
                -global_const.RD
                * (global_array.tw[j, k, :] + global_const.T0)
                * (global_array.ps[j, k + 1] - global_array.ps[j, k - 1])
                / global_array.dy
                / 2.0
                - global_array.psg[j, k]
                * (global_array.gp[j, k + 1, :] - global_array.gp[j, k - 1, :])
                / global_array.dy
                / 2.0
            )


def corioles_and_centrifugal_force():
    #
    # coriolis and centrifugal force
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.acx[j, k, :] = (
                global_array.f[k]
                + global_array.uw[j, k, :]
                * global_array.sn[k]
                / global_array.cs[k]
                / global_const.RE
            ) * global_array.v[j, k, :]
            global_array.acy[j, k, :] = (
                -(
                    global_array.f[k]
                    + global_array.uw[j, k, :]
                    * global_array.sn[k]
                    / global_array.cs[k]
                    / global_const.RE
                )
                * global_array.u[j, k, :]
            )


def div_zonal_impulse():
    #
    # zonal divergence of momentum
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.du[j, k, :] = (
                (
                    (global_array.u[j + 1, k, :] + global_array.u[j, k, :])
                    * (global_array.uw[j + 1, k, :] + global_array.uw[j, k, :])
                    - (global_array.u[j, k, :] + global_array.u[j - 1, k, :])
                    * (global_array.uw[j, k, :] + global_array.uw[j - 1, k, :])
                )
                / 4.0
                / global_array.dx[k]
            )
            +(
                (
                    global_array.v[j, k + 1, :] * global_array.cs[k + 1]
                    + global_array.v[j, k, :] * global_array.cs[k]
                )
                * (global_array.uw[j, k + 1, :] + global_array.uw[j, k, :])
                - (
                    global_array.v[j, k, :] * global_array.cs[k]
                    + global_array.v[j, k - 1, :] * global_array.cs[k - 1]
                )
                * (global_array.uw[j, k, :] + global_array.uw[j, k - 1, :])
            ) / 4.0 / global_array.dy / global_array.cs[k]


def div_meridional_impulse():
    #
    # meridional divergence of momentum
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.dv[j, k, :] = (
                (
                    (global_array.u[j + 1, k, :] + global_array.u[j, k, :])
                    * (global_array.vw[j + 1, k, :] + global_array.vw[j, k, :])
                    - (global_array.u[j, k, :] + global_array.u[j - 1, k, :])
                    * (global_array.vw[j, k, :] + global_array.vw[j - 1, k, :])
                )
                / 4.0
                / global_array.dx[k]
            )
            +(
                (
                    global_array.v[j, k + 1, :] * global_array.cs[k + 1]
                    + global_array.v[j, k, :] * global_array.cs[k]
                )
                * (global_array.vw[j, k + 1, :] + global_array.vw[j, k, :])
                - (
                    global_array.v[j, k, :] * global_array.cs[k]
                    + global_array.v[j, k - 1, :] * global_array.cs[k - 1]
                )
                * (global_array.vw[j, k, :] + global_array.vw[j, k - 1, :])
            ) / 4.0 / global_array.dy / global_array.cs[k]


def div_temperature_flow():
    #
    # divergence of temperature flow
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.dvt[j, k, :] = (
                (
                    (global_array.u[j + 1, k, :] + global_array.u[j, k, :])
                    * (global_array.tw[j + 1, k, :] + global_array.tw[j, k, :])
                    - (global_array.u[j, k, :] + global_array.u[j - 1, k, :])
                    * (global_array.tw[j, k, :] + global_array.tw[j - 1, k, :])
                )
                / 4.0
                / global_array.dx[k]
            )
            +(
                (
                    global_array.v[j, k + 1, :] * global_array.cs[k + 1]
                    + global_array.v[j, k, :] * global_array.cs[k]
                )
                * (global_array.tw[j, k + 1, :] + global_array.tw[j, k, :])
                - (
                    global_array.v[j, k, :] * global_array.cs[k]
                    + global_array.v[j, k - 1, :] * global_array.cs[k - 1]
                )
                * (global_array.tw[j, k, :] + global_array.tw[j, k - 1, :])
            ) / 4.0 / global_array.dy / global_array.cs[k]


def div_weighted_wind():
    #
    # divergence of mass-wighted wind
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.d[j, k, :] = (
                (global_array.u[j + 1, k, :] - global_array.u[j - 1, k, :])
                / 2.0
                / global_array.dx[k]
                + (
                    global_array.v[j, k + 1, :] * global_array.cs[k + 1]
                    - global_array.v[j, k - 1, :] * global_array.cs[k - 1]
                )
                / global_array.cs[k]
                / global_array.dy
                / 2.0
            )


def sigma_flow():
    #
    # divergence of mass flow between SIGMA=0 and SIGMA=SIGMA[l]
    #
    global_array.dm[:, :, 1] = global_array.d[:, :, 1] * global_int.dsig
    for l in range(2, global_const.NL):
        global_array.dm[:, :, l] = (
            global_array.dm[:, :, l - 1] + global_array.d[:, :, l] * global_int.dsig
        )


# maybe index error
def vert_speed_sigma():
    #
    # calculation of vertial velocity in SIGMA system
    #
    for l in range(0, global_const.NL - 2):
        global_array.dsdt[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l] = (
            global_array.sigma[l]
            / global_array.psg[1 : global_const.NJ + 1, 1 : global_const.NK + 1]
            * global_array.dm[
                1 : global_const.NJ + 1, 1 : global_const.NK + 1, global_const.NL - 1
            ]
            - global_array.dm[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            / global_array.psg[1 : global_const.NJ + 1, 1 : global_const.NK + 1]
        )


def adiabatic_heating():
    #
    # calculation of adiabatic compression heat
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.comp[j, k, :] = (
                global_const.FKAP
                * (global_array.tw[j, k, :] + global_const.T0)
                * (
                    global_array.uw[j, k, :]
                    * (global_array.ps[j + 1, k] - global_array.ps[j - 1, k])
                    / global_array.dx[k]
                    / 2.0
                    + global_array.vw[j, k, :]
                    * (global_array.ps[j, k + 1] - global_array.ps[j, k - 1])
                    / global_array.dy
                    / 2.0
                )
            )

    global_array.comp[1 : global_const.NJ + 1, 1 : global_const.NK + 1, 0] = (
        global_array.comp[1 : global_const.NJ + 1, 1 : global_const.NK + 1, 0]
        - global_const.FKAP
        * (
            global_array.tw[1 : global_const.NJ + 1, 1 : global_const.NK + 1, 1]
            + global_const.T0
        )
        * global_array.alp[1]
        * global_array.d[1 : global_const.NJ + 1, 1 : global_const.NK + 1, 1]
    )
    for l in range(1, global_const.NL):
        global_array.comp[
            1 : global_const.NJ + 1, 1 : global_const.NK + 1, l
        ] = global_array.comp[
            1 : global_const.NJ + 1, 1 : global_const.NK + 1, l
        ] - global_const.FKAP * (
            global_array.tw[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            + global_const.T0
        ) * (
            global_array.alp[l]
            * global_array.d[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            + (global_array.alp[l] + global_array.alp[l - 1])
            * global_array.dm[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l - 1]
            / global_int.dsig
        )


def div_vert_advection():
    #
    # calculation of divergence of vertical advective flow
    #
    for l in range(0, global_const.NL - 1):
        lp = l + 1
        lm = l - 1
        if lp > global_const.NL:
            lp = global_const.NL
        if lm < 1:
            lm = 1
        global_array.vdivu[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l] = (
            global_array.dsdt[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            * (
                global_array.u[1 : global_const.NJ + 1, 1 : global_const.NK + 1, lp]
                + global_array.u[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            )
            / 2.0
            - global_array.dsdt[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l - 1]
            * (
                global_array.u[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
                + global_array.u[1 : global_const.NJ + 1, 1 : global_const.NK + 1, lm]
            )
            / 2.0
        ) / global_int.dsig
        global_array.vdivv[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l] = (
            global_array.dsdt[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            * (
                global_array.v[1 : global_const.NJ + 1, 1 : global_const.NK + 1, lp]
                + global_array.v[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            )
            / 2.0
            - global_array.dsdt[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l - 1]
            * (
                global_array.v[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
                + global_array.v[1 : global_const.NJ + 1, 1 : global_const.NK + 1, lm]
            )
            / 2.0
        ) / global_int.dsig
        global_array.vdivt[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l] = (
            global_array.dsdt[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            * (
                global_array.t[1 : global_const.NJ + 1, 1 : global_const.NK + 1, lp]
                + global_array.t[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
            )
            / 2.0
            - global_array.dsdt[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l - 1]
            * (
                global_array.t[1 : global_const.NJ + 1, 1 : global_const.NK + 1, l]
                + global_array.t[1 : global_const.NJ + 1, 1 : global_const.NK + 1, lm]
            )
            / 2.0
        ) / global_int.dsig


def summarize_trends():
    #
    # summarize trends
    #
    for k in range(1, global_const.NK + 1):
        for j in range(1, global_const.NJ + 1):
            global_array.ut[j, k, :] = (
                global_array.acx[j, k, :]
                + global_array.apx[j, k, :]
                - global_array.du[j, k, :]
                - global_array.vdivu[j, k, :]
                + global_array.diffu[j, k, :]
            )
            global_array.vt[j, k, :] = (
                global_array.acy[j, k, :]
                + global_array.apy[j, k, :]
                - global_array.dv[j, k, :]
                - global_array.vdivv[j, k, :]
                + global_array.diffv[j, k, :]
            )
            global_array.tt[j, k, :] = (
                -global_array.dvt[j, k, :]
                + global_array.comp[j, k, :]
                - global_array.vdivt[j, k, :]
                + global_array.difft[j, k, :]
            )
            global_array.pst[j, k] = -global_array.dm[j, k, global_const.NL - 1]


#############################################################


def trend():
    true_wind_and_abs_ps()
    geopential()
    zonal_pressure_gradient_force()
    meridional_pressure_gradient_force()
    corioles_and_centrifugal_force()
    div_zonal_impulse()
    div_meridional_impulse()
    div_temperature_flow()
    div_weighted_wind()
    sigma_flow()
    vert_speed_sigma()
    adiabatic_heating()
    div_vert_advection()
    summarize_trends()


if __name__ == "__main__":
    import grid

    grid.grid()

    import init

    init.init_case()

    import boundary_conditions

    boundary_conditions.boundary_conditions()

    trend()
