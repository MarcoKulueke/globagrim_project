from .variables import global_const, global_int, global_array
from numba import njit # type: ignore

#                             #
# Berechnen der Zeittendenzen #
#                             #

@njit(parallel=True, fastmath=True)
def true_wind_and_abs_ps(NK, NJ, psg, ps, PS0, uw, u, vw, v, tw, t):
    #
    # calculation of true wind component and absolute surface pressure
    #
    for k in range(0, NK + 2):
        for j in range(0, NJ + 2):
            psg[j, k] = ps[j, k] + PS0
            uw[j, k, :] = u[j, k, :] / psg[j, k]
            vw[j, k, :] = v[j, k, :] / psg[j, k]
            tw[j, k, :] = t[j, k, :] / psg[j, k]

@njit(parallel=True, fastmath=True)
def geopential(gp, NL, phis, RD, tw, alp):
    #
    # calculation of geopotential with hydrostatic equation
    #
    gp[:, :, NL - 1] = (
        phis[:, :]
        + RD
        * tw[:, :, NL - 1]
        * alp[NL - 1]
    )
    for l in range(NL - 2, -1, -1):
        gp[:, :, l] = (
            gp[:, :, l + 1]
            + RD
            * (tw[:, :, l] + tw[:, :, l + 1])
            * alp[l]
        )

@njit(parallel=True, fastmath=True)
def zonal_pressure_gradient_force(NK, NJ, apx, RD, tw, T0, ps, dx, psg, gp):
    #
    # zonal pressure gradient force
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            apx[j, k, :] = (
                -RD
                * (tw[j, k, :] + T0)
                * (ps[j + 1, k] - ps[j - 1, k])
                / dx[k]
                / 2.0
                - psg[j, k]
                * (gp[j + 1, k, :] - gp[j - 1, k, :])
                / dx[k]
                / 2.0
            )

@njit(parallel=True, fastmath=True)
def meridional_pressure_gradient_force(NK, NJ, apy, RD, tw, T0, ps, dy, psg, gp):
    #
    # meridional pressure gradient force
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            apy[j, k, :] = (
                -RD
                * (tw[j, k, :] + T0)
                * (ps[j, k + 1] - ps[j, k - 1])
                / dy
                / 2.0
                - psg[j, k]
                * (gp[j, k + 1, :] - gp[j, k - 1, :])
                / dy
                / 2.0
            )

@njit(parallel=True, fastmath=True)
def corioles_and_centrifugal_force(NK, NJ, acx, f, uw, sn, cs, RE, v, acy, u):
    #
    # coriolis and centrifugal force
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            acx[j, k, :] = (
                f[k]
                + uw[j, k, :]
                * sn[k]
                / cs[k]
                / RE
            ) * v[j, k, :]
            acy[j, k, :] = (
                -(
                    f[k]
                    + uw[j, k, :]
                    * sn[k]
                    / cs[k]
                    / RE
                )
                * u[j, k, :]
            )

@njit(parallel=True, fastmath=True)
def div_zonal_impulse(NK, NJ, du, u, uw, dx, v, cs, dy):
    #
    # zonal divergence of momentum
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            du[j, k, :] = (
                (
                    (u[j + 1, k, :] + u[j, k, :])
                    * (uw[j + 1, k, :] + uw[j, k, :])
                    - (u[j, k, :] + u[j - 1, k, :])
                    * (uw[j, k, :] + uw[j - 1, k, :])
                )
                / 4.0
                / dx[k]
            )
            +(
                (
                    v[j, k + 1, :] * cs[k + 1]
                    + v[j, k, :] * cs[k]
                )
                * (uw[j, k + 1, :] + uw[j, k, :])
                - (
                    v[j, k, :] * cs[k]
                    + v[j, k - 1, :] * cs[k - 1]
                )
                * (uw[j, k, :] + uw[j, k - 1, :])
            ) / 4.0 / dy / cs[k]

@njit(parallel=True, fastmath=True)
def div_meridional_impulse(NK, NJ, dv, u, vw, dx, v, cs, dy):
    #
    # meridional divergence of momentum
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            dv[j, k, :] = (
                (
                    (u[j + 1, k, :] + u[j, k, :])
                    * (vw[j + 1, k, :] + vw[j, k, :])
                    - (u[j, k, :] + u[j - 1, k, :])
                    * (vw[j, k, :] + vw[j - 1, k, :])
                )
                / 4.0
                / dx[k]
            )
            +(
                (
                    v[j, k + 1, :] * cs[k + 1]
                    + v[j, k, :] * cs[k]
                )
                * (vw[j, k + 1, :] + vw[j, k, :])
                - (
                    v[j, k, :] * cs[k]
                    + v[j, k - 1, :] * cs[k - 1]
                )
                * (vw[j, k, :] + vw[j, k - 1, :])
            ) / 4.0 / dy / cs[k]

@njit(parallel=True, fastmath=True)
def div_temperature_flow(NK, NJ, dvt, u, tw, dx, v, cs, dy):
    #
    # divergence of temperature flow
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            dvt[j, k, :] = (
                (
                    (u[j + 1, k, :] + u[j, k, :])
                    * (tw[j + 1, k, :] + tw[j, k, :])
                    - (u[j, k, :] + u[j - 1, k, :])
                    * (tw[j, k, :] + tw[j - 1, k, :])
                )
                / 4.0
                / dx[k]
            )
            +(
                (
                    v[j, k + 1, :] * cs[k + 1]
                    + v[j, k, :] * cs[k]
                )
                * (tw[j, k + 1, :] + tw[j, k, :])
                - (
                    v[j, k, :] * cs[k]
                    + v[j, k - 1, :] * cs[k - 1]
                )
                * (tw[j, k, :] + tw[j, k - 1, :])
            ) / 4.0 / dy / cs[k]

@njit(parallel=True, fastmath=True)
def div_weighted_wind(NK, NJ, d, u, dx, v, cs, dy):
    #
    # divergence of mass-wighted wind
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            d[j, k, :] = (
                (u[j + 1, k, :] - u[j - 1, k, :])
                / 2.0
                / dx[k]
                + (
                    v[j, k + 1, :] * cs[k + 1]
                    - v[j, k - 1, :] * cs[k - 1]
                )
                / cs[k]
                / dy
                / 2.0
            )

@njit(parallel=True, fastmath=True)
def sigma_flow(dm, d, dsig, NL):
    #
    # divergence of mass flow between SIGMA=0 and SIGMA=SIGMA[l]
    #
    dm[:, :, 1] = d[:, :, 1] * dsig
    for l in range(2, NL):
        dm[:, :, l] = (
            dm[:, :, l - 1] + d[:, :, l] * dsig
        )

# maybe index error
@njit(parallel=True, fastmath=True)
def vert_speed_sigma(NL, dsdt, NJ, NK, sigma, psg, dm):
    #
    # calculation of vertial velocity in SIGMA system
    #
    for l in range(0, NL - 2):
        dsdt[1 : NJ + 1, 1 : NK + 1, l] = (
            sigma[l]
            / psg[1 : NJ + 1, 1 : NK + 1]
            * dm[
                1 : NJ + 1, 1 : NK + 1, NL - 1
            ]
            - dm[1 : NJ + 1, 1 : NK + 1, l]
            / psg[1 : NJ + 1, 1 : NK + 1]
        )

@njit(parallel=True, fastmath=True)
def adiabatic_heating(NK, NJ, comp, FKAP, tw, T0, uw, ps, dx, vw, dy, alp, d, NL, dm, dsig):
    #
    # calculation of adiabatic compression heat
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            comp[j, k, :] = (
                FKAP
                * (tw[j, k, :] + T0)
                * (
                    uw[j, k, :]
                    * (ps[j + 1, k] - ps[j - 1, k])
                    / dx[k]
                    / 2.0
                    + vw[j, k, :]
                    * (ps[j, k + 1] - ps[j, k - 1])
                    / dy
                    / 2.0
                )
            )

    comp[1 : NJ + 1, 1 : NK + 1, 0] = (
        comp[1 : NJ + 1, 1 : NK + 1, 0]
        - FKAP
        * (
            tw[1 : NJ + 1, 1 : NK + 1, 1]
            + T0
        )
        * alp[1]
        * d[1 : NJ + 1, 1 : NK + 1, 1]
    )
    for l in range(1, NL):
        comp[
            1 : NJ + 1, 1 : NK + 1, l
        ] = comp[
            1 : NJ + 1, 1 : NK + 1, l
        ] - FKAP * (
            tw[1 : NJ + 1, 1 : NK + 1, l]
            + T0
        ) * (
            alp[l]
            * d[1 : NJ + 1, 1 : NK + 1, l]
            + (alp[l] + alp[l - 1])
            * dm[1 : NJ + 1, 1 : NK + 1, l - 1]
            / dsig
        )

@njit(parallel=True, fastmath=True)
def div_vert_advection(NL, vdivu, NJ, NK, dsdt, u, vdivv, v, vdivt, t, dsig):
    #
    # calculation of divergence of vertical advective flow
    #
    for l in range(0, NL - 1):
        lp = l + 1
        lm = l - 1
        if lp > NL:
            lp = NL
        if lm < 1:
            lm = 1
        vdivu[1 : NJ + 1, 1 : NK + 1, l] = (
            dsdt[1 : NJ + 1, 1 : NK + 1, l]
            * (
                u[1 : NJ + 1, 1 : NK + 1, lp]
                + u[1 : NJ + 1, 1 : NK + 1, l]
            )
            / 2.0
            - dsdt[1 : NJ + 1, 1 : NK + 1, l - 1]
            * (
                u[1 : NJ + 1, 1 : NK + 1, l]
                + u[1 : NJ + 1, 1 : NK + 1, lm]
            )
            / 2.0
        ) / dsig
        vdivv[1 : NJ + 1, 1 : NK + 1, l] = (
            dsdt[1 : NJ + 1, 1 : NK + 1, l]
            * (
                v[1 : NJ + 1, 1 : NK + 1, lp]
                + v[1 : NJ + 1, 1 : NK + 1, l]
            )
            / 2.0
            - dsdt[1 : NJ + 1, 1 : NK + 1, l - 1]
            * (
                v[1 : NJ + 1, 1 : NK + 1, l]
                + v[1 : NJ + 1, 1 : NK + 1, lm]
            )
            / 2.0
        ) / dsig
        vdivt[1 : NJ + 1, 1 : NK + 1, l] = (
            dsdt[1 : NJ + 1, 1 : NK + 1, l]
            * (
                t[1 : NJ + 1, 1 : NK + 1, lp]
                + t[1 : NJ + 1, 1 : NK + 1, l]
            )
            / 2.0
            - dsdt[1 : NJ + 1, 1 : NK + 1, l - 1]
            * (
                t[1 : NJ + 1, 1 : NK + 1, l]
                + t[1 : NJ + 1, 1 : NK + 1, lm]
            )
            / 2.0
        ) / dsig

@njit(parallel=True, fastmath=True)
def summarize_trends(NK, NJ, ut, acx, apx, du, vdivu, diffu, vt, acy, apy, dv, vdivv, diffv, tt, dvt, comp, vdivt, difft, pst, dm, NL):
    #
    # summarize trends
    #
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            ut[j, k, :] = (
                acx[j, k, :]
                + apx[j, k, :]
                - du[j, k, :]
                - vdivu[j, k, :]
                + diffu[j, k, :]
            )
            vt[j, k, :] = (
                acy[j, k, :]
                + apy[j, k, :]
                - dv[j, k, :]
                - vdivv[j, k, :]
                + diffv[j, k, :]
            )
            tt[j, k, :] = (
                -dvt[j, k, :]
                + comp[j, k, :]
                - vdivt[j, k, :]
                + difft[j, k, :]
            )
            pst[j, k] = -dm[j, k, NL - 1]


#############################################################

def trend():
    true_wind_and_abs_ps(global_const.NK, global_const.NJ, global_array.psg, global_array.ps, global_const.PS0, global_array.uw, global_array.u, global_array.vw, global_array.v, global_array.tw, global_array.t)
    geopential(global_array.gp, global_const.NL, global_array.phis, global_const.RD, global_array.tw, global_array.alp)
    zonal_pressure_gradient_force(global_const.NK, global_const.NJ, global_array.apx, global_const.RD, global_array.tw, global_const.T0, global_array.ps, global_array.dx, global_array.psg, global_array.gp)
    meridional_pressure_gradient_force(global_const.NK, global_const.NJ, global_array.apy, global_const.RD, global_array.tw, global_const.T0, global_array.ps, global_int.dy, global_array.psg, global_array.gp)
    corioles_and_centrifugal_force(global_const.NK, global_const.NJ, global_array.acx, global_array.f, global_array.uw, global_array.sn, global_array.cs, global_const.RE, global_array.v, global_array.acy, global_array.u)
    div_zonal_impulse(global_const.NK, global_const.NJ, global_array.du, global_array.u, global_array.uw, global_array.dx, global_array.v, global_array.cs, global_int.dy)
    div_meridional_impulse(global_const.NK, global_const.NJ, global_array.dv, global_array.u, global_array.vw, global_array.dx, global_array.v, global_array.cs, global_int.dy)
    div_temperature_flow(global_const.NK, global_const.NJ, global_array.dvt, global_array.u, global_array.tw, global_array.dx, global_array.v, global_array.cs, global_int.dy)
    div_weighted_wind(global_const.NK, global_const.NJ, global_array.d, global_array.u, global_array.dx, global_array.v, global_array.cs, global_int.dy)
    sigma_flow(global_array.dm, global_array.d, global_int.dsig, global_const.NL)
    vert_speed_sigma(global_const.NL, global_array.dsdt, global_const.NJ, global_const.NK, global_array.sigma, global_array.psg, global_array.dm)
    adiabatic_heating(global_const.NK, global_const.NJ, global_array.comp, global_const.FKAP, global_array.tw, global_const.T0, global_array.uw, global_array.ps, global_array.dx, global_array.vw, global_int.dy, global_array.alp, global_array.d, global_const.NL, global_array.dm, global_int.dsig)
    div_vert_advection(global_const.NL, global_array.vdivu, global_const.NJ, global_const.NK, global_array.dsdt, global_array.u, global_array.vdivv, global_array.v, global_array.vdivt, global_array.t, global_int.dsig)
    summarize_trends(global_const.NK, global_const.NJ, global_array.ut, global_array.acx, global_array.apx, global_array.du, global_array.vdivu, global_array.diffu, global_array.vt, global_array.acy, global_array.apy, global_array.dv, global_array.vdivv, global_array.diffv, global_array.tt, global_array.dvt, global_array.comp, global_array.vdivt, global_array.difft, global_array.pst, global_array.dm, global_const.NL)


if __name__ == "__main__":
    import grid

    grid.grid()

    import init

    init.init_case()

    import boundary_conditions

    boundary_conditions.boundary_conditions()

    trend()
