function flowModelA(dVdtall, Vpall, p, t)
    (Λ, νd, T, d, n, rho_grid) = p
    k = Λ * exp(t)
    Vp = @view Vpall[1:(end)]
    dVdt = @view dVdtall[1:(end)]
    Vpfun = Spline1D(rho_grid, Vp; k=3, bc="nearest")
    rho0 = if Vpall[1] < 0.0
        roots(Vpfun)[1]
    else
        0.0
    end
    Vpp0 = derivative(Vpfun, rho0)
    Vpp = derivative(Vpfun, rho_grid)
    Vppp = derivative(Vpfun, rho_grid, 2)
    η = (4 * (2 + d) * k^(-2 + d) * T * Vpp0^2 * νd * rho0) / (k^2 + 2 * Vpp0 * rho0)^2
    @. dVdt =
        η * Vp + η * rho_grid * Vpp - (
            k^(2 + d) *
            T *
            (2 + d - η) *
            νd *
            (
                ((-1 + n) * Vpp) / (k^(2) + Vp)^2 +
                (3 * Vpp + 2 * Vppp * rho_grid) / (k^(2) + Vp + 2 * Vpp * rho_grid)^2
            )
        )
    return nothing
end