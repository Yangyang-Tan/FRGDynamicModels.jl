using Sundials,ProgressLogging
using DifferentialEquations
using SpecialFunctions
using Dierckx
using Plots
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
function solver_modelA(;
    d::T=3.0,
    n::T=4.0,
    ηpi::T=0.0,
    ησ::T=0.0,
    Zpi::T=1.0,
    Zσ::T=1.0,
    Tem::T=1.0,
    rhomin::T=0.1,
    rhomax::T=0.5,
    ngrid::Int=100,
    Λ::T=1000.0,
    kIR::T=1.0,
    upprho0ini::T=4.614236790859641,
    rho0ini::T=0.07147268020712277,
) where {T}
    νd = 1 / (2^(d + 1) * pi^(d / 2) * gamma(2 + d / 2))
    σ0 = 0.0:0.01:1.0
    rho0 = 0.5 * σ0 .^ 2
    Tc = 150.0
    rho_grid = rho0 * (1000^(d - 2) * Tc)
    u0ini =
        35.0 .* (rho_grid .- 0.358 * rho0ini * (1000^(d - 2) * Tc))
    u0 = u0ini
    p = (Λ, νd, Tem, d, n, rho_grid)
    tspan = (0.0, log(kIR / Λ))
    prob = ODEProblem(flowModelA, u0, tspan, p)
    sol = solve(
        prob,
        CVODE_BDF(; method=:Newton, linear_solver=:Dense),
        ;
        reltol=1e-8,
        dtmax=0.0005,
        progress=true,
        abstol=1e-8,
        saveat=-0.05,
        maxiters=2 * 10^4,
    )
    return [sol, rho_grid, Tem]
end

sol_lpap0 = solver_modelA_T_LPAP(;
    upprho0ini=upprho0ini0, rho0ini=rho0ini0, Λ=430.0, kIR=10.0, n=1.0, Tem=477.697135
)

let rho = sol_lpap0[2], k = 10.0, Tem = 477.697135, d = 3.0
    rhobar = rho ./ (k^(d - 2) * Tem)
    t = log(k / 430.0)
    plot(
        rho[1:5],
        sol_lpap0[1](t)[1:5];
        label="T=150,k=0.1",
        seriestype=:scatter,
        xlabel="ρ",
        ylabel="U'",
    )
end