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