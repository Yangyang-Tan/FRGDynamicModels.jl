using FRGDynamicModels
using Plots

sol_lpap0 = solver_modelA(;Λ=430.0, kIR=1.0, n=4.0, Tem=141.633955
)

let rho = sol_lpap0[2], k = 1.0, Tem = 141.633955, d = 3.0
    rhobar = rho ./ (k^(d - 2) * Tem)
    t = log(k / 430.0)
    plot(
        rho[1:end],
        sol_lpap0[1](t)[1:end];
        label="T=$Tem,k=$k",
        seriestype=:scatter,
        xlabel="ρ",
        ylabel="U'",
    )
end