using Plots
include("src/ConstantQuantities.jl")

function different_λ(λ::Int64, np::Int64, Λf::Float64 )
    ωf = 2.0√(1.0 + Λf)
    ω0 = 2.0
    f = ControlParameterFull(ω0, ωf, 0.3, np )
    i = ControlParameterInt(ω0, ωf, 0.3, np )
    tspan = range(0.05, 0.5, length=99) |> collect
    fid_hess_full = fidelity_time(f, tspan; nlambda = λ)
    fid_full = fidelity_time(f, tspan; hessian=false, nlambda = λ)
    fid_hess_int = fidelity_time(i, tspan; nlambda = λ)
    fid_int = fidelity_time(i, tspan; hessian=false, nlambda = λ)
    fid_sta = fidelity_time_sta(f, tspan)
    plot(tspan, fid_hess_full, label = "Hess full", title = "$λ coeffs $np particles") 
    plot!(tspan, fid_full, label = "No Hess full") 
    plot!(tspan, fid_hess_int, label = " Hess int") 
    plot!(tspan, fid_int, label = "No Hess int") 
    return plot!(tspan, fid_sta, label = "STA") 
end

Λ0 = 0.0
Λf = 200.0
different_λ(4, 30, 200.0)
