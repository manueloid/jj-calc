using Plots
include("src/ConstantQuantities.jl")

## Full Hamiltonian
f = ControlParameterFull(0.3, 20)
tspan = range(0.05, 0.5, length=99) |> collect
fid_hess = fidelity_time(f, tspan)
fid = fidelity_time(f, tspan; hessian=false)
plot(tspan, fid_hess, label = "full hess")
plot!(tspan, fid, label = "full no hess")

## Intermediate Hamiltonian
f = ControlParameterInt(0.3, 20)
fid_hess = fidelity_time(f, tspan)
fid = fidelity_time(f, tspan; hessian=false)
plot!(tspan, fid_hess, label = "int hess")
plot!(tspan, fid, label = "int no hess")
