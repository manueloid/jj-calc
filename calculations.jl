using JosephsonSTA
using QuantumOptics
include("src/ConstantQuantities.jl")

tspan = [range(0.05, 0.5, length=99);]
np = 50
f = ControlParameterFull(np, tspan[1])
i = ControlParameterInt(np, tspan[1])

fid_time_full = fidelity_time(f, tspan)
fid_time_int = fidelity_time(i, tspan)

using Plots
plot(tspan, fid_time_full)
plot!(tspan, fid_time_int)
