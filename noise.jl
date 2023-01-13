using Plots
include("./src/ConstantQuantities.jl")
## What I want is a fidelity that takes a control parameter and the control functions and spits out the  . 
##
f = ControlParameterFull(0.2, 40)
tspan = range(0.1, 0.5, length=99) |> collect

esta_rob = robustness_time(f, tspan)
sta_rob = robustness_time(f, tspan; esta=false)
##
plot(tspan, esta_rob)
plot!(tspan, sta_rob)

##
fid_esta_hess = fidelity_time(f, tspan)
fid_esta_no_hess = fidelity_time(f, tspan; hessian=false)

plot(tspan, fid_esta_hess)
plot!(tspan, fid_esta_no_hess)
