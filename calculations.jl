using Plots
include("./src/ConstantQuantities.jl")
function fidelity_all(np::Int64, tspan::Vector{Float64}, Λf::Float64 = 50.0; nlambda::Int64 = 5, ω0::Float64=2.0   )
    f = ControlParameterFull(ω0, 2.0√(Λf + 1.0), .2, np )
    i = ControlParameterInt(ω0, 2.0√(Λf + 1.0), .2, np )
    whole = [
        fidelity_time(f, tspan; nlambda = nlambda), # Fidelity of eSTA with the full Hamiltonian and hessian
        fidelity_time(f, tspan;  nlambda = nlambda, hessian = false), # Fidelity of eSTA with the full Hamilhonian and NO hessian
        fidelity_time(i, tspan;  nlambda = nlambda), # Fidelity of eSTA with the intermediate Hamiltonian and hessian
        fidelity_time(i, tspan;  nlambda = nlambda, hessian = false), # Fidelity of eSTA with the intermediate Hamilhonian and NO hessian
        fidelity_time(f, tspan; esta = false) # STA fidelity ( as esta is set to false)
    ]
    return whole
end
##
tspan = range(.04, 1.0, length = 10)|> collect
test = fidelity_all(50, tspan)
plot(tspan, test)
