include("./src/ConstantQuantities.jl")
f = ControlParameterFull(0.2, 40)
i = ControlParameterInt(0.2, 40)
tspan = range(0.1, 0.5, length=9) |> collect
nlambda = 3
robustness_time(f, tspan; esta=false) # Robustness of eSTA with the full Hamiltonian and hessian
##
using Plots
function robustness_all(np::Int64, tspan::Vector{Float64}, Λf::Float64=50.0; nlambda::Int64=5, ω0::Float64=2.0)
        f = ControlParameterFull(ω0, 2.0√(Λf + 1.0), 0.2, np)
        i = ControlParameterInt(ω0, 2.0√(Λf + 1.0), 0.2, np)
        whole = [
                robustness_time(f, tspan; nlambda=nlambda), # Robustness of eSTA with the full Hamiltonian and hessian
                robustness_time(f, tspan; nlambda=nlambda, hessian=false), # Robustness of eSTA with the full Hamilhonian and NO hessian
                robustness_time(i, tspan; nlambda=nlambda), # Robustness of eSTA with the intermediate Hamiltonian and hessian
                robustness_time(i, tspan; nlambda=nlambda, hessian=false), # Robustness of eSTA with the intermediate Hamilhonian and NO hessian
                robustness_time(f, tspan; esta=false) # STA robustness ( as esta is set to false)
        ]
        return whole
end
test = robustness_all(50, tspan; nlambda=2)
##
plot(tspan, test, label=["full" "full noh" "int" "int noh" "sta"])


