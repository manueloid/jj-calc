using Plots
include("./src/ConstantQuantities.jl")
# What I want is a fidelity that takes a control parameter and the control functions and spits out the  . 
f = ControlParameterFull(0.5, 30)
qts = ConstantQuantities(f)
corrs = corrections(f; hessian=false)
ω(t) = ω_esta(t, f, corrs)

using FiniteDifferences
"""
`robustness(cp::ControlParameter, qts::ConstantQuantities, ω::Function; order::Int64=5)`
return the derivative of the fidelity given an error of the form 1 + δ for δ = 0
# Arguments 
- cp 
- qts 
- ω this must be a function of only one parameter (ideally the time) and it has to be built outside the function declaration 
# Kwargs 
- order is the order of the finite difference method I need to use
"""
function robustness(cp::ControlParameter, qts::ConstantQuantities, ω::Function; order::Int64=5)
    function error(δ)
        ω_err(t) = (1 + δ) * ω(t)
        return fidelity(cp, qts, ω_err)
    end
    m = central_fdm(order, 1)
    return m(error, 0.0)
end

function robustness_time(cp::ControlParameter, tarr::Vector{Float64}; order::Int64=5)
    qts = ConstantQuantities(cp)
    ω(t) = ω_esta(t, f, corrs)
end
