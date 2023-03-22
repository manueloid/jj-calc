#= Here I need to evaluate the number squeezing parameter, i.e. ξ² = ΔJz² \ (N\4) 
What I want is a simulation of the Schrödinger equation that given the final time,
the boundary conditions and the number of particles, calculates the corrections and then saved the desired quantity.
I will only simulate the full Hessian and the STA Protocol to compare the two =#
using JosephsonSTA
include("src/ConstantQuantities.jl")
Λ_sta(t::Float64, cp::ControlParameter) = control_ω(t, cp) * 0.25 - 1.0
ω_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = control_ω(t, cp) - correction_poly(t, cp, corrections::Vector{Float64})
Λ_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = ω_esta(t, cp, corrections) * 0.25 - 1.0

f = ControlParameterFull(0.3, 30)
q = ConstantQuantities(f)

ξ(t, psi) = dagger(psi) * (q.Jz)^2 * psi - dagger(psi) * (q.Jz) * psi |> real

```
function squeezing(cp::ControlParameter,           qts::ConstantQuantities, f::Function)
    h = 2.0/cp.N
    tf = cp.final_time
# to finish
end
end