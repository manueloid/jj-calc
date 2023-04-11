include("../src/ConstantQuantities.jl")
# I need to backpropagate the wavefunction from the final state to the initial one
# First, I will define a function that will perform the time reversal of a generic function
"""
time_reverse(t::Float64, f::Function, cp::ControlParameter) -> f(cp.final_time - t) 
This function performs the time reversal of a generic function f(t,cp) with respect to the final time of the whole protocol.
The idea is to reverse the functions that are used to define the control parameters of the protocol, i.e. ω or Λ
"""
function time_reverse(t::Float64, f::Function, cp::ControlParameter)
    tf = cp.final_time
    return f(tf - t, cp)
end

cp = ControlParameterFull(0.2, 10)
qts = ConstantQuantities(cp)
"""
Hamiltonian(Jx, Jz, ω::Function, h::Float64) -> H(t, psi)
This function returns the time dependent Hamiltonian of the system.
Here is the breakdown of the arguments:
	- Jx, Jz: Angular momentum of the system
	- ω: control function
	- h: system size parameter
"""
function Hamiltonian(Jx, Jz, ω::Function, h::Float64)
    function H(t, psi)
        return (h * ω(t) * Jz^2 - 2.0 * Jx)
    end
    return H
end
# Multiple dispatch for the Hamiltonian function
Hamiltonian(Jx, Jz, ω::Function, NParticles::Int64) = Hamiltonian(Jx, Jz, ω, 2.0 / NParticles)
Hamiltonian(qts::ConstantQuantities, ω::Function, h::Float64) = Hamiltonian(qts.Jx, qts.Jz, ω, h)
"""
evolution(ψ0, H, tspan) -> ψ(t) 
Time evolution of a system starting from the initial state ψ0, under the Hamiltonian H, for a time span tspan.
"""
function evolution(ψ0, H, tspan)
    return timeevolution.schroedinger_dynamic(tspan, ψ0, H)
end

"""
evolution(cp::ControlParameter, qts::ConstantQuantities, ω::Function; steps=1000) -> ψ(t)
This function performs the time evolution of the system, where the conditions are defined 
by the control parameter cp, the constant quantities qts and the control function ω.
"""
function evolution(cp::ControlParameter, qts::ConstantQuantities, ω::Function; steps=1000)
    h = 2.0 / cp.NParticles
    tf = cp.final_time
    time = range(0.0, tf, length=steps) |> collect
    ψ0 = qts.ψ0
    H = Hamiltonian(qts, ω, h)
    return evolution(ψ0, H, time)
end

# Now I plot the evolution of the wavefunction
cp = ControlParameterFull(0.2, 10)
qts = ConstantQuantities(cp)
l(t) = Λ_esta(t, cp)
ψ = evolution(cp, qts, l)
ψ0 = qts.ψ0
