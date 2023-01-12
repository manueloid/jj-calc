using QuantumOptics
using JosephsonSTA

struct ConstantQuantities
    Jz::Operator
    Jx::Operator
    ψ0::Ket
    ψf::Ket
end

"""
constant_quantities(cp::ControlParameter) define the operators and the initial and final ground states of the protocol.
The result is a new ConstantQuantities variable that can be used later.
The idea is to evaluate the fidelity for different final times, so these values will remain the same and there is no need to define them every time.
"""
function ConstantQuantities(cp::ControlParameter)
    Jz = sigmaz(SpinBasis(cp.NParticles / 2)) / 2 |> dense
    Jx = sigmax(SpinBasis(cp.NParticles / 2)) / 2
    ψ0 = eigenstates(2.0 / cp.NParticles * (cp.ω0^2 / 4.0 - 1.0) * Jz^2 - 2.0 * Jx, 1)[2][1]
    ψf = eigenstates(2.0 / cp.NParticles * (cp.ωf^2 / 4.0 - 1.0) * Jz^2 - 2.0 * Jx, 1)[2][1]
    return ConstantQuantities(Jz, Jx, ψ0, ψf)
end

"""
    `rollout(qts::ConstantQuantities)` 
return all the values of the ConstantQuantities variable
"""
function rollout(qts::ConstantQuantities)
    return [getfield(qts, field) for field in fieldnames(typeof(qts))]
end

include("fidelities.jl")
