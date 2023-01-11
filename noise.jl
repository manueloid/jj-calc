using JosephsonSTA
using QuantumOptics
include("src/ConstantQuantities.jl")

function fidelity(cp::ControlParameter, qts::ConstantQuantities, ω::Function)
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
        return (2.0 / cp.NParticles * (ω(t) / 4.0 - 1.0) * qts.Jz^2 - 2.0 * qts.Jx)
    end
    fidelity(t, psi) = abs2.(dagger(qts.ψf) * psi)
    timeevolution.schroedinger_dynamic([0.0, cp.final_time], qts.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
end

##
cparam = ControlParameterFull()
ctsq = ConstantQuantities(cparam)
corrs = corrections(cparam)
ω(t) = ω_esta(t, cparam, corrs)
ω_noise(t) = ω(t) + 1
fidelity(cparam, ctsq, ω_noise)
