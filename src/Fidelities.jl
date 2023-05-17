# Oneliner functions to produce the control function of the system, both for Λ and ω
Λ_sta(t::Float64, cp::ControlParameter) = control_ω(t, cp) * 0.25 - 1.0
ω_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = control_ω(t, cp) - correction_poly(t, cp, corrections::Vector{Float64})
function ω_esta(t::Float64, cp::ControlParameter)
    corrs = corrections(cp)
    return ω_esta(t, cp, corrs)
end
Λ_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = ω_esta(t, cp, corrections) * 0.25 - 1.0
Λ_esta(t::Float64, cp::ControlParameter) = ω_esta(t, cp) * 0.25 - 1.0
Λ_ad(t::Float64, cp::ControlParameter) = control_ad(t, cp) * 0.25 - 1.0

# Function that will be used by other functions
"""
fidelity(cp::ControlParameter, qts::ConstantQuantities, ω::Function) return the fidelity of a process with desired boundary conditions, constant quantities and the control function ω.
Useful to include it in other functions to iteratively calculate the fidelity
"""
function fidelity(cp::ControlParameter, qts::ConstantQuantities, Λ::Function)
    h = 2.0 / cp.NParticles
    tf = cp.final_time
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
        return (h * Λ(t) * qts.Jz^2 - 2.0 * qts.Jx)
    end
    fidelity(t, psi) = abs2.(dagger(qts.ψf) * psi)
    return timeevolution.schroedinger_dynamic([0.0, tf], qts.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
end
"""
fidelity_evo(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, npoints=100)
return the time evolution of the fidelity of the state with the target state, given the control parameter, the constant quantities, the control function and the number of points to evaluate the fidelity.
The initial time is 0 and the final time is the final time of the control parameter.
"""
function fidelity(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, npoints)
    h = 2.0 / cp.NParticles
    tf = cp.final_time
    tspan = range(0.0, tf, length=npoints)
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian	
        return (h * Λ(t) * qts.Jz^2 - 2.0 * qts.Jx)
    end
    f(t, psi) =  dagger(qts.ψf) *  psi |> abs2
    return timeevolution.schroedinger_dynamic(tspan, qts.ψ0, H_eSTA; fout=f)[2]
end
