# Oneliner functions to produce the control function of the system, both for Λ and ω
Λ_sta(t::Float64, cp::ControlParameter) = control_ω(t, cp) * 0.25 - 1.0
ω_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = control_ω(t, cp) - correction_poly(t, cp, corrections::Vector{Float64})
function ω_esta(t::Float64, cp::ControlParameter)
	corrs = corrections(cp)
	return ω_esta(t, cp, corrs)
end
Λ_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = ω_esta(t, cp, corrections) * 0.25 - 1.0
Λ_esta(t::Float64, cp::ControlParameter) = ω_esta(t, cp) * 0.25 - 1.0

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
fidelities(cp::ControlParameter, final_times)
Return the fidelity of a process with desired boundary conditions and constant quantities for a range of final times, for both the eSTA and STA protocol.
The output is a tuple with the fidelities for the eSTA protocol in the first position and the fidelities for the STA protocol in the second position.
"""
function fidelities(cp::ControlParameter, final_times::Vector{Float64})
		qts = ConstantQuantities(cp) # Constant quantities not depending on the final time
		fidelities_esta = zeros(length(final_times)) # Array to store the fidelities for the eSTA protocol
		fidelities_sta = zeros(length(final_times)) # Array to store the fidelities for the STA protocol
		Threads.@threads for (index, tf) in enumerate(final_times)
			local cp = cp_time(cp, tf) # Control parameter with the new final time
			corrs = corrections(cp) # Corrections for the eSTA protocol for the new final time
			esta(t) = Λ_esta(t, cp, corrs) # eSTA control function for the new final time and corrections
			sta(t) = Λ_sta(t, cp) # STA control function for the new final time 
			fidelities_esta[index] = fidelity(cp, qts, esta) # Fidelity for the eSTA protocol
			println("Fidelity eSTA: ", fidelities_esta[index])
			fidelities_sta[index] = fidelity(cp, qts, sta) # Fidelity for the STA protocol
			println("Fidelity STA: ", fidelities_sta[index])
		end
		return fidelities_esta, fidelities_sta # Return the fidelities for both protocols as a tuple
end
