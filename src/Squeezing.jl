# squeezing return the squeezing parameter ξ with respect of time. ξ(t)
"""
squeezing(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, npoints=100) -> (tspan, ξ)
This function returns the squeezing parameter ξ² with respect of time as a tuple (tspan, ξ²)
"""
function squeezing(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, npoints=100)
    h = 2.0 / cp.NParticles
    tf = cp.final_time
    tspan = range(0.0, tf, length=npoints)
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian	
        return (h * Λ(t) * qts.Jz^2 - 2.0 * qts.Jx)
    end
    ΔJ(t, psi) = dagger(psi) * (qts.Jz)^2 * psi - (dagger(psi) * (qts.Jz) * psi)^2 |> real
    ξ = timeevolution.schroedinger_dynamic(tspan, qts.ψ0, H_eSTA; fout=ΔJ)[2] .* 2.0h
    return tspan, ξ
end
"""
squeezings(cp::ControlParameter, final_times)
This function returns the squeezing parameter ξ² for different final times for the STA protocol and sthe eSTA one
"""
function squeezings(cp::ControlParameter, final_times)
    squeezings_esta = zeros(length(final_times))
    squeezings_sta = zeros(length(final_times))
    qts = ConstantQuantities(cp)
    Threads.@threads for i in 1:length(final_times)
        cparam = cp_time(cp, final_times[i])
        corrs = corrections(cparam)
        esta(t) = Λ_esta(t, cparam, corrs)
        sta(t) = Λ_sta(t, cparam)
        squeezings_esta[i] = squeezing(cparam, qts, esta, 2)[2][2]
        println("Squeezing for eSTA: ", squeezings_esta[i], " for final time: ", final_times[i])
        squeezings_sta[i] = squeezing(cparam, qts, sta, 2)[2][2]
        println("Squeezing for STA: ", squeezings_sta[i], " for final time: ", final_times[i])
    end
    return squeezings_esta, squeezings_sta
end

"""
squeezing(cp::ControlParameter) -> Float64 
This function returns the squeezing parameter ξ² for the ideal case.
"""
function squeezing(cp::ControlParameter)
    qts = ConstantQuantities(cp)
    ΔJ = dagger(qts.ψf) * (qts.Jz)^2 * qts.ψf - (dagger(qts.ψf) * (qts.Jz) * qts.ψf)^2 |> real
    h = 2.0 / cp.NParticles
    return ΔJ * 2.0h
end
