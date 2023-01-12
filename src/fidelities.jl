ω_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = control_ω(t, cp) - correction_poly(t, cp, corrections::Vector{Float64})

function fidelity_time(cp::ControlParameter, tarr::Vector{Float64}; hessian::Bool=true, nlambda::Int64=5, maxbra::Int64=4)
    qts = ConstantQuantities(cp)
    fidelity_array = zeros(length(tarr))
    Threads.@threads for (index, tf) in enumerate(tarr) |> collect
        # for (index, tf) in enumerate(tarr) |> collect
        cparam = cp_time(cp, tf)
        corrs = corrections(cparam; hessian=hessian, nlambda=nlambda, maxbra=maxbra)
        ω(t) = ω_esta(t, cparam, corrs)
        function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
            return (2.0 / cparam.NParticles * (ω(t) / 4.0 - 1.0) * qts.Jz^2 - 2.0 * qts.Jx)
        end
        fidelity(t, psi) = abs2.(dagger(qts.ψf) * psi)
        fidelity_array[index] = timeevolution.schroedinger_dynamic([0.0, tf], qts.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
        println("Calculating fidelity for final time $tf ")
    end
    return fidelity_array
end

function fidelity_time_sta(cp::ControlParameter, tarr::Vector{Float64}; hessian::Bool=true, nlambda::Int64=5, maxbra::Int64=4)
    qts = ConstantQuantities(cp)
    fidelity_array = zeros(length(tarr))
    Threads.@threads for (index, tf) in enumerate(tarr) |> collect
        # for (index, tf) in enumerate(tarr) |> collect
        cparam = cp_time(cp, tf)
        ω(t) = control_ω(t, cparam)
        function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
            return (2.0 / cparam.NParticles * (ω(t) / 4.0 - 1.0) * qts.Jz^2 - 2.0 * qts.Jx)
        end
        fidelity(t, psi) = abs2.(dagger(qts.ψf) * psi)
        fidelity_array[index] = timeevolution.schroedinger_dynamic([0.0, tf], qts.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
        println("Calculating fidelity for final time $tf ")
    end
    return fidelity_array
end

function fidelity_search(cp::ControlParameter, epsilons::Vector{Float64}; hessian::Bool=true, nlambda::Int64=5, maxbra::Int64=4)
    corrs = corrections(cp; hessian=hessian, nlambda=nlambda, maxbra=maxbra)
    qts = ConstantQuantities(cp)
    fidelity_array = zeros(length(epsilons))
    Threads.@threads for (index, ϵ) in enumerate(epsilons) |> collect
        # for (index, ϵ) in enumerate(epsilons) |> collect
        ω(t) = ω_esta(t, cp, ϵ * corrs)
        function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
            return (2.0 / cp.NParticles * (ω(t) / 4.0 - 1.0) * qts.Jz^2 - 2.0 * qts.Jx)
        end
        fidelity(t, psi) = abs2.(dagger(qts.ψf) * psi)
        println("Calculating fidelity for ϵ $ϵ")
        fidelity_array[index] = timeevolution.schroedinger_dynamic([0.0, cp.final_time], qts.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
    end
    return fidelity_array
end

"""
fidelity(cp::ControlParameter, qts::ConstantQuantities, ω::Function) return the fidelity of a process with desired boundary conditions, constant quantities and the control function ω.
Useful to include it in other functions to iteratively calculate the fidelity
"""
function fidelity(cp::ControlParameter, qts::ConstantQuantities, ω::Function; hessian::Bool=true, nlambda::Int64=5, maxbra::Int64=4)
    Jz, Jx, ψ0, ψf = rollout(qts)
    h = 2.0 / cp.NParticles
    tf = cp.final_time
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
        return (h * (ω(t) / 4.0 - 1.0) * Jz^2 - 2.0 * Jx)
    end
    fidelity(t, psi) = abs2.(dagger(ψf) * psi)
    return timeevolution.schroedinger_dynamic([0.0, tf], ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
end
