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


ω_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = control_ω(t, cp) - correction_poly(t, cp, corrections::Vector{Float64})

function fidelity(cp::ControlParameter; nlambda::Int64=5, maxbra::Int64=4)
    corrs = corrections(cp; nlambda=nlambda, maxbra=maxbra)
    qts = ConstantQuantities(cp)
    fidelity_col = zeros
    ω(t) = ω_esta(t, cp, corrs)
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
        return (2.0 / cp.NParticles * (ω(t) / 4.0 - 1.0) * qts.Jz^2 - 2.0 * qts.Jx)
    end
    fidelity(t, psi) = abs2.(dagger(qts.ψf) * psi)
    timeevolution.schroedinger_dynamic([0.0, cp.final_time], qts.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
end

function fidelity(cp::ControlParameter, corrs::Vector{Float64}; nlambda::Int64=5, maxbra::Int64=4)
    qts = ConstantQuantities(cp)
    ω(t) = ω_esta(t, cp, corrs)
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
        return (2.0 / cp.NParticles * (ω(t) / 4.0 - 1.0) * qts.Jz^2 - 2.0 * qts.Jx)
    end
    fidelity(t, psi) = abs2.(dagger(qts.ψf) * psi)
    timeevolution.schroedinger_dynamic([0.0, cp.final_time], qts.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
end

function fidelity_time(cp::ControlParameter, tarr::Vector{Float64}; hessian::Bool=true, nlambda::Int64=5, maxbra::Int64=4)
    qts = ConstantQuantities(cp)
    fidelity_array = zeros(length(tarr))
    Threads.@threads for (index, tf) in enumerate(tarr) |> collect
        # for (index, tf) in enumerate(tarr) |> collect
        cparam = cp_time(cp, tf)
        if hessian == true
            corrs = corrections_hess(cparam, nlambda=nlambda)
        else
            corrs = corrections(cparam, nlambda=nlambda, maxbra=maxbra)
        end
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

function fidelity_search(cp::ControlParameter, epsilons::Vector{Float64}; nlambda::Int64=5, maxbra::Int64=4)
    corrs = corrections(cp, nlambda=nlambda, maxbra=maxbra)
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

