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
function fidelity(cp::ControlParameter, qts::ConstantQuantities, Λ::Function; hessian::Bool=true, nlambda::Int64=5, maxbra::Int64=4)
        h = 2.0 / cp.NParticles
        tf = cp.final_time
        function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
                return (h * Λ(t) * qts.Jz^2 - 2.0 * qts.Jx)
        end
        fidelity(t, psi) = abs2.(dagger(qts.ψf) * psi)
        return timeevolution.schroedinger_dynamic([0.0, tf], qts.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
end

function fidelity_time_esta(cp::ControlParameter, tarr::Vector{Float64}; hessian::Bool=true, nlambda::Int64=5, maxbra::Int64=4)
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

function fidelity_time(cp::ControlParameter, tarr::Vector{Float64}; esta::Bool=true, hessian::Bool=true, nlambda::Int64=5, maxbra::Int64=4)
        if esta == true
                return fidelity_time_esta(cp, tarr; hessian=hessian, nlambda=nlambda, maxbra=maxbra)
        else
                return fidelity_time_sta(cp, tarr)
        end
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

using FiniteDifferences
"""
`robustness(cp::ControlParameter, qts::ConstantQuantities, ω::Function; order::Int64=5)`
return the derivative of the fidelity given an error in the modulation of the control parameter of the form Λδ(t) = (1 + δ)Λ(t)
# Arguments 
- `cp` 
- `qts` 
- `Λ` must be a function of only one parameter (ideally the time) and it has to be built outside the function declaration 
# Kwargs 
- `order` is the order of the finite difference method I need to use
"""
function robustness(cp::ControlParameter, qts::ConstantQuantities, Λ::Function; order::Int64=5)
        function error(δ)
                Λ_err(t) = (1 + δ) * Λ(t)
                return fidelity(cp, qts, Λ_err)
        end
        m = central_fdm(order, 1)
        return m(error, 0.0) |> abs
end
"""
`robustness_timenoise(cp::ControlParameter, qts::ConstantQuantities, ω::Function; order::Int64=5)`
return the derivative of the fidelity given an offset in the time of the control parameter of the form Λδ(t) = Λ(t + δ)
# Arguments 
- `cp` 
- `qts` 
- `Λ` must be a function of only one parameter (ideally the time) and it has to be built outside the function declaration 
# Kwargs 
- `order` is the order of the finite difference method I need to use
"""
function robustness_timenoise(cp::ControlParameter, qts::ConstantQuantities, Λ::Function; order::Int64=5)
        function error(δ)
                Λ_err(t) = Λ(δ + t)
                return fidelity(cp, qts, Λ_err)
        end
        m = central_fdm(order, 1)
        return m(error, 0.0) |> abs
end

function robustness_time_esta(cp::ControlParameter, tarr::Vector{Float64}; nlambda::Int64=5, hessian::Bool=true)
        qts = ConstantQuantities(cp)
        robustness_arr = zeros(length(tarr))
        Threads.@threads for (index, tf) in enumerate(tarr) |> collect
                cparam = cp_time(cp, tf)
                corrs = corrections(cparam; hessian=hessian, nlambda=nlambda)
                esta(t) = Λ_esta(t, cparam, corrs)
                println("Calculating robustness for time $tf")
                robustness_arr[index] = robustness(cparam, qts, esta)
        end
        return robustness_arr
end

function robustness_time_sta(cp::ControlParameter, tarr::Vector{Float64})
        qts = ConstantQuantities(cp)
        robustness_arr = zeros(length(tarr))
        Threads.@threads for (index, tf) in enumerate(tarr) |> collect
                cparam = cp_time(cp, tf)
                sta(t) = Λ_sta(t, cparam)
                println("Calculating robustness for time $tf")
                robustness_arr[index] = robustness(cparam, qts, sta)
        end
        return robustness_arr
end

function robustness_time(cp::ControlParameter, tarr::Vector{Float64}; nlambda::Int64=5, hessian::Bool=true, esta::Bool=true)
        if esta == true
                return robustness_time_esta(cp, tarr; nlambda=nlambda, hessian=hessian)
        else
                return robustness_time_sta(cp, tarr)
        end
end
