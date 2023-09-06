# I need to backpropagate the wavefunction from the final state to the initial one
# First, I will define a function that will perform the time reversal of a generic function
"""
time_reverse(t::Float64, f::Function, cp::ControlParameter) -> f(cp.final_time - t) 
This function performs the time reversal of a generic function f(t,cp) with respect to the final time of the whole protocol.
The idea is to reverse the functions that are used to define the control parameters of the protocol, i.e. ω or Λ
"""
function time_reverse(t::Float64, f::Function, cp::ControlParameter)
    tf = cp.final_time
    return f(tf - t)
end
"""
Hamiltonian(Jx, Jz, ω::Function, h::Float64) -> H(t, psi)
This function returns the time dependent Hamiltonian of the system.
Here is the breakdown of the arguments:
	- Jx, Jz: Angular momentum operators of the system
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
"""
H1(t::Float64, qts::ConstantQuantities, ω::Function) -> H1(t)
Return the Jz part of the Hamiltonian at time t.
It is important to point out that the scaling factor is not included here.
We need then to multiply the result of this function by h.
"""
function H1(t::Float64, qts::ConstantQuantities, ω)
    # return qts.Jz^2 * ω(t)
    return qts.Jz
    # try to do the same thing but only with Jz
    # return qts.Jz
end
tspan(cp::ControlParameter, steps::Int64) = range(0.0, cp.final_time, length=steps) |> collect
"""
bf_evolution(cp::ControlParameter, ω::Function; t=tspan(cp, 1000)) -> ψ0(t), ψT(t), H1(t)
Short for "Backward Forward Evolution".
This function performs the time evolution of the system, both forward and backward in time.
It returns the wavefunctions and the time dependent Hamiltonian.
They are returned in the following order:
	- ψ0: array containing the wavefunction of the forward evolution at each time step
	- ψT: array containing the wavefunction of the backward evolution at each time step
	- H1: array containing the time dependent Hamiltonian at each time step. This is the Jz^2 part of the Hamiltonian.
"""
function bf_evolution(cp::ControlParameter, ω::Function; t=tspan(cp, 1000))
    qts = ConstantQuantities(cp) # Constant quantities of the system: operators, initial and final states
    ω_rev(t) = time_reverse(t, ω, cp) # Time reversed control function
    h = 2.0 / cp.NParticles # Scaling factor
    H_fwd = Hamiltonian(qts, ω, h) # Forward Hamiltonian
    H_rev = Hamiltonian(qts, ω_rev, h) # Time reversed Hamiltonian
    ψ0 = qts.ψ0 # Initial state for the forward evolution
    ψT = qts.ψf # Initial state for the backward evolution
    tout, ψ0 = evolution(ψ0, H_fwd, t) # Time span and wavefunction with respect of time for the forward evolution
    tout, ψT = evolution(ψT, H_rev, t) # Time span and wavefunction with respect of time for the backward evolution
	reverse!(ψT) # I reverse the wavefunction of the backward evolution
    H1_fwd(t) = h * H1(t, qts, ω) # Forward Hamiltonian with respect to time, that I will use to compute the sensitivity
    return ψ0, ψT, H1_fwd.(t) # I return the wavefunctions and the time dependent Hamiltonian
end
# trapezoidal rule for integration
"""
trapezoidal_rule(f::Array{Float64,1}, t::Array{Float64,1}) -> I
This function performs the integration of a function f(t) using the trapezoidal rule assuming that both f and t are arrays
"""
function trapezoidal_rule(f::Array{Float64,1}, t::Array{Float64,1})
    n = length(f)
    h = t[2] - t[1]
    return h * (sum(f[2:n-1]) + 0.5 * f[n] + 0.5 * f[1])
end
# simpson rule for integration
"""
simpson_rule(f::Array{Float64,1}, t::Array{Float64,1}) -> I
This function performs the integration of a function f(t) using the Simpson rule assuming that both f and t are arrays
"""
function simpson_rule(f::Array{Float64,1}, t::Array{Float64,1})
    n = length(f)
    h = t[2] - t[1]
    return h / 3.0 * (f[1] + 4.0 * sum(f[2:2:n-1]) + 2.0 * sum(f[3:2:n-2]) + f[n])
end
function sensitivity(cp::ControlParameter, ω::Function; t=tspan(cp, 1000))
    ψ0, ψf, H1_in = bf_evolution(cp, ω) # Results of the backward forward evolution
    ψ0d = dagger.(ψ0) # Conjugate transpose of the wavefunction of the forward evolution
    ψfd = dagger.(ψf) # Conjugate transpose of the wavefunction of the backward evolution
    ψ0dψf = ψ0d .* ψf # <ψ0|ψf> term
    ψfH12ψ0 = ψfd .* H1_in .^ 2 .* ψ0 # <ψf|H1^2|ψ0> term
    first = real.(ψ0dψf .* ψfH12ψ0) # real part of the <ψ0|ψf><ψf|H1^2|ψ0> term
    ψ0dH1ψf = ψ0d .* H1_in .* ψf # <ψ0|H1|ψf> term that will be subtracted from the previous term
    second = abs2.(ψ0dH1ψf) # |<ψ0|H1|ψf>|^2 term
    h = 2.0 / cp.NParticles # Scaling factor
    integrand =first .- second ./ h^2 # I compute the integrand
	return abs(simpson_rule(integrand, t)) # I integrate the integrand using the Simpson rule
end
# I need to integrate the integrand array over the time.
# I will use the trapezoidal rule
"""
sensitivities(n::Int64, final_times) -> sensitivities
This function computes the sensitivities for a given number of particles n and a given array of final times.
It returns an array containing the sensitivities for each final time for both the esta and sta methods.
"""
function sensitivities(n::Int64, final_times)
    sensitivities_esta = zeros(length(final_times))
    sensitivities_sta = zeros(length(final_times))
	p = Progress(length(final_times), 1, "Computing sensitivities")
    Threads.@threads for i in 1:length(final_times)
        cp = ControlParameterFull(final_times[i], n)
        corrs = corrections(cp)
        esta(t) = Λ_esta(t, cp, corrs)
        sta(t) = Λ_sta(t, cp)
        sensitivities_esta[i] = sensitivity(cp, esta; t=tspan(cp, 10000))
        sensitivities_sta[i] = sensitivity(cp, sta, t=tspan(cp, 10000))
        next!(p)
    end
    return sensitivities_esta, sensitivities_sta
end
