#= Here I need to evaluate the number squeezing parameter, i.e. ξ² = ΔJz² \ (N\4) 
What I want is a simulation of the Schrödinger equation that given the final time,
the boundary conditions and the number of particles, calculates the corrections and then saved the desired quantity.
I will only simulate the full Hessian and the STA Protocol to compare the two =#
using JosephsonSTA
include("src/ConstantQuantities.jl")
Λ_sta(t::Float64, cp::ControlParameter) = control_ω(t, cp) * 0.25 - 1.0
ω_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = control_ω(t, cp) - correction_poly(t, cp, corrections::Vector{Float64})
Λ_esta(t::Float64, cp::ControlParameter, corrections::Vector{Float64}) = ω_esta(t, cp, corrections) * 0.25 - 1.0
# squeezing return the squeezing parameter ξ with respect of time. ξ(t)
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
function squeezing(cp::ControlParameter)
    qts = ConstantQuantities(cp)
    ΔJ = dagger(qts.ψf) * (qts.Jz)^2 * qts.ψf - (dagger(qts.ψf) * (qts.Jz) * qts.ψf)^2 |> real
    h = 2.0 / cp.NParticles
    return ΔJ * 2.0h
end
#The two following functions are different in the sense that they define a different control paramater depending on the use of eSTA or not. 
function squeezing_esta(cp::ControlParameter, npoints=100)
    qts = ConstantQuantities(cp)
    corrs = corrections(cp)
    Λ(t) = Λ_esta(t, cp, corrs)
    return squeezing(cp, qts, Λ, npoints)
end
function squeezing_sta(cp::ControlParameter, npoints=100)
    qts = ConstantQuantities(cp)
    corrs = corrections(cp)
    Λ(t) = Λ_sta(t, cp)
    return squeezing(cp, qts, Λ, npoints)
end

# finding and plotting the evolution of squeezing depending on time
f = ControlParameterFull(0.2, 30)
t, ξ = squeezing_esta(f, 100)
t, ξs = squeezing_sta(f, 100)
using Plots
plot(t, ξ, label="eSTA", xlabel="t", ylabel="ξ", title="time evolution ξ")
plot!(t, ξs, label="STA")

# finding the squeezing at final time for different final times 
final_times = range(0.2, 0.8, length=100)
ξf = zeros(length(final_times))
ξfs = zeros(length(final_times))
ξfi = zeros(length(final_times))
for (index, final_time) in enumerate(final_times)
    local f = ControlParameterFull(final_time, 20)
    ξf[index] = squeezing_esta(f, 2)[2][2]
    ξfs[index] = squeezing_sta(f, 2)[2][2]
    ξfi[index] = squeezing(f)
end

#Plotting the results comparing the final squeezing for different final times
plot(final_times, ξf, label="eSTA", xlabel="t", ylabel="ξ", title="final ξ")
plot!(final_times, ξfs, label="STA")
plot!(final_times, ξfi, label="Ideal")
