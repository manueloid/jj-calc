using Plots
include("./src/ConstantQuantities.jl")
# What I want is a fidelity that takes a control parameter and the control functions and spits out the  . 
function rollout(qts::ConstantQuantities)
    return [getfield(qts, field) for field in fieldnames(typeof(qts))]
end

function fidelity(cp::ControlParameter, qts::ConstantQuantities, ω::Function; nlambda::Int64=5, maxbra::Int64=4, hessian::Bool=true)
    Jz, Jx, ψ0, ψf = rollout(qts)
    tf = cp.final_time
    N = cp.NParticles
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
        return (2.0 / N * (ω(t) / 4.0 - 1.0) * Jz^2 - 2.0 * Jx)
    end
    fidelity(t, psi) = abs2.(dagger(ψf) * psi)
    return timeevolution.schroedinger_dynamic([0.0, tf], ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
end

function fidelity_noise(cp::ControlParameter, ϵ)
    qts = ConstantQuantities(cp)
    corrs = corrections_hess(cp)
    ω(t) = control_ω(t, cp) * (1 + ϵ) - correction_poly(t, cp, corrs)
    return fidelity(cp, qts, ω)
end


function fidelity_noise(cp::ControlParameter, epsilon::Vector{Float64})
    qts = ConstantQuantities(cp)
    corrs = corrections_hess(cp)
    fidelity_array = zeros(length(epsilon))
    Threads.@threads for (index, ϵ) in enumerate(epsilon) |> collect
        ω(t) = control_ω(t, cp) * (1 + ϵ) - correction_poly(t, cp, corrs)
        fidelity_array[index] = fidelity(cp, qts, ω)
    end
    return fidelity_array
end

f = ControlParameterFull(0.3, 20)
epsilon = range(-0.1, 0.1, length=99) |> collect
fid = fidelity_noise(f, epsilon)
plot(epsilon, fid)
