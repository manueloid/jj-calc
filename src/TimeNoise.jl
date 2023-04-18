# Will use the subscript 'tn' to refer to the time noise
time_noise(δ::Float64, Λ::Function, t::Float64) = Λ(t + δ) # Time noise Λ needs to be dependent by only one parameter
"""
`fidelity_tn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, δ::Float64)`
Given the boundary conditions, and the control function Λ, return the fidelity of the process with a time noise of shift δ.
"""
function fidelity_tn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, δ::Float64)
    Λ_tn(t) = time_noise(δ, Λ, t) # Time noise function for the control function Λ and shift δ
    return fidelity(cp, qts, Λ_tn) # Fidelity of the process with the time noise
end

"""
robustness_tn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, δ::Float64)
Return the derivative of the fidelity given a time shift in the control parameter of the form Λδ(t) = Λ(t + δ)
"""
function robustness_tn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, δ::Float64)
    f(δ) = fidelity_tn(cp, qts, Λ, δ) # Fidelity function with the time noise of shift δ
    return finite_difference(f, δ) |> abs # Derivative of the fidelity with respect to the shift of the time noise
end

"""
robustnesses_tn(cp::ControlParameter, final_times, δ::Float64)
Return the derivative of the fidelity given a time shift in the control parameter of the form Λδ(t) = Λ(t + δ)
The derivative is calculated for a range of final times, for both the eSTA and STA protocols.
"""
function robustnesses_tn(cp::ControlParameter, final_times, δ::Float64)
    robustnesses_esta = zeros(length(final_times))
    robustnesses_sta = zeros(length(final_times))
    qts = ConstantQuantities(cp)
    Threads.@threads for i in 1:length(final_times)
        cparam = cp_time(cp, final_times[i])
        corrs = corrections(cparam)
        esta(t) = Λ_esta(t, cparam, corrs)
        sta(t) = Λ_sta(t, cparam)
        robustnesses_esta[i] = robustness_tn(cparam, qts, esta, δ)
        robustnesses_sta[i] = robustness_tn(cparam, qts, sta, δ)
    end
    return robustnesses_esta, robustnesses_sta
end
# up until now, the robustnesses are calculated for a single value of the modulation noise amplitude δ, where δ is fixed.
# From now on, the derivative will be calculated automatically

using FiniteDifferences
"""
`robustness_tn(cp::ControlParameter, qts::ConstantQuantities, ω::Function; order::Int64=5)`
return the derivative of the fidelity given an error in the modulation of the control parameter of the form Λδ(t) = (1 + δ)Λ(t)
# Arguments 
- `cp` 
- `qts` 
- `Λ` must be a function of only one parameter (ideally the time) and it has to be built outside the function declaration 
# Kwargs 
- `order` is the order of the finite difference method I need to use
"""
function robustness_tn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function; order::Int64=5)
    function error(δ)
        Λ_err(t) = Λ(t + δ)
        return fidelity(cp, qts, Λ_err)
    end
    m = central_fdm(order, 1)
    return m(error, 0.0) |> abs
end

"""
`robustnesses_tn(cp::ControlParameter, final_times; order::Int64=5)`
Return the derivative of the fidelity given an error in the modulation of the control parameter of the form Λδ(t) = (1 + δ)Λ(t).
The derivative is calculated for different final times and for both the eSTA and STA protocols using the finite difference method with an automatic grid definition
The output is a tuple of two arrays, the first one containing the robustnesses of the eSTA protocol and the second one containing the robustnesses of the STA protocol.
"""
function robustnesses_tn(cp::ControlParameter, final_times; order::Int64=5)
    qts = ConstantQuantities(cp) # Constant quantities for the control parameter cp
    robustnesses_esta = zeros(length(final_times)) # Array to store the robustnesses of the eSTA protocol
    robustnesses_sta = zeros(length(final_times)) # Array to store the robustnesses of the STA protocol
    Threads.@threads for i in 1:length(final_times) # Loop over the final times
        cparam = cp_time(cp, final_times[i]) # Control parameter with the new final time
        corrs = corrections(cparam) # Corrections for the control parameter cp with the new final time
        esta(t) = Λ_esta(t, cparam, corrs) # eSTA control function for the control parameter cp with the new final time
        sta(t) = Λ_sta(t, cparam) # STA control function for the control parameter cp with the new final time
        robustnesses_esta[i] = robustness_tn(cparam, qts, esta; order=order) # Robustness of the eSTA protocol
        println("Robustness of eSTA: ", robustnesses_esta[i], " for final time ", final_times[i])
        robustnesses_sta[i] = robustness_tn(cparam, qts, sta; order=order) # Robustness of the STA protocol
        println("Robustness of STA: ", robustnesses_sta[i], " for final time ", final_times[i])
    end
    return robustnesses_esta, robustnesses_sta # Return the robustnesses of the eSTA and STA protocols
end
