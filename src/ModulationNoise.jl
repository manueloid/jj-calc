# Will use the subscript 'mn' to refer to the modulation noise
modulation_noise(δ::Float64, Λ::Function, t::Float64) = (1 + δ) * Λ(t) # Modulation noise Λ needs to be dependent by only one parameter
"""
fidelity_mn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, δ::Float64)
Given the boundary conditions, and the control function Λ, return the fidelity of the process with a modulation noise of amplitude δ.
"""
function fidelity_mn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, δ::Float64)
    Λ_mn(t) = modulation_noise(δ, Λ, t) # Modulation noise function for the control function Λ and amplitude δ
    return fidelity(cp, qts, Λ_mn) # Fidelity of the process with the modulation noise
end

"""
finite_difference(f::Function, δ::Float64; x::Float64=0.0)
Return the derivative of the function f at the point x using the finite difference method
"""
function finite_difference(f::Function, δ::Float64; x::Float64=0.0)
    return (f(x + δ) - f(x - δ)) / (2 * δ)
end

"""
robustness_mn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, δ::Float64)
Return the derivative of the fidelity given an error in the modulation of the control parameter of the form Λδ(t) = (1 + δ)Λ(t)
"""
function robustness_mn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function, δ::Float64)
    f(δ) = fidelity_mn(cp, qts, Λ, δ) # Fidelity function with the modulation noise of amplitude δ
    return finite_difference(f, δ) |> abs # Derivative of the fidelity with respect to the amplitude of the modulation noise
end

"""
robustnesses_mn(cp::ControlParameter, final_times, δ::Float64)
Return the derivative of the fidelity given an error in the modulation of the control parameter of the form Λδ(t) = (1 + δ)Λ(t).
The derivative is calculated for different final times and for both the eSTA and STA protocols.
The output is a tuple of two arrays, the first one containing the robustnesses of the eSTA protocol and the second one containing the robustnesses of the STA protocol.
"""
function robustnesses_mn(cp::ControlParameter, final_times, δ::Float64)
    qts = ConstantQuantities(cp) # Constant quantities for the control parameter cp
    robustnesses_esta = zeros(length(final_times)) # Array to store the robustnesses of the eSTA protocol
    robustnesses_sta = zeros(length(final_times)) # Array to store the robustnesses of the STA protocol
    p = Progress(length(final_times), 1, "Calculating robustnesses")
    Threads.@threads for i in 1:length(final_times) # Loop over the final times
        cparam = cp_time(cp, final_times[i]) # Control parameter with the new final time
        corrs = corrections(cparam) # Corrections for the control parameter cp with the new final time
        esta(t) = Λ_esta(t, cparam, corrs) # eSTA control function for the control parameter cp with the new final time
        sta(t) = Λ_sta(t, cparam) # STA control function for the control parameter cp with the new final time
        robustnesses_esta[i] = robustness_mn(cparam, qts, esta, δ) # Robustness of the eSTA protocol
        robustnesses_sta[i] = robustness_mn(cparam, qts, sta, δ) # Robustness of the STA protocol
        next!(p)
    end
    return robustnesses_esta, robustnesses_sta # Return the robustnesses of the eSTA and STA protocols
end
# up until now, the robustnesses are calculated for a single value of the modulation noise amplitude δ, where δ is fixed.
# From now on, the derivative will be calculated automatically

using FiniteDifferences
"""
`robustness_mn(cp::ControlParameter, qts::ConstantQuantities, ω::Function; order::Int64=5)`
return the derivative of the fidelity given an error in the modulation of the control parameter of the form Λδ(t) = (1 + δ)Λ(t)
# Arguments 
- `cp` 
- `qts` 
- `Λ` must be a function of only one parameter (ideally the time) and it has to be built outside the function declaration 
# Kwargs 
- `order` is the order of the finite difference method I need to use
"""
function robustness_mn(cp::ControlParameter, qts::ConstantQuantities, Λ::Function; order::Int64=5)
    function error(δ)
        Λ_err(t) = (1 + δ) * Λ(t)
        return fidelity(cp, qts, Λ_err)
    end
    m = central_fdm(order, 1)
    return m(error, 0.0) |> abs
end

"""
`robustnesses_mn(cp::ControlParameter, final_times; order::Int64=5)`
Return the derivative of the fidelity given an error in the modulation of the control parameter of the form Λδ(t) = (1 + δ)Λ(t).
The derivative is calculated for different final times and for both the eSTA and STA protocols using the finite difference method with an automatic grid definition
The output is a tuple of two arrays, the first one containing the robustnesses of the eSTA protocol and the second one containing the robustnesses of the STA protocol.
"""
function robustnesses_mn(cp::ControlParameter, final_times; order::Int64=5)
    qts = ConstantQuantities(cp) # Constant quantities for the control parameter cp
    robustnesses_esta = zeros(length(final_times)) # Array to store the robustnesses of the eSTA protocol
    robustnesses_sta = zeros(length(final_times)) # Array to store the robustnesses of the STA protocol
    p = Progress(length(final_times), 1, "Calculating robustnesses")
    Threads.@threads for i in 1:length(final_times) # Loop over the final times
        cparam = cp_time(cp, final_times[i]) # Control parameter with the new final time
        corrs = corrections(cparam) # Corrections for the control parameter cp with the new final time
        esta(t) = Λ_esta(t, cparam, corrs) # eSTA control function for the control parameter cp with the new final time
        sta(t) = Λ_sta(t, cparam) # STA control function for the control parameter cp with the new final time
        robustnesses_esta[i] = robustness_mn(cparam, qts, esta; order=order) # Robustness of the eSTA protocol
        robustnesses_sta[i] = robustness_mn(cparam, qts, sta; order=order) # Robustness of the STA protocol
        next!(p)
    end
    return robustnesses_esta, robustnesses_sta # Return the robustnesses of the eSTA and STA protocols
end
