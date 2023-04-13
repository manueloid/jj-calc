include("../src/ConstantQuantities.jl")
"""
robustness_timenoise(cp::ControlParameter, qts::ConstantQuantities, Λ::Function; order::Int64=5) -> Float64
Calculates the robustness of the control parameter cp to time noise, using the function Λ.
It does this by calculating the derivative of the fidelity with respect to the time noise δ, and then calculating the absolute value of that derivative.
It does this by using the central finite difference method, with order `order`.
"""
function robustness_timenoise(cp::ControlParameter, qts::ConstantQuantities, Λ::Function; order::Int64=5)
    function error(δ) # Evaluate the fidelity with shift δ in time
        Λ_err(t) = Λ(δ + t)
        return fidelity(cp, qts, Λ_err)
    end
    m = central_fdm(order, 1) # Create a central finite difference method
    return m(error, 0.0) |> abs # Calculate the derivative of the fidelity with respect to δ, and take the absolute value
end
"""
robustnesses_time(cp::ControlParameter; t::Vector{Float64}; nlambda::Int64=5) -> Vector{Float64}
Evaluates the robustness of the control parameter cp to time noise, using the eSTA corrections and the STA corrections.
"""
function robustnesses_time(cp::ControlParameter, final_times; nlambda::Int64=5)
	qts = ConstantQuantities(cp)
	robustness_sta = zeros(length(final_times))
	robustness_esta = zeros(length(final_times))
	Threads.@threads for (index, tf) in enumerate(final_times) |> collect
		cparam = cp_time(cp, tf)
		corrs = corrections(cparam; nlambda=nlambda)
		esta(t) = Λ_esta(t, cparam, corrs)
		println("Calculating robustness for time $tf, eSTA")
		sta(t) = Λ_sta(t, cparam)
		println("Calculating robustness for time $tf, STA")
		robustness_sta[index] = robustness_timenoise(cparam, qts, sta)
		robustness_esta[index] = robustness_timenoise(cparam, qts, esta)
	end
	return robustness_sta, robustness_esta
end
"""
robustnesses_time(n::Int64, final_times; nlambda::Int64=5) -> Vector{Float64}
Given a number of particles n, evaluates the robustness for both the STA and eSTA corrections, for a range of times.
"""
function robustnesses_time(n::Int64, final_times; nlambda::Int64=5)
	cp = ControlParameterFull(0.2, n)
	return robustnesses_time(cp, final_times; nlambda=nlambda)
end

tspan = range(0.2, 0.85, length=100) |> collect
sta, esta = robustnesses_time(cp, tspan)
using Plots
plot(tspan, esta, label="eSTA")
plot!(tspan, sta, label="STA")
