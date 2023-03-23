include("src/ConstantQuantities.jl")

function robustness_timenoise(cp::ControlParameter, qts::ConstantQuantities, Λ::Function; order::Int64=5)
        function error(δ)
                Λ_err(t) = Λ(δ + t)
                return fidelity(cp, qts, Λ_err)
        end
        m = central_fdm(order, 2)
        return m(error, 0.0) |> abs
end
function robustness_timenoise_esta(cp::ControlParameter, tarr::Vector{Float64}; nlambda::Int64=5, hessian::Bool=true)
        qts = ConstantQuantities(cp)
        robustness_arr = zeros(length(tarr))
        Threads.@threads for (index, tf) in enumerate(tarr) |> collect
                cparam = cp_time(cp, tf)
                corrs = corrections(cparam; hessian=hessian, nlambda=nlambda)
                esta(t) = Λ_esta(t, cparam, corrs)
                println("Calculating robustness for time $tf")
                robustness_arr[index] = robustness_timenoise(cparam, qts, esta)
        end
        return robustness_arr
end
function robustness_timenoise_sta(cp::ControlParameter, tarr::Vector{Float64})
        qts = ConstantQuantities(cp)
        robustness_arr = zeros(length(tarr))
        Threads.@threads for (index, tf) in enumerate(tarr) |> collect
                cparam = cp_time(cp, tf)
                sta(t) = Λ_sta(t, cparam)
                println("Calculating robustness for time $tf")
                robustness_arr[index] = robustness_timenoise(cparam, qts, sta)
        end
        return robustness_arr
end
function robustness_timenoise(cp::ControlParameter, tarr::Vector{Float64}; nlambda::Int64=5, hessian::Bool=true, esta::Bool=true)
        if esta == true
                return robustness_timenoise_esta(cp, tarr; nlambda=nlambda, hessian=hessian)
        else
                return robustness_timenoise_sta(cp, tarr)
        end
end

cp = ControlParameterFull(.2,40)
tspan = range(0.2,0.85, length = 100)|> collect
sta = robustness_timenoise(cp, tspan; esta = false)
esta = robustness_timenoise(cp, tspan)
using Plots
plot(tspan, esta, label = "eSTA")  
plot!(tspan, sta, label = "STA")  

cp = ControlParameterFull(.2,40)
qts = ConstantQuantities(cp)
corrs = corrections(cp) 
Λ(t)=Λ_esta(t+.1, cp, corrs) 
Λ(.2)

fidelity(cp, qts, Λ)
