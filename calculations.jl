include("./src/ConstantQuantities.jl")
using DataFrames, CSV

filename(np::Int64) = "./Documents/Notes/Chapters/Results/gfx/data/fidelity_np$np.dat"
data(np::Int64) = CSV.read(filename(np), DataFrame; header=true)
time_filter(time_input::Vector{Float64}, target::Vector{Float64}) = time_input[time_input.∉[target]]
time_filter(data::DataFrame, time_input::Vector{Float64}) = time_filter(time_input, data.tf)

function fidelity_all(np::Int64, tspan::Vector{Float64}, Λf::Float64=50.0; nlambda::Int64=5, ω0::Float64=2.0)
    f = ControlParameterFull(ω0, 2.0√(Λf + 1.0), 0.2, np)
    i = ControlParameterInt(ω0, 2.0√(Λf + 1.0), 0.2, np)
    df = DataFrame(
        tf=tspan,
        full=fidelity_time(f, tspan; nlambda=nlambda), # Fidelity of eSTA with the full Hamiltonian and hessian
        full_orig=fidelity_time(f, tspan; nlambda=nlambda, hessian=false), # Fidelity of eSTA with the full Hamilhonian and NO hessian
        interm=fidelity_time(i, tspan; nlambda=nlambda), # Fidelity of eSTA with the intermediate Hamiltonian and hessian
        interm_orig=fidelity_time(i, tspan; nlambda=nlambda, hessian=false), # Fidelity of eSTA with the intermediate Hamilhonian and NO hessian
        sta=fidelity_time(f, tspan; esta=false) # STA fidelity ( as esta is set to false)
    )
    return df
end
function time_check(np::Int64, time_input::Vector{Float64})
    reduced_time = time_filter(data(np), time_input)
    isempty(reduced_time) ? println("nothing to do") : return fidelity_all(np, reduced_time)
end

fidelity_write(np::Int64, df::DataFrame) = isfile(filename(np)) ? CSV.write(filename(np), df, append=true) : CSV.write(filename(np), df)
fidelity_write(np::Int64, df::Nothing) = "Nothing to be done here, your file already has the fidelities for the requested final times"
fidelity(np::Int64, time_input::Vector{Float64}) = isfile(filename(np)) ? time_check(np, time_input) : fidelity_all(np, time_input)


fidelity_write(20, fidelity(20, range(0.05, 0.5, length=100) |> collect))

# I'm not going to split between eSTA and STA. Everything will be stored in the same DataFrame.
# Let's set up the names of the columns
# tf = final time -> the final time at which the feature is calculated
# N = Number of particles -> the number of particles used in the simulation
# λ = Number of correction points -> the number of correction points used in the eSTA method
# F_eSTA = fidelity of eSTA -> the fidelity of the eSTA method for the given final time
# F_STA = fidelity of STA -> the fidelity of the STA method for the given final time
# F_ad = fidelity of adiabatic -> the fidelity of the adiabatic method for the given final time
# Wn_eSTA = Sensitivity of eSTA -> the sensitivity of the eSTA method with respect to the white noise
# Wn_STA = Sensitivity of STA -> the sensitivity of the STA method with respect to the white noise
# Tn_eSTA = Sensitivity of eSTA -> the sensitivity of the eSTA method with respect to the time noise
# Tn_STA = Sensitivity of STA -> the sensitivity of the STA method with respect to the time noise
# Mn_eSTA = Sensitivity of eSTA -> the sensitivity of the eSTA method with respect to the modulation noise
# Mn_STA = Sensitivity of STA -> the sensitivity of the STA method with respect to the modulation noise
# Sq_eSTA = Squeezing of eSTA -> the Squeezing of the eSTA method 
# Sq_STA = Squeezing of STA -> the Squeezing of the STA method
# Λf = Final value of the control parameter -> the final value of the control parameter at the end of the simulation
# Λ0 = Initial value of the control parameter -> the initial value of the control parameter at the beginning of the simulation 
# The file will be called `whole_data.dat` and will be stored in the `data` folder
data_name = "./data/whole_data.dat"
# Let's create the empty DataFrame
empty() = DataFrame(tf = Float64[],N = Int64[], λ = Int64[], F_eSTA = Float64[], F_STA = Float64[],F_ad=Float64[], Wn_eSTA = Float64[], Wn_STA = Float64[], Tn_eSTA = Float64[], Tn_STA = Float64[], Mn_eSTA = Float64[], Mn_STA = Float64[], Sq_eSTA = Float64[], Sq_STA = Float64[], Λf = Float64[], Λ0 = Float64[])
Λfromω(ω) = 0.25 * ω - 1.0
ωfromΛ(Λ) = 2.0√(Λ + 1.0)
# Here I will define a function to calculate all the features, and will use the same name of the name of the respective column in the dataframe 
"""
calculate_all!(cp::ControlParameter, final_times; nlambda::Int64=5)
Calculates the fidelity, sensitivity and squeezing for the STA and eSTA methods for the given final times.
It returns a new DataFrame with the calculated features, so that it can be appended to an existing DataFrame or saved in a new file.
It takes as input:
	- cp: the control parameter used in the simulation
	- final_times: the final times at which the features are calculated. It can be a vector or a range.
It also supports the following keyword argument:
	- nlambda: the number of correction points used in the eSTA method
"""
function calculate_all(cp::ControlParameter, final_times::AbstractVector; nlambda::Int64=5)
    qts = ConstantQuantities(cp) # Constant quantities not depending on the final time
	df = empty()
	p=Progress(length(final_times), 1, "Calculating features for all final times")
	Threads.@threads for tf in final_times
        cparam = cp_time(cp,tf) # Control parameter with the new final time
        corrs = corrections(cparam; nlambda=nlambda) # Corrections for the eSTA protocol for the new final time
        esta(t) = Λ_esta(t, cparam, corrs) # eSTA control function for the new final time and corrections
        sta(t) = Λ_sta(t, cparam) # STA control function for the new final time 
		ad(t) = Λ_ad(t, cparam) # Adiabatic control function for the new final time
		row=Dict(
			:tf => tf,
			:N => cp.NParticles,
			:λ => nlambda,
			:F_eSTA => fidelity(cparam, qts,esta),
			:F_STA => fidelity(cparam, qts,sta),
			:F_ad => fidelity(cparam, qts,ad),
			:Wn_eSTA => sensitivity(cparam, esta),
			:Wn_STA => sensitivity(cparam, sta),
			:Tn_eSTA => robustness_tn(cparam, qts, esta, 1e-7),
			:Tn_STA => robustness_tn(cparam, qts, sta, 1e-7),
			:Mn_eSTA => robustness_mn(cparam, qts, esta),
			:Mn_STA => robustness_mn(cparam, qts, sta),
			:Sq_eSTA => squeezing(cparam, qts, esta,2)[2][2],
			:Sq_STA => squeezing(cparam, qts, sta,2)[2][2],
			:Λf => Λfromω(cparam.ωf),
			:Λ0 => Λfromω(cparam.ω0)
		)
	df=	vcat(df, DataFrame(row))
	next!(p)
	end
	return df
end

final_times = range(0.05, 0.5, length=400) |> collect
cp=ControlParameterFull(.1,10)
df = calculate_all(cp, final_times)

CSV.write(data_name, df)
CSV.write(data_name, df, append=true)

# Up until this point, I have calculated the features for all the final times and stored them in a DataFrame.
# Now I will focus on those quantitities that need to plotted as a function of time, not only depending on the final time.
# Those are the squeezing ξs and the control function Λ(t)
evo_name = "./data/evo_data.dat"
# Again, I will store all the relevant features and the constant quantities in a DataFrame
# The only big difference in this case is that I will call the time t instead of tf as it is not the final time anymore
empty_evo() = DataFrame(t = Float64[],N = Int64[], λ = Int64[], ξs_eSTA = Float64[], ξs_STA = Float64[], Λ_eSTA = Float64[], Λ_STA = Float64[], Λ_ad = Float64[], Λ0 = Float64[], Λf = Float64[])

"""
calculate_all(cp::ControlParameter, final_time::Float64; nlambda::Int64=5, steps::Int64=1000)
Calculates the squeezing and the control function for the STA and eSTA methods for the given final time, as well as the control function for the adiabatic approach.
This function has the same name as the other one, but it takes a single final time as input.
Again, it returns a new DataFrame with the calculated features, so that it can be appended to an existing DataFrame or saved in a new file.
It takes as input:
	- `cp`: the control parameter used in the simulation
	- `final_time`: the final time at which the features are calculated.
It also supports the following keyword argument:
	- nlambda: the number of correction points used in the eSTA method
	- steps: the number of steps used in the simulation
"""
function calculate_all(cp::ControlParameter, final_time::Float64; nlambda::Int64=5, steps::Int64=1000)
	qts = ConstantQuantities(cp) # Constant quantities not depending on the final time
	df = empty_evo()
	corrs = corrections(cp; nlambda=nlambda) # Corrections for the eSTA protocol for a given final time
	times = range(0.0, final_time, length=steps) # Times at which the features are calculated
	cparam = cp_time(cp,final_time) # Control parameter with the new final time
	corrs = corrections(cparam; nlambda=nlambda) # Corrections for the eSTA protocol for the new final time
		esta(t) = Λ_esta(t, cparam, corrs) # eSTA control function for the new final time and corrections
		sta(t) = Λ_sta(t, cparam) # STA control function for the new final time 
		ad(t) = Λ_ad(t, cparam) # Adiabatic control function for the new final time
	df=	vcat(df, DataFrame(row))
	next!(p)
	return df
end


