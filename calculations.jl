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
# df = DataFrame(tf = Float64[], N = Int64[], λ = Int64[], F_eSTA = Float64[], F_STA = Float64[], Wn_eSTA = Float64[], Wn_STA = Float64[], Tn_eSTA = Float64[], Tn_STA = Float64[], Mn_eSTA = Float64[], Mn_STA = Float64[], Sq_eSTA = Float64[], Sq_STA = Float64[], Λf = Float64[], Λ0 = Float64[])

final_time = range(0.05, 0.5, length=100)
np = 10
nlambda = 5
Λfromω(ω) = 0.25 * ω - 1.0
ωfromΛ(Λ) = 2.0√(Λ + 1.0)
Λ0 = 0.0
Λf = 50.0
ω0 = ωfromΛ(Λ0)
ωf = ωfromΛ(Λf)
cparam = ControlParameterFull(ω0, ωf, final_time[1], np)
