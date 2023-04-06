include("./src/ConstantQuantities.jl")
using DataFrames, CSV
data_dir(nlambda::Int64) = "./Documents/Notes/Chapters/Results/gfx/data/nlambda$(nlambda)/" # Saving directory for different lambda
filename(np::Int64, nlambda::Int64=5) = data_dir(nlambda) * "fidelity_np$np.dat" # Saving filename for given number of lambda and given number of particles

data(np::Int64, nlambda::Int64=5) = CSV.read(filename(np, nlambda), DataFrame; header=true) # Reading data from directory and turning it into a DataFrame
time_filter(time_input::Vector{Float64}, target::Vector{Float64}) = time_input[time_input.∉[target]] # Given two vectors an `input` and a `target`, returns a vector the elements of which are the elements of `input` that are not into `target`
time_filter(time_input::Vector{Float64}, data::DataFrame) = time_filter(time_input, data.tf) # Same as the line above but in this case the `target` is the column of final times of the loaded DataFrame

"""
function fidelity_all(np::Int64, tspan::Vector{Float64}, Λf::Float64=50.0; kwargs) For a given number of particles `np` and an initial value `Λf` give the fidelity of all the types of protocol for each time in `tspan`.
        It returns all the fidelities and the corresponding final times by storing them in a DataFrame object that is easy to manipulate
        kwargs:
        - `nlambda::Int64` number of points for which the corrections are calculated, default is 5
        - `ω0::Float64` initial value for ω0, which is usually set to 2.0
"""
function fidelity_all(np::Int64, tspan::Vector{Float64}, Λf::Float64=50.0; nlambda::Int64=5, ω0::Float64=2.0)
    f = ControlParameterFull(ω0, 2.0√(Λf + 1.0), 0.2, np)
    i = ControlParameterInt(ω0, 2.0√(Λf + 1.0), 0.2, np)
    df = DataFrame(
        tf=tspan,
        full=fidelity_time(f, tspan; nlambda=nlambda), # fidelity of eSTA with the full Hamiltonian and hessian
        full_orig=fidelity_time(f, tspan; nlambda=nlambda, hessian=false), # fidelity of eSTA with the full Hamilhonian and NO hessian
        interm=fidelity_time(i, tspan; nlambda=nlambda), # fidelity of eSTA with the intermediate Hamiltonian and hessian
        interm_orig=fidelity_time(i, tspan; nlambda=nlambda, hessian=false), # fidelity of eSTA with the intermediate Hamilhonian and NO hessian
        sta=fidelity_time(f, tspan; esta=false) # STA fidelity ( as esta is set to false)
    )
    return df
end

"""
function time_check(np::Int64, time_input::Vector{Float64}) given a number of particles `np` and an input vector of times `time_input`, checks in the corresponding data file if there is anything that needs to be calculated.
    It checks how many elements of `time_input` are in the `tf` column of the corresponding DataFrame. It then evaluates the fidelities for the resulting set and does nothing if the set is empty.
"""
function time_check(np::Int64, time_input::Vector{Float64}, nlambda::Int64=5)
    reduced_time = time_filter(time_input, data(np, nlambda))
    isempty(reduced_time) ? println("nothing to do") : return fidelity_all(np, reduced_time; nlambda=nlambda)
end

"""
fidelity(np::Int64, time_input::Vector{Float64}, nlambda::Int64 =5) It first checks if there is any data file for the given number of particles and the given number of lambdas. 
    If that is the case, it will check how many times we need to calculate the fidelity by filtering the values in `time_input` and the `tf` column in the corresponding DataFrame
    It returns an a DataFrame object that can then be saved in an other file
"""
fidelity(np::Int64, time_input::Vector{Float64}, nlambda::Int64=5) = isfile(filename(np, nlambda)) ? time_check(np, time_input, nlambda) : fidelity_all(np, time_input; nlambda=nlambda)

function fidelity_write(np::Int64, df::DataFrame, nlambda::Int64=5)
    if isfile(filename(np, nlambda))
        CSV.write(filename(np, nlambda), df, append=true)  # Appends a DataFrame object to a file if the file already exists, else it will create it and start writing to it
    else
        mkpath(data_dir(nlambda))
        CSV.write(filename(np, nlambda), df)
    end
end

fidelity_write(np::Int64, df::Nothing, nlambda::Int64=5) = "Nothing to be done here, your file already has the fidelities for the requested final times" # If the DataFrame object is empty, it means all the values of the `time_input` are alredy in the `tf` column of the DataFrame object, so there is nothing to do

# Actual Calculation
λs = 8
npmax = 10

[fidelity_write(np, fidelity(np, range(0.05, 0.5, length=100) |> collect, λs), λs) for np in 10:10:npmax]
