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
