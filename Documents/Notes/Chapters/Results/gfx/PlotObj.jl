using DataFrames, CSV
data_dir(nlambda::Int64) = "./data/nlambda$(nlambda)/" # Saving directory for different lambda
filename(np::Int64, feature="fidelity", nlambda::Int64=5) = data_dir(nlambda) * "$(feature)_np$np.dat" # Saving filename for given number of lambda and given number of particles
data(np::Int64, feature="fidelity", nlambda::Int64=5) = CSV.read(filename(np, feature, nlambda), DataFrame; header=true) |> s -> sort(s, :tf)
using Colors, PGFPlotsX
# Colors for the plots
colors = (
    colorant"#FF0000", # eSTA Full Hamiltonian with Hessian 
    colorant"#FFa000", # eSTA Full Hamiltonian with original version
    colorant"#0000FF", # eSTA Intermediate Hamiltonian with Hessian
    colorant"#00FF00",# eSTA Intermediate Hamiltonian original version
    colorant"#000000" # STA
)

# Array where all the line styles are stored
styles = (
    solid="solid",                         # eSTA Intermediate Hamiltonian with Hessian
    ldash="dash pattern={on 2pt off 2pt}",# eSTA Full Hamiltonian with Hessian 
    dash=" dashed",                        # eSTA Full Hamiltonian with original version
    dot="dash pattern={on 1pt off 1pt}",  # eSTA Intermediate Hamiltonian original version
    dot_dash="dash pattern={on 4pt off 1pt on 1pt off 1pt}" #STA
)

names = (
    "Full Hessian",
    "Full Original",
    "Intermediate Hessian",
    "Intermediate Original",
    "STA"
)

struct Options # Structure to hold the values of the Option
    name::String
    color::RGB
    style::String
end
pick(number::Int64, array) = array[number%length(array)+1]
# Options utility to create the options for different lambdas given the arrays `colors` and `styles`
Options(nlambda::Int64) = Options("$nlambda lambdas", pick(nlambda, colors), pick(nlambda, styles))

abstract type PlotObj end

struct PlotObjFid <: PlotObj
    options::Options
    feature::Vector{Float64}
    timespan::Vector{Float64}
end

struct PlotObjRob <: PlotObj
    options::Options
    feature::Vector{Float64}
    timespan::Vector{Float64}
end

PlotObjFid(options::Options, protocol::Symbol, input_data::DataFrame) = PlotObjFid(options, input_data[!, protocol], input_data[!, :tf])
PlotObjRob(options::Options, protocol::Symbol, input_data::DataFrame) = PlotObjRob(options, input_data[!, protocol], input_data[!, :tf])


function PlotObjFid(options::Options, protocol::Symbol, np::Int64, nlambda::Int64=5)
    input_data = data(np, "fidelity", nlambda)
    return PlotObjFid(options, protocol, input_data)
end

function PlotObjRob(options::Options, protocol::Symbol, np::Int64, nlambda::Int64=5)
    input_data = data(np, "robustness", nlambda)
    return PlotObjRob(options, protocol, input_data)
end

PlotObjFid(protocol::Symbol, np::Int64, nlambda::Int64) = PlotObjFid(Options(nlambda), protocol, np, nlambda)
PlotObjRob(protocol::Symbol, np::Int64, nlambda::Int64) = PlotObjRob(Options(nlambda), protocol, np, nlambda)

function PlotObjectsFid(np::Int64, nlambda::Int64=5)
    fidelities = data(np, "fidelity", nlambda)
    protocols = propertynames(fidelities)[2:end]
    objects = []
    for (index, protocol) in enumerate(protocols)
        options = Options(names[index], colors[index], styles[index])
        object = PlotObjFid(options, protocol, fidelities)
        push!(objects, object)
    end
    return (; zip(protocols, objects)...)
end

function PlotObjectsRob(np::Int64, nlambda::Int64=5)
    fidelities = data(np, "robustness", nlambda)
    protocols = propertynames(fidelities)[2:end]
    objects = []
    for (index, protocol) in enumerate(protocols)
        options = Options(names[index], colors[index], styles[index])
        object = PlotObjRob(options, protocol, fidelities)
        push!(objects, object)
    end
    return (; zip(protocols, objects)...)
end
