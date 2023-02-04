using DataFrames, CSV
filename(np::Int64, feature="fidelity") = "./Documents/Notes/Chapters/Results/gfx/data/$(feature)_np$np.dat"
data(np::Int64, feature="fidelity") = CSV.read(filename(np, feature), DataFrame; header=true)

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

struct Options # temporary structure to hold the values of the Option
    name::String
    color::RGB
    style::String
end

struct PlotObj
    options::Options
    fidelity::Vector{Float64}
    robustness::Vector{Float64}
    timespan::Vector{Float64}
end

# Assuming the two datasets have already been loaded?
# I only want to pass the color and the style?
# I think so. I will load the data outside the declaration of the PlotObj object and only then I will create it
function PlotObj(options::Options, protocol::Symbol, fidelity::DataFrame, robustness::DataFrame)
    return PlotObj(
        options,
        fidelity[!, protocol],
        robustness[!, protocol],
        robustness[!, :tf]
    )
end

function PlotObj(options::Options, protocol::Symbol, np::Int64)
    return PlotObj(
        options,
        data(np)[!, protocol],
        data(np, "robustness")[!, protocol],
        data(np)[!, :tf],
    )
end

function PlotObj(np::Int64)
    fidelities = data(np, "fidelity")
    robustnesses = data(np, "robustness")
    protocols = propertynames(fidelities)[2:end]
    objects = []
    for (index, protocol) in enumerate(protocols)
        options = Options(names[index], colors[index], styles[index])
        object = PlotObj(options, protocol, fidelities, robustnesses)
        push!(objects, object)
    end
    return (; zip(protocols, objects)...)
end

