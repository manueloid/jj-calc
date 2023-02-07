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
include("PlotObj.jl")
options(opt_obj::Options) = @pgf {color = opt_obj.color, style = opt_obj.style, thick}
options(plt_obj::PlotObj) = options(plt_obj.options)

plot(plt_obj::PlotObj) = @pgf Plot(options(plt_obj), Table([plt_obj.timespan, plt_obj.feature]))
legend(plt_obj::PlotObj) = @pgf LegendEntry(plt_obj.options.name)

function compare_robustness_plot(np::Int64, protocols::Vector=[:full, :full_orig, :interm, :interm_orig, :sta], nlambda::Int64=5)
    whole = PlotObjectsRob(np, nlambda)
    plots = []
    legends = []
    feature_string = "robustness"
    title = uppercasefirst(feature_string)
    for symbols in protocols
        push!(plots, plot(whole[symbols]))
        push!(legends, legend(whole[symbols]))
    end
    ax = @pgf Axis({title = "$(title) $np particles", legend_pos = "north east"}, plots, legends)
    filename_save = "$(feature_string)_np$(np)_nlambda$(nlambda)"
    pgfsave(filename_save * ".pdf", ax)
    println("Plot saved in $filename_save.pdf")
end

function compare_fidelity_plot(np::Int64, protocols::Vector=[:full, :full_orig, :interm, :interm_orig, :sta], nlambda::Int64=5)
    whole = PlotObjectsRob(np, nlambda)
    plots = []
    legends = []
    feature_string = "fidelity"
    title = uppercasefirst(feature_string)
    for symbols in protocols
        push!(plots, plot(whole[symbols]))
        push!(legends, legend(whole[symbols]))
    end
    ax = @pgf Axis({title = "$(title) $np particles, $(nlambda) \\lambda", legend_pos = "south east"}, plots, legends)
    filename_save = "$(feature_string)_np$(np)_nlambda$(nlambda)"
    pgfsave(filename_save * ".pdf", ax)
    println("Plot saved in $filename_save.pdf")
end
protocols = [:full, :full_orig, :interm, :interm_orig, :sta]
compare_fidelity_plot(30, protocols, 1)


function PlotObjectsRob(np::Int64, nlambda::Int64=5)
    robustnesses = data(np, "robustness", nlambda)
    protocols = propertynames(robustnesses)[2:end]
    objects = []
    for (index, protocol) in enumerate(protocols)
        options = Options(names[index], colors[index], styles[index])
        object = PlotObjFid(options, protocol, robustnesses)
        push!(objects, object)
    end
    return (; zip(protocols, objects)...)
end
