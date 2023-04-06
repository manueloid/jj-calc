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
        ax = @pgf Axis({title = "$(title) $np particles", legend_pos = "north east", xlabel = "t"}, plots, legends)
        filename_save = "$(feature_string)_np$(np)_nlambda$(nlambda)"
        pgfsave(filename_save * ".pdf", ax)
        println("Plot saved in $filename_save.pdf")
end
function compare_fidelity_plot(np::Int64, protocols::Vector=[:full, :full_orig, :interm, :interm_orig, :sta], nlambda::Int64=5)
        whole = PlotObjectsFid(np, nlambda)
        plots = []
        legends = []
        feature_string = "fidelity"
        title = uppercasefirst(feature_string)
        for symbols in protocols
                push!(plots, plot(whole[symbols]))
                push!(legends, legend(whole[symbols]))
        end
        ax = @pgf Axis({title = "$(title) $np particles", legend_pos = "south east", xlabel = "\$t/\\tau\$"}, plots, legends)
        filename_save = "$(feature_string)_np$(np)_nlambda$(nlambda)"
        pgfsave(filename_save * ".pdf", ax)
        println("Plot saved in $filename_save.pdf")
end
function compare_sensitivity_plot(np::Int64, protocols::Vector=[:full, :full_orig, :interm, :interm_orig, :sta], nlambda::Int64=5)
        whole = PlotObjectsSens(np, nlambda)
        plots = []
        legends = []
        feature_string = "sensitivity"
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
##
protocols = [:full, :sta]
compare_robustness_plot(10, protocols, 5)
compare_fidelity_plot(10, protocols, 5)

## Function definition to plot the comparison with different λ for:
# - Fidelity
# - Robustness
# - Sensitivity
function compare_fidelity_lambda(protocol::Symbol, np::Int64, nlambda::Vector{Int64})
        plots = []
        leg = []
        for λ in λs
                push!(plots, plot(PlotObjFid(protocol, np, λ))) # I will only focus on the :full protocol
                push!(leg, LegendEntry("\$\\lambda $λ\$"))
        end
        ax = @pgf Axis({legend_pos = "south east", title = "Fidelity $np particles"}, plots, leg)
        filename = "fidelity_compare$np.pdf"
        pgfsave("fidelity_compare$np.pdf", ax)
        println("Plot saved in $filename")
end
function compare_robustness_lambda(protocol::Symbol, np::Int64, nlambda::Vector{Int64})
        plots = []
        leg = []
        for λ in λs
                push!(plots, plot(PlotObjRob(protocol, np, λ))) # I will only focus on the :full protocol
                push!(leg, LegendEntry("\$\\lambda $λ\$"))
        end
        local ax = @pgf Axis({legend_pos = "north east", title = "Robustness $np particles"}, plots, leg)
        filename = "robustness_compare$np.pdf"
        pgfsave(filename, ax)
        println("Plot saved in $filename")
end
function compare_sensitivity_lambda(protocol::Symbol, np::Int64, nlambda::Vector{Int64})
        plots = []
        leg = []
        for λ in λs
                F = PlotObjFid(protocol, np, λ).feature # Fidelity for given protocol, number of particle and λ
                R = PlotObjRob(protocol, np, λ).feature  # Robustness for given protocol, number of particle and λ
                η = sqrt.((1 .- F) .^ 2 + R .^ 2) # Calculation of the sensitivity
                tf = PlotObjFid(protocol, np, λ).timespan # Fidelity for given protocol, number of particle and λ
                # Now I have the sensitivity and I can plot it.
                # I need to get the timespan vector and I will create a new PlotObjFid object (I could have used a PlotObjRob)
                object = PlotObjFid(Options(λ), η, tf)
                push!(plots, plot(object))
                push!(leg, LegendEntry("\$\\lambda $λ\$"))
        end
        local ax = @pgf Axis({legend_pos = "north east", title = "Sensitivity $np particles"}, plots, leg)
        filename = "sensitivity_compare$np.pdf"
        pgfsave(filename, ax)
        println("Plot saved in $filename")
end
## Actual plotting
λs = [1, 5, 8]
np = 10
compare_robustness_lambda(:full, 50, λs)
compare_fidelity_lambda(:full, 50, λs)
compare_sensitivity_lambda(:full, 30, [1, 5, 8])
##

compare_sensitivity_plot(50, protocols, 5)
