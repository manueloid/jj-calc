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
    whole = PlotObjectsFid(np, nlambda)
    plots = []
    legends = []
    feature_string = "fidelity"
    title = uppercasefirst(feature_string)
    for symbols in protocols
        push!(plots, plot(whole[symbols]))
        push!(legends, legend(whole[symbols]))
    end
    ax = @pgf Axis({title = "$(title) $np particles", legend_pos = "south east"}, plots, legends)
    filename_save = "$(feature_string)_np$(np)_nlambda$(nlambda)"
    pgfsave(filename_save * ".pdf", ax)
    println("Plot saved in $filename_save.pdf")
end

##
protocols = [:full, :full_orig, :interm, :interm_orig, :sta]
compare_robustness_plot(30, protocols, 5)
compare_fidelity_plot(30, protocols, 5)
