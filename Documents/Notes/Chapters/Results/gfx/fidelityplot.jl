include("PlotObj.jl")
options(opt_obj::Options) = @pgf {color = opt_obj.color, style = opt_obj.style, thick}
options(plt_obj::PlotObj) = options(plt_obj.options)

plot_feature(plt_obj::PlotObj, feature::Symbol) = @pgf Plot(options(plt_obj), Table([plt_obj.timespan, getfield(plt_obj, feature)]))
legend(plt_obj::PlotObj) = @pgf LegendEntry(plt_obj.options.name)

function compare_feature_plot(np::Int64, protocols::Vector=[:full, :full_orig, :interm, :interm_orig, :sta], feature::Symbol=:fidelity, nlambda::Int64=5)
    whole = PlotObj(np, nlambda)
    plots = []
    legends = []
    feature_string = String(feature)
    for symbols in protocols
        push!(plots, plot_feature(whole[symbols], feature))
        push!(legends, legend(whole[symbols]))
    end
    ax = @pgf Axis({title = "$(feature_string) $np particles", legend_pos = "south east"}, plots, legends)
    filename_save = "./Documents/Notes/Chapters/Results/gfx/$(feature_string)_np$(np)_nlambda$(nlambda)"
    pgfsave(filename_save * ".pdf", ax)
    println("Plot saved in $filename_save.pdf")
end

protocols = [:full, :full_orig, :interm, :interm_orig, :sta]
# compare_feature_plot(30, [:full], :robustness)

# I want a function that given a λ returns a new plotting object. I can do that bc I have defined multiple constructors

@pgf Axis(
    [plot_feature(PlotObj(:interm, 50, 1), :fidelity),
    plot_feature(PlotObj(:interm, 50, 8), :fidelity)]
)
