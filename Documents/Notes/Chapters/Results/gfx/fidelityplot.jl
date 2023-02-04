include("./Documents/Notes/Chapters/Results/gfx/PlotObj.jl")
options(opt_obj::Options) = @pgf {color = opt_obj.color, style = opt_obj.style, thick}
options(plt_obj::PlotObj) = options(plt_obj.options)

plot_fidelity(plt_obj::PlotObj) = @pgf Plot(options(plt_obj), Table([plt_obj.timespan, plt_obj.fidelity]))
plot_robustness(plt_obj::PlotObj) = @pgf Plot(options(plt_obj), Table([plt_obj.timespan, plt_obj.robustness]))
legend(plt_obj::PlotObj) = @pgf LegendEntry(plt_obj.options.name)

function export_fidelity(np::Int64, protocols::Vector=[:full, :full_orig, :interm, :interm_orig, :sta])
    whole = PlotObj(np)
    plots = []
    legends = []
    for symbols in protocols
        push!(plots, plot_fidelity(whole[symbols]))
        push!(legends, legend(whole[symbols]))
    end
    ax = @pgf Axis({title = "Fidelity $np particles", legend_pos = "south east"}, plots, legends)
    filename_save = "./Documents/Notes/Chapters/Results/gfx/fidelity$np.pdf"
    pgfsave(filename_save, ax)
    println("Plot saved in $filename_save")
end

function export_robustness(np::Int64, protocols::Vector=[:full, :full_orig, :interm, :interm_orig, :sta])
    whole = PlotObj(np)
    plots = []
    legends = []
    for symbols in protocols
        push!(plots, plot_robustness(whole[symbols]))
        push!(legends, legend(whole[symbols]))
    end
    ax = @pgf Axis({title = "Fidelity $np particles"}, plots, legends)
    filename_save = "./Documents/Notes/Chapters/Results/gfx/robustness$np.pdf"
    pgfsave(filename_save, ax)
    println("Plot saved in $filename_save")
end

export_fidelity(20)
