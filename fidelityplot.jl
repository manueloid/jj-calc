include("./src/ConstantQuantities.jl")
include("./src/PlotObj.jl")
"""
plot_fidelity(n::Int64, final_times) -> Axis
Calculate the fidelity of the eSTA and STA for a given number of particles.
It returns an Axis object that can be saved with the pgfsave function.
"""
function plot_fidelity(n::Int64, final_times)
    # Calculate the fidelity
    cp = ControlParameterFull(0.1, n)
    esta_fid, sta_fid = fidelities(cp, final_times)
    esta_fid1 = fidelities(cp, final_times; nlambda=1)[1]
    # Plotting with PGFPlotsX
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles[1]}
    esta_opt1 = @pgf {color = colors.blue, line_width = 1.5, style = styles.[3]}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles[2]}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "F", legend_pos = "south east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_fid)),
        Plot(sta_opt, Coordinates(final_times, sta_fid)),
        Plot(esta_opt1, Coordinates(final_times, esta_fid1)),
        # legend
    )
    return ax
end
save_fidelity(n, final_times) = pgfsave(plot_folder * "fidelity_$(n).pdf", plot_fidelity(n, final_times))
# Plotting the sensitivity with respect to the modulation noise
"""
plot_sensitivity_mn(n::Int64, final_times) -> Axis
Calculate the sensitivity of the eSTA and STA for a given number of particles with respect to the modulation noise.
It returns an Axis object that can be saved with the pgfsave function.
"""
function plot_sensitivity_mn(n::Int64, final_times)
    # Calculate the fidelity
    cp = ControlParameterFull(0.1, n)
    esta_mn, sta_mn = robustnesses_mn(cp, final_times)
    # Plotting with PGFPlotsX
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles[1]}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles[2]}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "S", legend_pos = "north east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_mn)),
        Plot(sta_opt, Coordinates(final_times, sta_mn)),
        # legend
    )
    return ax
end
save_sensitivity_mn(n, final_times) = pgfsave(plot_folder * "sensitivity_mn_$(n).pdf", plot_sensitivity_mn(n, final_times))
# Plotting the sensitivity with respect to the timing noise
"""
plot_sensitivity_tn(n::Int64, final_times) -> Axis
Calculate the sensitivity of the eSTA and STA for a given number of particles with respect to the timing noise.
It returns an Axis object that can be saved with the pgfsave function.
"""
function plot_sensitivity_tn(n::Int64, final_times, eps::Float64)
    # Calculate the fidelity
    cp = ControlParameterFull(0.1, n)
    esta_tn, sta_tn = robustnesses_tn(cp, final_times, eps)
    # Plotting with PGFPlotsX
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles[1]}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles[2]}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "S", legend_pos = "north east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_tn)),
        Plot(sta_opt, Coordinates(final_times, sta_tn)),
        # legend
    )
    return ax
end
# Plotting the squeezing
"""
plot_squeezing(n::Int64, final_times) -> Axis
Calculate the squeezing of the eSTA and STA for a given number of particles.
It returns an Axis object that can be saved with the pgfsave function.
"""
save_sensitivity_tn(n, final_times, eps) = pgfsave(plot_folder * "sensitivity_tn_$(n).pdf", plot_sensitivity_tn(n, final_times, eps))
function plot_squeezing(n::Int64, final_times)
    # Calculate the squeezing of the eSTA and STA at final time
    cp = ControlParameterFull(0.1, n)
    esta_sq, sta_sq, id_sq = squeezings(cp, final_times)
    # Plotting with PGFPlotsX
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles[1]}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles[2]}
    id_opt = @pgf {color = colors.green, line_width = 1.5, style = styles[2]}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA"), LegendEntry("Ideal")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "\$\\xi_s^2\$", legend_pos = "north east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_sq)),
        Plot(sta_opt, Coordinates(final_times, sta_sq)),
        Plot(id_opt, Coordinates(final_times, id_sq)),
        # legend
    )
    return ax
end
save_squeezing(n, final_times) = pgfsave(plot_folder * "squeezing_$(n).pdf", plot_squeezing(n, final_times))
# Plotting the sensitivity with respect to the white noise
"""
plot_sensitivity_wn(n::Int64, final_times) -> Axis
Calculate the sensitivity of the eSTA and STA for a given number of particles with respect to the white noise.
It returns an Axis object that can be saved with the pgfsave function.
"""
function plot_sensitivity_wn(n::Int64, final_times)
    # Calculate the fidelity
    esta_wn, sta_wn = sensitivities(n, final_times)
    # Plotting with PGFPlotsX
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles[1]}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles[2]}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "S", legend_pos = "north east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_wn)),
        Plot(sta_opt, Coordinates(final_times, sta_wn)),
        # legend
    )
    return ax
end
# Here I will sum up all the functions that I have defined above and save the plots in the folder `gfx`
save_sensitivity_wn(n, final_times) = pgfsave(plot_folder * "sensitivity_wn_$(n).pdf", plot_sensitivity_wn(n, final_times))

# Here I will save the plots
final_times = range(0.15, 0.8, length=100) |> collect
# Plotting with PGFPlotsX
using PGFPlotsX
# Shared options for all plots
plot_folder = "./gfx/"
cp = ControlParameterFull(0.1, 10)

include("./src/ConstantQuantities.jl")

save_fidelity(10, final_times)
save_sensitivity_mn(10, final_times)
save_sensitivity_tn(10, final_times, 0.0001)
save_squeezing(10, final_times)
save_sensitivity_wn(10, final_times)
