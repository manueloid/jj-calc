include("./src/ConstantQuantities.jl")
include("./src/PlotObj.jl")

cp = ControlParameterFull(0.1, 50) # I will use 10 particles to start with
final_times = range(0.05, 0.5, length=100) |> collect

# Plotting with PGFPlotsX
using PGFPlotsX
# Shared options for all plots
plot_folder = "./gfx/"
"""
plot_fidelity(n::Int64, final_times) -> Axis
Calculate the fidelity of the eSTA and STA for a given number of particles.
It returns an Axis object that can be saved with the pgfsave function.
"""
function plot_fidelity(n::Int64, final_times)
    # Calculate the fidelity
    cp = ControlParameterFull(0.1, n)
    esta_fid, sta_fid = fidelities(cp, final_times)
    # Plotting with PGFPlotsX
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles.dash}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles.solid}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "F", legend_pos = "south east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_fid)),
        Plot(sta_opt, Coordinates(final_times, sta_fid)),
        legend
    )
    return ax
end
plot_fidelity(10, final_times)

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
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles.dash}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles.solid}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "S", legend_pos = "north east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_mn)),
        Plot(sta_opt, Coordinates(final_times, sta_mn)),
        legend
    )
    return ax
end

plot_sensitivity_mn(10, final_times)

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
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles.dash}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles.solid}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "S", legend_pos = "north east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_tn)),
        Plot(sta_opt, Coordinates(final_times, sta_tn)),
        legend
    )
    return ax
end

plot_sensitivity_tn(10, final_times, 1e-6)

# Plotting the squeezing
"""
plot_squeezing(n::Int64, final_times) -> Axis
Calculate the squeezing of the eSTA and STA for a given number of particles.
It returns an Axis object that can be saved with the pgfsave function.
"""
function plot_squeezing(n::Int64, final_times)
    # Calculate the squeezing of the eSTA and STA at final time
    cp = ControlParameterFull(0.1, n)
    esta_sq, sta_sq = squeezings(cp, final_times)
    # Plotting with PGFPlotsX
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles.dash}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles.solid}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "S", legend_pos = "north east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_sq)),
        Plot(sta_opt, Coordinates(final_times, sta_sq)),
        legend
    )
    return ax
end

plot_squeezing(10, final_times)

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
    esta_opt = @pgf {color = colors.red, line_width = 1.5, style = styles.dash}
    sta_opt = @pgf {color = colors.black, line_width = 1.5, style = styles.solid}
    legend = @pgf [LegendEntry("eSTA"), LegendEntry("STA")]
    ax = @pgf Axis(
        # Setting the plot axis
        {xlabel = raw"$t_f/\tau$", ylabel = "S", legend_pos = "north east"},
        # Plotting the fidelity
        Plot(esta_opt, Coordinates(final_times, esta_wn)),
        Plot(sta_opt, Coordinates(final_times, sta_wn)),
        legend
    )
    return ax
end

plot_sensitivity_wn(10, final_times)

fidelity
