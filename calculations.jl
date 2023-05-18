include("./src/ConstantQuantities.jl")
"""
fidelities(np::Int64, final_times::AbstractVector)
returns  all the feature I am interested in plotting for a given number of particles and a vector of final times.
The returned values are a list of vectors, as it follows:
1. fid_esta: fidelity of the ESTA protocol
2. fid_esta1: fidelity of the ESTA protocol with only one correction
3. fid_esta_cont: fidelity of the ESTA protocol with continuous corrections
4. fid_sta: fidelity of the STA protocol
5. fid_ad: fidelity of the adiabatic protocol
6. senstn_esta: sensitivity of the ESTA protocol
7. senstn_sta: sensitivity of the STA protocol
8. sensmn_esta: sensitivity of the ESTA protocol
9. sensmn_sta: sensitivity of the STA protocol
10. final_times: the same vector passed as argument
"""
function fidelities(np::Int64, final_times::AbstractVector)
    cparam = ControlParameterFull(np, final_times[1])
    qts = ConstantQuantities(cparam)
    fid_esta = zeros(length(final_times))
    fid_sta = zeros(length(final_times))
    fid_esta1 = zeros(length(final_times))
    fid_esta_cont = zeros(length(final_times))
    fid_ad = zeros(length(final_times))
    senstn_esta = zeros(length(final_times))
    senstn_sta = zeros(length(final_times))
    sensmn_esta = zeros(length(final_times))
    sensmn_sta = zeros(length(final_times))
    Threads.@threads for index in 1:length(final_times)
        tf = final_times[index]
        cparam = cp_time(cparam, tf)
        cparam_cont = ControlParameterInt(np, tf)
        corrs = corrections(cparam)
        corrs1 = corrections(cparam; nlambda=1)
        corrs_cont = corrections(cparam_cont)
        esta(t) = Λ_esta(t, cparam, corrs)
        esta1(t) = Λ_esta(t, cparam, corrs1)
        esta_cont(t) = Λ_esta(t, cparam_cont, corrs_cont)
        sta(t) = Λ_sta(t, cparam)
        ad(t) = Λ_ad(t, cparam)
        fid_esta[index] = fidelity(cparam, qts, esta)
        fid_esta1[index] = fidelity(cparam, qts, esta1)
        fid_sta[index] = fidelity(cparam, qts, sta)
        fid_esta_cont[index] = fidelity(cparam_cont, qts, esta_cont)
        fid_ad[index] = fidelity(cparam, qts, ad)
        senstn_esta[index] = robustness_tn(cparam, qts, esta, 1e-7)
        senstn_sta[index] = robustness_tn(cparam, qts, sta, 1e-7)
        sensmn_esta[index] = robustness_mn(cparam, qts, esta, 1e-7)
        sensmn_sta[index] = robustness_mn(cparam, qts, sta, 1e-7)
    end
    return fid_esta, fid_esta1, fid_esta_cont, fid_sta, fid_ad, senstn_esta, senstn_sta, sensmn_esta, sensmn_sta, final_times ./ pi
end
final_times = range(0.01pi, 0.2pi, length=100)
ft_10 = fidelities(10, final_times)
ft_30 = fidelities(30, final_times)

using PGFPlotsX
PGFPlotsX.enable_interactive(false)
colors = (
    black=colorant"#000000", # STA	
    red=colorant"#FF0000", # eSTA Full Hamiltonian with Hessian 
    blue=colorant"#0000FF", # eSTA Intermediate Hamiltonian with Hessian
    yellow=colorant"#FFa000", # eSTA Full Hamiltonian with original version
    green=colorant"#008000"# eSTA Intermediate Hamiltonian original version
)
styles = (
    solid="solid",                         # eSTA Intermediate Hamiltonian with Hessian
    dot_dash="dash pattern={on 4pt off 1pt on 1pt off 1pt}",#STA
    ldash="dash pattern={on 2pt off 2pt}",# eSTA Full Hamiltonian with Hessian 
    dash=" dashed",                        # eSTA Full Hamiltonian with original version
    dot="dash pattern={on 1pt off 1pt}",  # eSTA Intermediate Hamiltonian original version
)
common_style = @pgf {group_size = "1 by 2", vertical_sep = "0pt", xticklabels_at = "edge bottom", xlabels_at = "edge bottom"}
esta_opt = @pgf {color = colors.red, line_width = 1, style = styles.solid}
esta1_opt = @pgf {color = colors.blue, line_width = 1, style = styles.dot}
esta_cont_opt = @pgf {color = colors.yellow, line_width = 1, style = styles.ldash}
sta_opt = @pgf{color = colors.black, line_width = 1, style = styles.dash}
ad_opt = @pgf {color = colors.green, line_width = 1, style = styles.dot_dash}

# Here I will plot all the fidelities
fidelities_plot = @pgf GroupPlot(
    {
        group_style = common_style,
        enlarge_x_limits = "false",
        ymin = 0.7, ymax = 1.02,
        ylabel = "F", xlabel = "\$\\mathrm{t_f/t_R}\$",
        ticklabel_style = "/pgf/number format/fixed",
    },
    {}, # Plot for 10 particles
    Plot(esta_opt, Table(ft_10[10], ft_10[1])),
    Plot(esta1_opt, Table(ft_10[10], ft_10[2])),
    Plot(esta_cont_opt, Table(ft_10[10], ft_10[3])),
    Plot(sta_opt, Table(ft_10[10], ft_10[4])),
    Plot(ad_opt, Table(ft_10[10], ft_10[5])),
    ["\\node[anchor=north west] at (rel axis cs:0.7,0.2) {N = 10};"],
    {}, # Plot for 30 particles
    Plot(esta_opt, Table(ft_30[10], ft_30[1])),
    Plot(esta1_opt, Table(ft_30[10], ft_30[2])),
    Plot(esta_cont_opt, Table(ft_30[10], ft_30[3])),
    Plot(sta_opt, Table(ft_30[10], ft_30[4])),
    Plot(ad_opt, Table(ft_30[10], ft_30[5])),
    ["\\node[anchor=north west] at (rel axis cs:0.7,0.2) {N = 30};"],
)
display("gfx/fid_plot.pdf", fidelities_plot)

# Here I will plot all the sensitivities, starting with the TN
tn_plot = @pgf GroupPlot(
    {
        group_style = common_style,
        enlarge_x_limits = "false",
        xmin = ft_10[10][9], xmax = ft_10[10][end],
        ymin = 0.0, ymax = 0.22,
        ylabel = "\$S_{t}\$", xlabel = "\$\\mathrm{t_f/t_R}\$",
        ticklabel_style = "/pgf/number format/fixed",
        max_space_between_ticks = "40pt",
    },
    {}, # Plot for 10 particles
    Plot(esta_opt, Table(ft_10[10], ft_10[6])),
    Plot(sta_opt, Table(ft_10[10], ft_10[7])),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 10};"],
    {}, # Plot for 30 particles
    Plot(esta_opt, Table(ft_30[10], ft_30[6])),
    Plot(sta_opt, Table(ft_30[10], ft_30[7])),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 30};"],
)
display("gfx/tn_plot.pdf", tn_plot)

# Now it's time for the MN
mn_plot = @pgf GroupPlot(
    {
        group_style = common_style,
        enlarge_x_limits = "false",
        xmin = ft_10[10][9], xmax = ft_10[10][end],
        ymin = 0.0, ymax = 0.18,
        ylabel = "\$S_{m}\$", xlabel = "\$\\mathrm{t_f/t_R}\$",
        ticklabel_style = "/pgf/number format/fixed",
        max_space_between_ticks = "40pt",
    },
    {}, # Plot for 10 particles
    Plot(esta_opt, Table(ft_10[10], ft_10[8])),
    Plot(sta_opt, Table(ft_10[10], ft_10[9])),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 10};"],
    {}, # Plot for 30 particles
    Plot(esta_opt, Table(ft_30[10], ft_30[8])),
    Plot(sta_opt, Table(ft_30[10], ft_30[9])),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 30};"],
)
display("gfx/mn_plot.pdf", mn_plot)

# And now for the effectiveness, η evaluted as ( (1-F^2) + (S_t^2 + S_m^2) )^(1/2)
# First for 10 particles
η_esta10 = sqrt.((1 .- ft_10[1] .^ 2) .+ (ft_10[6] .^ 2 .+ ft_10[8] .^ 2))
η_sta10 = sqrt.((1 .- ft_10[4] .^ 2) .+ (ft_10[7] .^ 2 .+ ft_10[9] .^ 2))
# And now for 30 particles
η_esta30 = sqrt.((1 .- ft_30[1] .^ 2) .+ (ft_30[6] .^ 2 .+ ft_30[8] .^ 2))
η_sta30 = sqrt.((1 .- ft_30[4] .^ 2) .+ (ft_30[7] .^ 2 .+ ft_30[9] .^ 2))
eff_plot = @pgf GroupPlot(
    {
        group_style = common_style,
        enlargelimits = "false",
        # xmin = ft_10[10][9], xmax = ft_10[10][end],
        ymax = 1.05,
        ylabel = "\$\\eta\$", xlabel = "\$\\mathrm{t_f/t_R}\$",
        ticklabel_style = "/pgf/number format/fixed",
        max_space_between_ticks = "40pt",
    },
    {}, # Plot for 10 particles
    Plot(esta_opt, Table(ft_10[10], η_esta10)),
    Plot(sta_opt, Table(ft_10[10], η_sta10)),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 10};"],
    {}, # Plot for 30 particles
    Plot(esta_opt, Table(ft_30[10], η_esta30)),
    Plot(sta_opt, Table(ft_30[10], η_sta30)),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 30};"],
)
display("gfx/eff_plot.pdf", eff_plot)

# What I need to do now is to find the time for which the fidelity reaches 0.99 for both eSTA and STA.
# I will do this for both 10 and 30 particles.
function fidelity_time(array, np::Int64)
    t_esta = array[end][findfirst(x -> x .> 0.98, array[1])]
    t_sta = array[end][findfirst(x -> x .> 0.98, array[4])]
    cp_esta = ControlParameterFull(np, t_esta * pi)
    cp_sta = ControlParameterFull(np, t_sta * pi)
    cp_ad = ControlParameterFull(np, t_sta * 5.0 * pi)
    qts = ConstantQuantities(cp_esta)
    corrs = corrections(cp_esta)
    esta(t) = Λ_esta(t, cp_esta, corrs)
    sta(t) = Λ_sta(t, cp_sta)
    ad(t) = Λ_ad(t, cp_ad)
    tab_esta = @pgf Table(range(0, stop=t_esta, length=100), fidelity(cp_esta, qts, esta, 100))
    tab_sta = @pgf Table(range(0, stop=t_sta, length=100), fidelity(cp_sta, qts, sta, 100))
    tab_ad = @pgf Table(range(0, stop=5.0 * t_sta, length=100), fidelity(cp_ad, qts, ad, 100))
    return @pgf Axis({
            xlabel = "\$\\mathrm{t/t_R}\$",
            ylabel = "F",
            ymin = 0.6, ymax = 1.0,
            xmin = 0.0, xmax = t_sta,
            ticklabel_style = "/pgf/number format/fixed",
        },
        Plot(esta_opt, tab_esta),
        Plot(sta_opt, tab_sta),
        Plot(ad_opt, tab_ad),
    )
end
fidelity_time_10 = fidelity_time(ft_10, 10)
display("gfx/fidelity_time.pdf", fidelity_time_10)

tf = 0.1pi
nsteps = 100
t = range(0, stop=tf, length=nsteps)
cparam = ControlParameterFull(10, tf)
qts = ConstantQuantities(cparam)
corrs = corrections(cparam)
esta(t) = Λ_esta(t, cparam, corrs)
sta(t) = Λ_sta(t, cparam)
ad(t) = Λ_ad(t, cparam)
fid_esta = fidelity(cparam, qts, esta, nsteps)
fid_sta = fidelity(cparam, qts, sta, nsteps)
fid_ad = fidelity(cparam, qts, ad, nsteps)
α_esta = alpha(cparam, qts, esta, nsteps)
α_sta = alpha(cparam, qts, sta, nsteps)
α_ad = alpha(cparam, qts, ad, nsteps)
ξN_esta = squeezing(cparam, qts, esta, nsteps)
ξN_sta = squeezing(cparam, qts, sta, nsteps)
ξN_ad = squeezing(cparam, qts, ad, nsteps)
todecibel(x) = 10 * log10(x)
ξs_esta = todecibel.(ξN_esta .^ 2 ./ (α_esta .^ 2))
ξs_sta = todecibel.(ξN_sta .^ 2 ./ (α_sta .^ 2))
ξs_ad = todecibel.(ξN_ad .^ 2 ./ (α_ad .^ 2))

evo_plot = @pgf GroupPlot({
        group_style = {
            group_size = "1 by 4",
            xlabels_at = "edge bottom",
            xticklabels_at = "edge bottom",
            vertical_sep = "10pt",
        },
        enlarge_y_limits = "0.01",
        enlarge_x_limits = "false",
        ticklabel_style = "/pgf/number format/fixed",
        max_space_between_ticks = "40pt",
        height = "5cm",
        width = "10cm",
        ylabel_style = "at ={(rel axis cs: -0.08,0.5)}",
    },
    {
        ylabel = "\$\\Lambda\$(t)"
    },
    Plot(esta_opt, Table(t ./ pi, esta.(t))),
    Plot(sta_opt, Table(t ./ pi, sta.(t))),
    Plot(ad_opt, Table(t ./ pi, ad.(t))),
    {
        ylabel = "F",
    },
    Plot(esta_opt, Table(t ./ pi, fid_esta)),
    Plot(sta_opt, Table(t ./ pi, fid_sta)),
    Plot(ad_opt, Table(t ./ pi, fid_ad)),
    # ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 10};" ],
    # Number squeezing
    {
        ylabel = "\$\\xi_N^2\$ "
    },
    Plot(esta_opt, Table(t ./ pi, ξN_esta)),
    Plot(sta_opt, Table(t ./ pi, ξN_sta)),
    Plot(ad_opt, Table(t ./ pi, ξN_ad)),
    # Coherent squeezing
    {
        ylabel = "\$\\xi_s^2\$ [dB]",
        xlabel = "\$\\mathrm{t/t_R}\$"
    },
    Plot(esta_opt, Table(t ./ pi, ξs_esta)),
    Plot(sta_opt, Table(t ./ pi, ξs_sta)),
    Plot(ad_opt, Table(t ./ pi, ξs_ad)),
)
display("gfx/evo_plot.pdf", evo_plot)

#==
`fid_plot.pdf` is the plot of the fidelities for 10 and 30 particles. In particular, we have
- red line: eSTA with 5 λ applied to the discrete Hamiltonian
- blue line: eSTA with 1 λ applied to the discrete Hamiltonian
- yellow line: eSTA with 5 λ applied to the continuous Hamiltonian
- black line: STA protocol
- green line: adiabatic protocol 

`tn_plot.pdf` is the plot for the sensitivities with respect to the time noise for 10 and 30 particles.
`mn_plot.pdf` is the plot for the sensitivities with respect to the modulation noise for 10 and 30 particles.
`eff_plot.pdf` is the plot for the effectiveness of the protocols for 10 and 30 particles. It is evaluated as ( (1-F^2) + (S_t^2 + S_m^2) )^(1/2)
`fidelity_time.pdf` is the plot of the fidelity that shows how lont it takes to reach a given fidelity for 10 particles. In this case I picked F = 0.99.
`evo_plot.pdf` is the plot whit all the time evolution of the relevant quantities for 10 particles and for final time 0.1 t_rabi. In particular, we have
- first row: time evolution of the control functions Λ(t)
- second row: time evolution of the fidelity F(t)
- third row: time evolution of the number squeezing ξ_N^2(t)
- fourth row: time evolution of the coherent squeezing ξ_s^2(t) in dB
==#
