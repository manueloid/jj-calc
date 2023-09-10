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
6. tn_esta: sensitivity of the ESTA protocol
7. tn_sta: sensitivity of the STA protocol
8. tn_esta1: sensitivity of the ESTA protocol with only one correction
9. tn_esta_cont: sensitivity of the ESTA protocol with continuous corrections
10. mn_esta: sensitivity of the ESTA protocol
11. mn_sta: sensitivity of the STA protocol
12. mn_esta1: sensitivity of the ESTA protocol with only one correction
13. mn_esta_cont: sensitivity of the ESTA protocol with continuous corrections
14. pl_esta: particle loss of the ESTA protocol
15. pl_esta: particle loss of the ESTA protocol
16. pl_sta: particle loss of the STA protocol
17. pl_esta1: particle loss of the ESTA protocol with only one correction
18. tfs: final times 

"""
function fidelities(np::Int64, final_times::AbstractVector)
    cparam = ControlParameterFull(np, final_times[1])
    qts = ConstantQuantities(cparam)
    fid_esta = zeros(length(final_times))
    fid_sta = zeros(length(final_times))
    fid_esta1 = zeros(length(final_times))
    fid_esta_cont = zeros(length(final_times))
    fid_ad = zeros(length(final_times))
    tn_esta = zeros(length(final_times))
    tn_sta = zeros(length(final_times))
    tn_esta1 = zeros(length(final_times))
    tn_esta_cont = zeros(length(final_times))
    mn_esta = zeros(length(final_times))
    mn_sta = zeros(length(final_times))
    mn_esta1 = zeros(length(final_times))
    mn_esta_cont = zeros(length(final_times))
    pl_esta = zeros(length(final_times))
    pl_sta = zeros(length(final_times))
    pl_esta1 = zeros(length(final_times))
    pl_esta_cont = zeros(length(final_times))
    tfs = zeros(length(final_times))
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
        fid_esta_cont[index] = fidelity(cparam_cont, qts, esta_cont)
        fid_sta[index] = fidelity(cparam, qts, sta)
        fid_ad[index] = fidelity(cparam, qts, ad)
        tn_esta[index] = robustness_tn(cparam, qts, esta, 1e-7)
        tn_esta1[index] = robustness_tn(cparam, qts, esta1, 1e-7)
        tn_esta_cont[index] = robustness_tn(cparam_cont, qts, esta_cont, 1e-7)
        tn_sta[index] = robustness_tn(cparam, qts, sta, 1e-7)
        mn_esta[index] = robustness_mn(cparam, qts, esta, 1e-7)
        mn_esta1[index] = robustness_mn(cparam, qts, esta1, 1e-7)
        mn_esta_cont[index] = robustness_mn(cparam_cont, qts, esta_cont, 1e-7)
        mn_sta[index] = robustness_mn(cparam, qts, sta, 1e-7)
        pl_esta[index] = sensitivity(cparam, esta)
        pl_esta1[index] = sensitivity(cparam, esta1)
        pl_esta_cont[index] = sensitivity(cparam, esta_cont)
        pl_sta[index] = sensitivity(cparam, sta)
    end
    return fid_esta, fid_esta1, fid_esta_cont, fid_sta, fid_ad, tn_esta, tn_sta, tn_esta1, tn_esta_cont, mn_esta, mn_sta, mn_esta1, mn_esta_cont, pl_esta, pl_esta1, pl_esta_cont, pl_sta, final_times ./ pi
end
final_times = range(0.01pi, 0.2pi, length=150)
ft_10 = fidelities(10, final_times)
ft_100 = fidelities(100, final_times)

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
        ylabel = "\$F\$", xlabel = "\$t_f/t_R\$",
        ticklabel_style = "/pgf/number format/fixed",
    },
    {}, # Plot for 10 particles
    Plot(esta_opt, Table(ft_10[end], ft_10[1])),
    Plot(esta1_opt, Table(ft_10[end], ft_10[2])),
    Plot(esta_cont_opt, Table(ft_10[end], ft_10[3])),
    Plot(sta_opt, Table(ft_10[end], ft_10[4])),
    Plot(ad_opt, Table(ft_10[end], ft_10[5])),
    ["\\node[anchor=north west] at (rel axis cs:0.7,0.2) {N = 10};"],
    {}, # Plot for 100 particles
    Plot(esta_opt, Table(ft_100[end], ft_100[1])),
    Plot(esta1_opt, Table(ft_100[end], ft_100[2])),
    Plot(esta_cont_opt, Table(ft_100[end], ft_100[3])),
    Plot(sta_opt, Table(ft_100[end], ft_100[4])),
    Plot(ad_opt, Table(ft_100[end], ft_100[5])),
    ["\\node[anchor=north west] at (rel axis cs:0.7,0.2) {N = 100};"],
)
tp = @pgf TikzPicture({font = "\\fontsize{11}{11}\\selectfont"})
push!(tp, fidelities_plot)
# display("gfx/fid_plot.pdf", fidelities_plot)
display("/home/manuel/pCloudDrive/PhD/QuantumThermo/JJeSTA/Paper/fig_fid.pdf", tp)
# And now for the effectiveness, η evaluted as ( (1-F^2) + (S_t^2 + S_m^2) )^(1/2)
# First for 10 particles
η_esta10 = sqrt.((1 .- ft_10[1] .^ 2) .+ (ft_10[6] .^ 2 .+ ft_10[10] .^ 2))
η_sta10 = sqrt.((1 .- ft_10[4] .^ 2) .+ (ft_10[7] .^ 2 .+ ft_10[11] .^ 2))
η_esta1_10 = sqrt.((1 .- ft_10[2] .^ 2) .+ (ft_10[8] .^ 2 .+ ft_10[12] .^ 2))
η_esta_cont10 = sqrt.((1 .- ft_10[3] .^ 2) .+ (ft_10[9] .^ 2 .+ ft_10[13] .^ 2))
# And now for 100 particles
η_esta100 = sqrt.((1 .- ft_100[1] .^ 2) .+ (ft_100[6] .^ 2 .+ ft_100[10] .^ 2))
η_sta100 = sqrt.((1 .- ft_100[4] .^ 2) .+ (ft_100[7] .^ 2 .+ ft_100[11] .^ 2))
η_esta1_100 = sqrt.((1 .- ft_100[2] .^ 2) .+ (ft_100[8] .^ 2 .+ ft_100[12] .^ 2))
η_esta_cont100 = sqrt.((1 .- ft_100[3] .^ 2) .+ (ft_100[9] .^ 2 .+ ft_100[13] .^ 2))
# Plot where all the robustness parameters are plotted together. The layout will be 3 x 2.
robustness_plot = @pgf GroupPlot(
    {
        group_style = {
            group_size = "3 by 2",
            vertical_sep = "0.0cm",
            horizontal_sep = "1.4cm",
            xlabels_at = "edge bottom",
            xticklabels_at = "edge bottom",
        },
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.02",
        xmin = ft_10[end][9], xmax = ft_10[end][75],
        width = "0.5\\textwidth",
        ticklabel_style = "/pgf/number format/fixed",
        xlabel = "\$t_f/t_R\$",
    },
    #
    # Plots for 10 particles,
    #
    {
        ylabel = "\$S_{m}\$",
    }, # robustness against modulation noise
    Plot(esta_opt, Table(ft_10[end], ft_10[10])),
    Plot(sta_opt, Table(ft_10[end], ft_10[11])),
    Plot(esta1_opt, Table(ft_10[end], ft_10[12])),
    Plot(esta_cont_opt, Table(ft_10[end], ft_10[13])),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 10};"],
    {
        ylabel = "\$S_{t}\$",
    }, # robustness against time noise
    Plot(esta_opt, Table(ft_10[end], ft_10[6])),
    Plot(sta_opt, Table(ft_10[end], ft_10[7])),
    Plot(esta1_opt, Table(ft_10[end], ft_10[8])),
    Plot(esta_cont_opt, Table(ft_10[end], ft_10[9])),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 10};"],
    {
        ylabel = "\$\\eta\$",
    }, # effectiveness Plots
    Plot(esta_opt, Table(ft_10[end], η_esta10)),
    Plot(sta_opt, Table(ft_10[end], η_sta10)),
    Plot(esta1_opt, Table(ft_10[end], η_esta1_10)),
    Plot(esta_cont_opt, Table(ft_10[end], η_esta_cont10)),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 10};"],
    #
    # Plots for 100 particles,
    # robustness against modulation noise
    {
        ylabel = "\$S_{m}\$",
        "extra_description/.code = { \\node at (rel axis cs: 0,-0.2) {(b)};}",
    },
    Plot(esta_opt, Table(ft_100[end], ft_100[10])),
    Plot(sta_opt, Table(ft_100[end], ft_100[11])),
    Plot(esta1_opt, Table(ft_100[end], ft_100[12])),
    Plot(esta_cont_opt, Table(ft_100[end], ft_100[13])),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 100};"],
    #  robustness against time noise
    {
        ylabel = "\$S_{t}\$",
        "extra_description/.code = { \\node at (rel axis cs: 0,-0.2) {(a)};}",
    },
    Plot(esta_opt, Table(ft_100[end], ft_100[6])),
    Plot(sta_opt, Table(ft_100[end], ft_100[7])),
    Plot(esta1_opt, Table(ft_100[end], ft_100[8])),
    Plot(esta_cont_opt, Table(ft_100[end], ft_100[9])),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 100};"],
    # effectiveness Plots
    {
        ylabel = "\$\\eta\$",
        "extra_description/.code = { \\node at (rel axis cs: 0,-0.2) {(c)};}",
    },
    Plot(esta_opt, Table(ft_100[end], η_esta100)),
    Plot(sta_opt, Table(ft_100[end], η_sta100)),
    Plot(esta1_opt, Table(ft_100[end], η_esta1_100)),
    Plot(esta_cont_opt, Table(ft_100[end], η_esta_cont100)),
    ["\\node[anchor = north east] at (rel axis cs: 0.9,0.9) {N = 100};"],
)
tp = @pgf TikzPicture({font = "\\fontsize{8}{8}\\selectfont"})
push!(tp, robustness_plot)
# display("gfx/robustness_plot.pdf", tp)
display("/home/manuel/pCloudDrive/PhD/QuantumThermo/JJeSTA/Paper/robustness_plot.pdf", tp)

# What I need to do now is to find the time for which the fidelity reaches 0.99 for both eSTA and STA.
# I will do this for both 10 and 100 particles.
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
            xlabel = "\$t/t_R\$",
            ylabel = "\$F\$",
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
        try_min_ticks = 3,
        xtick_distance = 0.05,
        width = "0.5\\textwidth",
        height = "0.25\\textwidth",
        ylabel_style = "at ={(rel axis cs: -0.12,0.5)}",
    },
    {
        ylabel = "\$\\Lambda\$(t)",
        "extra_description/.code = { \\node at (rel axis cs: -0.16,0.94) {(a)};}",
    },
    Plot(esta_opt, Table(t ./ pi, esta.(t))),
    Plot(sta_opt, Table(t ./ pi, sta.(t))),
    Plot(ad_opt, Table(t ./ pi, ad.(t))),
    # Number squeezing
    {
        ylabel = "\$\\xi_N^2\$ ",
        "extra_description/.code = { \\node at (rel axis cs: -0.16,0.94) {(b)};}",
    },
    Plot(esta_opt, Table(t ./ pi, ξN_esta)),
    Plot(sta_opt, Table(t ./ pi, ξN_sta)),
    Plot(ad_opt, Table(t ./ pi, ξN_ad)),
    # Coherent squeezing
    {
        ylabel = "\$\\xi_s^2\$ ",
        "extra_description/.code = { \\node at (rel axis cs: -0.16,0.94) {(c)};}",
    },
    Plot(esta_opt, Table(t ./ pi, ξs_esta)),
    Plot(sta_opt, Table(t ./ pi, ξs_sta)),
    Plot(ad_opt, Table(t ./ pi, ξs_ad)),
    # Fidelity
    {
        ylabel = "\$F\$",
        xlabel = "\$t/t_R\$",
        ymin = 0.6, ymax = 1.0,
        "extra_description/.code = { \\node at (rel axis cs: -0.16,0.94) {(d)};}",
    },
    Plot(esta_opt, Table(t ./ pi, fid_esta)),
    Plot(sta_opt, Table(t ./ pi, fid_sta)),
    Plot(ad_opt, Table(t ./ pi, fid_ad)),
)
tp = @pgf TikzPicture({font = "\\fontsize{8}{8}\\selectfont"})
push!(tp, evo_plot)
# display("gfx/evo_plot.pdf", tp)
display("/home/manuel/pCloudDrive/PhD/QuantumThermo/JJeSTA/Paper/fig_evo.pdf", tp)

#==
`fid_plot.pdf` is the plot of the fidelities for 10 and 100 particles. In particular, we have
- red line: eSTA with 5 λ applied to the discrete Hamiltonian
- blue line: eSTA with 1 λ applied to the discrete Hamiltonian
- yellow line: eSTA with 5 λ applied to the continuous Hamiltonian
- black line: STA protocol
- green line: adiabatic protocol 

`tn_plot.pdf` is the plot for the sensitivities with respect to the time noise for 10 and 100 particles.
`mn_plot.pdf` is the plot for the sensitivities with respect to the modulation noise for 10 and 100 particles.
`eff_plot.pdf` is the plot for the effectiveness of the protocols for 10 and 100 particles. It is evaluated as ( (1-F^2) + (S_t^2 + S_m^2) )^(1/2)
`fidelity_time.pdf` is the plot of the fidelity that shows how lont it takes to reach a given fidelity for 10 particles. In this case I picked F = 0.99.
`evo_plot.pdf` is the plot whit all the time evolution of the relevant quantities for 10 particles and for final time 0.1 t_rabi. In particular, we have
- first row: time evolution of the control functions Λ(t)
- second row:time evolution of the coherent squeezing ξ_s^2(t) in dB
- third row: time evolution of the number squeezing ξ_N^2(t)
- fourth row:  time evolution of the fidelity F(t)
==#
