using JosephsonSTA
include("./src/ConstantQuantities.jl")
ts = range(0.05, 1.0, length=100)
# esta10, sta10 = sensitivities(10, ts)
# esta100, sta100 = sensitivities(100, ts)
esta, sta = sensitivities(10, ts; λs=1, cont=true)

using DelimitedFiles
writedlm("./data/ploss_10.dat", [ts esta10 sta10])
writedlm("./data/ploss_100.dat", [ts esta100 sta100])

data10 = readdlm("./data/ploss_10.dat")
data100 = readdlm("./data/ploss_100.dat")

using PGFPlotsX, Colors
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
sta_opt = @pgf{color = colors.black, line_width = 1, style = styles.dash}

# PGFPlotsX.enable_interactive("false")
ploss_plot = @pgf GroupPlot(
    {
        group_style = common_style,
        enlarge_x_limits = false,
        enlarge_y_limits = false,
        ylabel = raw"$P_l$",
    },
    {
        legend_pos = "north west",
    },
    # Plot(esta_opt, Table(data10[:, 1] ./ pi, data10[:, 2])),
    # Plot(sta_opt, Table(data10[:, 1] ./ pi, data10[:, 3])),
    Plot(esta_opt, Table(ts, esta))
    Plot(sta_opt, Table(ts, sta))
    Legend(["esta", "sta"]),
    {},
    # Plot(esta_opt, Table(data100[:, 1] ./ pi, data100[:, 2])),
    # Plot(sta_opt, Table(data100[:, 1] ./ pi, data100[:, 3])),
)
display("./gfx/ploss_group.pdf", ploss_plot)
