using CSV, DataFrames, PGFPlotsX, Colors
filter_val(df::DataFrame, col_name::Symbol, value) = filter(row -> row[col_name] == value, df) # filter by value
filter_val(df::DataFrame, col_name::Symbol, min, max) = filter(row -> row[col_name] > min && row[col_name] < max, df) # filter by range
filter_val(df::DataFrame, N::Int64, λ::Int64) = filter(row -> row.N == N && row.λ == λ, df)
df = CSV.read("./data/whole_data.dat", DataFrame)

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
ad_opt = @pgf {color = colors.green, line_width = 1, style = styles.dot_dash}
extra_opt = @pgf {color = colors.yellow, line_width = 1, style = styles.ldash}
feat(column::String, df::DataFrame) = Table(df.tf, df[!, column])

PGFPlotsX.enable_interactive(false)
@pgf gp_fid = GroupPlot(
    {
        xmin = 0.1, xmax = 0.6,
        ymin = 0.00, ymax = 0.2,
        group_style = common_style,
        xlabel = raw"$t_f$",
        ylabel = raw"$\mathcal{S}$"
    },
    {},
    Plot(esta_opt, feat("Mn_eSTA", df_10)),
    Plot(sta_opt, feat("Mn_STA", df_10)),
    Plot(ad_opt, feat("Tn_eSTA", df_10)),
    Plot(extra_opt, feat("Tn_STA", df_10)),
    {},
    Plot(esta_opt, feat("Mn_eSTA", df_30)),
    Plot(sta_opt, feat("Mn_STA", df_30)),
    Plot(ad_opt, feat("Tn_eSTA", df_30)),
    Plot(extra_opt, feat("Tn_STA", df_30)),
)
display("./data/fidelity.pdf", gp_fid)

effectiveness(fidelity::Float64, sensitivity::Float64) = sqrt((1 - fidelity)^2 + sensitivity^2)
broadcast(effectiveness, df_10.F_eSTA, df_10.Tn_eSTA)
