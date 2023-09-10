#=
# Plotting for the paper

In this file I am going to load the file for the fidelity and all the different sensitivities and plot them.
The data is stored into two different kind of files: one containing only the sensitivities against particle loss, and the other one containing all the rest.

## 1 Loading the data
The name of the first type of file is `./data/pl_<n>.csv`, where `<n>` is the number of particles.
The name of the second type of file is `./data/fid_<n>.csv`, where `<n>` is the number of particles.

The columns for the first type of file are:

1. pl_esta: particle loss of the esta protocol for 5 corrections applied to the original Hamiltonian
2. pl_esta1: partilcle loss of the ESTA protocol with only one correction applied to the original Hamiltonian
3. pl_esta_cont: particle loss of the ESTA protocol with 5 corrections applied to the continous approximation of the original Hamiltonian
4. pl_sta: particle loss of the STA protocol

The columns for the second type of file are:

1. fid_esta: fidelity of the ESTA protocol
2. fid_esta1: fidelity of the ESTA protocol with only one correction
3. fid_esta_cont: fidelity of the ESTA protocol with continuous corrections
4. fid_sta: fidelity of the STA protocol
5. fid_ad: fidelity of the adiabatic protocol
6. tn_esta: sensitivity against time error of the ESTA protocol with 5 corrections
7. tn_sta: sensitivity against time error of the STA protocol
8. tn_esta1: sensitivity against time error of the ESTA protocol with only one correction
9. tn_esta_cont: sensitivity against time error of the ESTA protocol with continuous corrections
10. mn_esta: sensitivity for the modulation noise of the ESTA protocol
11. mn_sta: sensitivity for the modulation noise of the STA protocol
12. mn_esta1: sensitivity for the modulation noise of the ESTA protocol with only one correction
13. mn_esta_cont: sensitivity for the modulation noise of the ESTA protocol with continuous corrections
14. final_times: the same vector passed as argument

I will load the files through the `DataFrames` and `CSV` packages and then plot using `PGFPlotsX`.
The available number of particles are 10 and 100
=#

using DataFrames, CSV, PGFPlotsX, Colors
PGFPlotsX.enable_interactive(false)
pl_10 = DataFrame(CSV.File("./data/pl_10.csv"))
pl_100 = DataFrame(CSV.File("./data/pl_100.csv"))
fid_10 = DataFrame(CSV.File("./data/fid_10.csv"))
fid_100 = DataFrame(CSV.File("./data/fid_100.csv"))
whole10 = innerjoin(pl_10, fid_10, on=:final_times)
whole100 = innerjoin(pl_100, fid_100, on=:final_times)

col(feature, df) = df[!, feature]
types = ["esta", "esta1", "esta_cont", "sta"]
features = ["tn", "mn", "pl"]
col(feature, type, df) = col(string(feature, "_", type), df)
# Plot all the modulation noise sensitivities

col("mn", "esta", whole10)

#=
## 2 Plotting style

Here I will define the style and the options for all the different plots.
I will copy them from other files, but the structure is always the same.

There will be a `line_style` variable that will store all the relevant styles regarding the lines.
There will also be a `color_style` variable that will store all the relevant styles regarding the colors.

The types of each variable are the one that PGFPlotsX accepts as argument, so they will be `String` for the style and `Color` for the colors.
The data structure will be a `Dict` so that it will be a little bit more flexible.

Once I defined all the different styles, I will then assign to each column the data I want.
=#

color_style = (
    black=colorant"#000000", # STA	
    red=colorant"#FF0000", # eSTA Full Hamiltonian with Hessian 5 corrections
    blue=colorant"#0000FF", # eSTA Intermediate Hamiltonian with one correction
    yellow=colorant"#FFa000", # eSTA intermediate Hamiltonian with Hessian 5 corrections
    green=colorant"#008000"# adiabatic
)
line_style = (
    solid="solid",
    dot_dash="dash pattern={on 4pt off 1pt on 1pt off 1pt}",
    ldash="dash pattern={on 2pt off 2pt}",
    dash=" dashed",
    dot="dash pattern={on 1pt off 1pt}",
)
esta_style = @pgf {color = color_style.red, line_width = 1, style = line_style.solid}
esta1_style = @pgf {color = color_style.blue, line_width = 1, style = line_style.dot}
esta_cont_style = @pgf {color = color_style.yellow, line_width = 1, style = line_style.ldash}
sta_style = @pgf{color = color_style.black, line_width = 1, style = line_style.dash}
ad_style = @pgf {color = color_style.green, line_width = 1, style = line_style.dot_dash}
styles = [esta_style, esta1_style, esta_cont_style, sta_style, ad_style]



#=
## 3 Calculating the effectiveness

The effectiveness is defined as $( (1-F^2) + (S_t^2 + S_m^2+ S_p) )^(1/2)$, where $F$ is the fidelity and $S_t$, $S_m$ and $S_p$ are the sensitivities against time error, modulation error and particle loss respectively.
I am going to define a new array starting from the columns of the dataframe I already have.
=#

eff10_data = [sqrt.(1 .- col("fid", type, whole10) .^ 2 .+ col("tn", type, whole10) .^ 2 .+ col("mn", type, whole10) .^ 2 .+ col("pl", type, whole10) .^ 2) for type in types]
eff100_data = [sqrt.(1 .- col("fid", type, whole100) .^ 2 .+ col("tn", type, whole100) .^ 2 .+ col("mn", type, whole100) .^ 2 .+ col("pl", type, whole100) .^ 2) for type in types]

#=
## 3 Single plots

In this section I am going to define a series of plots that I am only lately going to combine in a single one.
The plots for the fidelities have already been made, here I am going to focus on the diffent sensitivities.
The plots will be:
1. Sensitivity for error in modulation, for 10 and 100 particles
2. Sensitivity for error in time, for 10 and 100 particles
3. Sensitivity for particle loss, for 10 and 100 particles
4. The so called *effectiveness*, which is defined as $( (1-F^2) + (S_t^2 + S_m^2+ S_p) )^(1/2)$

Each of the plots will have four lines in it, respectively for the STA, ESTA with 5 corrections, ESTA with 1 correction and ESTA with continuous corrections.

In the end we will have a total of 8 plots, 4 for each number of particles.
Then we will decide if it will make more sense to put everything together in a single group plot, or it will make more sense to use an inset.

The plots will be an array of `Plot` objects, so that we can easily combine them later.
=#

mn10 = @pgf [Plot(style, Table(whole10.final_times, col("mn", type, whole10))) for (style, type) in zip(styles, types)]
mn100 = @pgf [Plot(style, Table(whole100.final_times, col("mn", type, whole100))) for (style, type) in zip(styles, types)]
tn10 = @pgf [Plot(style, Table(whole10.final_times, col("tn", type, whole10))) for (style, type) in zip(styles, types)]
tn100 = @pgf [Plot(style, Table(whole100.final_times, col("tn", type, whole100))) for (style, type) in zip(styles, types)]
pl10 = @pgf [Plot(style, Table(whole10.final_times, col("pl", type, whole10))) for (style, type) in zip(styles, types)]
pl100 = @pgf [Plot(style, Table(whole100.final_times, col("pl", type, whole100))) for (style, type) in zip(styles, types)]
eff10 = @pgf [Plot(style, Table(whole10.final_times, eff)) for (style, eff) in zip(styles, eff10_data)]
eff100 = @pgf [Plot(style, Table(whole100.final_times, eff)) for (style, eff) in zip(styles, eff100_data)]

#=
## 4 Group plots

Now it is time to plot all together in a single group plot.
I will try a 2x4 grid where all the quantities have their own plots, as well as a 2x3 grid, where the particle loss is not included but it is only an inset.
=#

gp = @pgf GroupPlot({
        group_style = {
            group_size = "4 by 2",
            vertical_sep = "0.0cm",
            horizontal_sep = "1.4cm",
            xlabels_at = "edge bottom",
            xticklabels_at = "edge bottom",
        },
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.02",
        xmin = whole10.final_times[9], xmax = whole10.final_times[75],
        width = "0.5\\textwidth",
        ticklabel_style = "/pgf/number format/fixed",
        xlabel = "\$t_f/t_R\$",
    },
    Axis({}, mn10),
    Axis({}, tn10),
    Axis({}, pl10),
    Axis({}, eff10),
    Axis({}, mn100),
    Axis({}, tn100),
    Axis({}, pl100),
    Axis({}, eff100),
)
display("gfx/test.pdf", gp)
