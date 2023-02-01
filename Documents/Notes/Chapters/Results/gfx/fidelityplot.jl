# Plot styles and features

# Array where the colors are defined and stored 
using Colors, PGFPlotsX
colors = (
    red=colorant"#FF0000",
    yellow=colorant"#FFa000",
    blue=colorant"#0000FF",
    green=colorant"#00FF00",
    black=colorant"#000000"
)

# Array where all the line styles are stored
styles = (
    ldash="loosely dashed",
    dash=" dashed",
    solid="solid",
    dot="dash pattern={on 1pt off 1pt}",
    dot_dash="dash pattern={on 4pt off 1pt on 1pt off 1pt}"
)
opt = [@pgf {style = styles[index], color = colors[index], mark = "", thick} for index in 1:5]
# Reading from file
using DataFrames, CSV
filename(np::Int64, feature="fidelity") = "./Documents/Notes/Chapters/Results/gfx/data/$(feature)_np$np.dat"
data(np::Int64, feature="fidelity") = CSV.read(filename(np, feature), DataFrame; header=true)
df = data(30, "robustness")
sort!(df, [:tf])
pgfplot(data::DataFrame, options, type::Symbol) = @pgf Plot(options, Table(data.tf, data[!, type]))
## Creating plots
p1 = pgfplot(df, opt[2], :sta)
p2 = @pgf PlotInc(opt[2], Table(df.tf, df.full))
for opts in zip(styles, colors)

end

# Setting up axes and actual plottin
ax = Axis([p1, p2])
pgfsave("./Documents/Notes/Chapters/Results/gfx/test.pdf", ax)
