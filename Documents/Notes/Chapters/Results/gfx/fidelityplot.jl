##
# Plot styles and features
##

# Array where the colors are defined and stored 
using Colors, PGFPlotsX
colors = (
    red=colorant"#FF0000",
    yellow=colorant"#FFFF00",
    blue=colorant"#0000FF",
    green=colorant"#00FF00",
    black=colorant"#000000"
)

# Array where all the line styles are stored
styles = (
    ldash="loosely dashed",
    dash=" dashed",
    solid="solid",
    dot="dotted",
    dot_dashed="dash pattern={on 4pt off 1pt on 1pt off 1pt}"
)

##
# Reading from file
##
using DataFrames, CSV
filename(np::Int64, feature = "fidelity") = "./Documents/Notes/Chapters/Results/gfx/data/$(feature)_np$np.dat"
data(np::Int64) = CSV.read(filename(np), DataFrame; header = true)
