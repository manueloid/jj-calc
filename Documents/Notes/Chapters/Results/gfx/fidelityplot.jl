gfplot(data::DataFrame, options, type::Symbol) = @pgf Plot(options, Table(data.tf, data[!, type]))
## Creating plots
protocols = propertynames(df)[2:end]
ax = Axis()
[ push!(ax, pgfplot(df, opt[index], protocols[index])) for index in 1:5 ]
# Setting up axes and actual plottin
pgfsave("./Documents/Notes/Chapters/Results/gfx/test.pdf", ax)
run(`zathura ./Documents/Notes/Chapters/Results/gfx/test.pdf`)
