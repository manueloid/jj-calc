include("./src/ConstantQuantities.jl")
function fidelity_all(np::Int64, tspan::Vector{Float64}, Λf::Float64=50.0; nlambda::Int64=5, ω0::Float64=2.0)
    f = ControlParameterFull(ω0, 2.0√(Λf + 1.0), 0.2, np)
    i = ControlParameterInt(ω0, 2.0√(Λf + 1.0), 0.2, np)
    whole = [
        fidelity_time(f, tspan; nlambda=nlambda), # Fidelity of eSTA with the full Hamiltonian and hessian
        fidelity_time(f, tspan; nlambda=nlambda, hessian=false), # Fidelity of eSTA with the full Hamilhonian and NO hessian
        fidelity_time(i, tspan; nlambda=nlambda), # Fidelity of eSTA with the intermediate Hamiltonian and hessian
        fidelity_time(i, tspan; nlambda=nlambda, hessian=false), # Fidelity of eSTA with the intermediate Hamilhonian and NO hessian
        fidelity_time(f, tspan; esta=false) # STA fidelity ( as esta is set to false)
    ]
    return whole
end

##
#Calculation and plotting
##
tspan = range(0.04, 1.0, length=100) |> collect
test = fidelity_all(50, tspan)

##
# Plotting stuff: I need style and colors, I will store them in two different arrays
# I decided to use named tuples as they should give me finer control on the plotting options
# but at the same type they give me the ability to loop through them if needed
##
using PGFPlotsX, Colors

# Array where the colors are defined and stored 
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

for index in axes(styles, 1)
    println(index)
end


ax = @pgf Axis()

for (index, style) in enumerate(styles)
    @pgf opt = {style = style, "thick", color = colors[index]}
    push!(ax, @pgf Plot(opt, Table(tspan, test[index])))
end
ax

@pgf opt = {style = styles.solid, "thick", color = colors.red}
push!(ax, @pgf Plot(opt, Table(tspan, test[1])))
@pgf opt = {style = styles.solid, "thick", color = colors.blue}
push!(ax, @pgf Plot(opt, Table(tspan, test[5])))
