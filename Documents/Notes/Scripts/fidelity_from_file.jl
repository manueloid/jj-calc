##{{{Fidelity plot eSTA vs STA in Josephson Junction ##src
using CSV, DataFrames, PGFPlotsX, Colors
data = CSV.read("Data/fidelities_full.dat", DataFrame, header = 1 )
reduced = data[!,[ :tf,:eSTA,:STA,:N ]]
sort!(reduced,[:N,:tf])
##
sta_color = Colors.JULIA_LOGO_COLORS.blue
adiab_color = Colors.JULIA_LOGO_COLORS.red
select_particle = n -> filter(row->row.N == n , reduced) # Select the rows relative to the number of particles 
sta_data = N -> @pgf Table([select_particle(N).tf, select_particle(N).STA]) # Returns a function to create the STA fidelity over time table for N particles 
esta_data = N-> @pgf Table([select_particle(N).tf, select_particle(N).eSTA])# Returns a function to create the eSTA fidelity over time table for N particles 
ax = @pgf Axis({xtick_pos = "left",ytick_pos = "left",
                ylabel = "F", xlabel = raw"$t_f$",
                xmin = 0.045, xmax =0.5, ymax = 1.01,
                ylabel_style = "{xshift = 3.2cm ,yshift = -1cm, rotate = -90}", 
                xlabel_style = "{xshift  = 2.8cm, yshift = .2cm}",
                legend_style = "{at = {( 1.0,0.5 ) anchor = south}, draw = none}"
               }
              )
N1, N2 = 10,30 # variables to set the total number of particles
@pgf push!(ax, Plot({color = sta_color, very_thick}, esta_data(N1))) 
@pgf push!(ax, Plot({color = sta_color, very_thick, dashed}, sta_data(N1)))
@pgf push!(ax, Plot({color = adiab_color, very_thick}, esta_data(N2)))
@pgf push!(ax, Plot({color = adiab_color, very_thick, dashed}, sta_data(N2)))
@pgf push!(ax, Legend("N = $N1", "","N = $N2"))
@pgf push!(ax, HLine({dashed, color = "black"}, 1.0))
pgfsave("../gfx/enhanced_JJ.pdf", ax)
##}}} ##src
