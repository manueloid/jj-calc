using PGFPlotsX, DelimitedFiles, Colors 
##{{{Fidelity plot eSTA vs STA in anharmonic trap ##src
ho_adiab = readdlm("adiabatic_ho.dat")[:,2]
ho_sta = readdlm("adiabatic_sta.dat")[:,2]
time = readdlm("adiabatic_sta.dat")[:,1]
time = time / time[end]
adiab_data  = @pgf Table([ time, ho_adiab ])
sta_data = @pgf Table([ time, ho_sta ])
delta = .1
ax = @pgf Axis({xtick_pos = "left",ytick_pos = "left",
                ylabel = "F", xlabel = raw"$t_f$",
                xmin = 0-.05, xmax = 1 + .05 , ymin = 0.2, ymax = 1.04,
                ylabel_style = "{xshift = 3.2cm ,yshift = -1cm, rotate = -90}", 
                xlabel_style = "{xshift  = 2.8cm, yshift = .2cm}",
                legend_style = "{at = {( 1.0,0.5 ) anchor = south}, draw = none}"
               }
              )

@pgf push!(ax, Plot({color = sta_color, very_thick}, sta_data))
@pgf push!(ax, Plot({color = adiab_color, very_thick, dashed}, adiab_data))
@pgf push!(ax, Legend("STA", "Adiabatic"))
pgfsave("sta_ho.pdf", ax)
##}}} ##src
##{{{Shape of the control parameter for the STA protocol and the adiabatic scheme for the harmonic trap ##src
ax = @pgf Axis({xtick_pos = "left",ytick_pos = "left",
                ylabel = raw"$\omega^2(t)/\omega_0^2$", xlabel = raw"$t/t_f$",
                xmin = -.05, xmax = 1.02, 
                #ylabel_style = "{ rotate = -90}", 
                xlabel_style = "{xshift  = 2.8cm, yshift = .2cm}",
                legend_style = "{at = {( 1.0,0.9 ) anchor = south}, draw = none}"
               }
              )
function control_ω(initial_omega::Float64, final_omega::Float64,final_time::Float64)
    γ = √(initial_omega/ final_omega)
    b = t -> 6.0*(γ - 1.0)*t^5 - 15.0*(γ -1.0)*t^4 + 10.0*(γ -1.0)*t^3 +1 # Polynomial obtained from the paper
    db = t ->ForwardDiff.derivative(b,t)
    ddb = t -> ForwardDiff.derivative(db,t)
    ω2 = t -> -ddb(t)/( b(t)*final_time^2 ) + initial_omega^2/b(t)^4
    return ω2 
end
sta_color = Colors.JULIA_LOGO_COLORS.blue
adiab_color = Colors.JULIA_LOGO_COLORS.red
time = 0.0:.001:1.0
ω0 = 250.0*pi
ωf = 2.5*pi
ω2 = control_ω(ω0,ωf, .2)
sta_control = Table([time,ω2.(time) / ω0^2 ])
adiab_control = Table(time,range(ω0^2, ωf^2, length = length(time))/ω0^2)
@pgf push!(ax, Plot({color = sta_color, ultra_thick}, sta_control))
@pgf push!(ax, Plot({color = adiab_color, ultra_thick, dashed}, adiab_control))
@pgf push!(ax, Legend("STA", "adiabatic"))
pgfsave("sta_ho_parameters.pdf", ax)
##}}} ##src
##{{{Fidelity plot eSTA vs STA in anharmonic trap ##src
ho_adiab = readdlm("sta_fidelity_anharmonic.dat")[:,2]
ho_sta = readdlm("enhanced_fidelity_anharmonc.dat")[:,2]
time = readdlm("sta_fidelity_anharmonic.dat")[:,1]
adiab_data  = @pgf Table([ time, ho_adiab ])
sta_data = @pgf Table([ time, ho_sta ])
delta = .1
x_inf = time[1] - delta; x_sup = time[end] + delta
y_inf = ho_sta[1]; y_sup = ho_sta[end] 
ax = @pgf Axis({xtick_pos = "left",ytick_pos = "left",
                ylabel = "F", xlabel = raw"$t_f$",
                xmin = x_inf - .1, xmax = x_sup + .1, ymin = y_inf - .05, ymax = y_sup + .02,
                ylabel_style = "{xshift = 3.2cm ,yshift = -1cm, rotate = -90}", 
                xlabel_style = "{xshift  = 2.8cm, yshift = .2cm}",
                legend_style = "{at = {( 1.0,0.5 ) anchor = south}, draw = none}"
               }
              )

@pgf push!(ax, Plot({color = sta_color, very_thick}, sta_data))
@pgf push!(ax, Plot({color = adiab_color, very_thick, dashed}, adiab_data))
@pgf push!(ax, Legend("eSTA", "STA"))
@pgf push!(ax, HLine({dashed, color = "black"}, 1.0))
pgfsave("enhanced_ho.pdf", ax)
##}}} ##src
##{{{Fidelity plot eSTA vs STA in Josephson Junction ##src
using CSV, DataFrames
data = CSV.read("fidelities_full.dat", DataFrame, header = 1 )
reduced = data[!,[ :tf,:eSTA,:STA,:N ]]
sort!(reduced,[:N,:tf])
##
select_particle = n -> filter(row->row.N == n , reduced)
sta_data = N -> @pgf Table([select_particle(N).tf, select_particle(N).STA])
esta_data =N-> @pgf Table([select_particle(N).tf, select_particle(N).eSTA])
ax = @pgf Axis({xtick_pos = "left",ytick_pos = "left",
                ylabel = "F", xlabel = raw"$t_f$",
                xmin = 0.045, xmax =.5, ymax = 1.01,
                ylabel_style = "{xshift = 3.2cm ,yshift = -1cm, rotate = -90}", 
                xlabel_style = "{xshift  = 2.8cm, yshift = .2cm}",
                legend_style = "{at = {( 1.0,0.5 ) anchor = south}, draw = none}"
               }
              )

@pgf push!(ax, Plot({color = sta_color, very_thick}, esta_data(30)))
@pgf push!(ax, Plot({color = adiab_color, very_thick, dashed}, sta_data(30)))
@pgf push!(ax, Plot({color = sta_color, very_thick}, esta_data(10)))
@pgf push!(ax, Plot({color = adiab_color, very_thick, dashed}, sta_data(10)))
@pgf push!(ax, Legend("eSTA", "STA"))
@pgf push!(ax, HLine({dashed, color = "black"}, 1.0))
pgfsave("enhanced_JJ.pdf", ax)
##}}} ##src
##{{{Control parameters for the Josephson Junction eSTA vs STA ##src
include("../../definitions.jl")
ax = @pgf Axis({xtick_pos = "left",ytick_pos = "left",
                ylabel = raw"$\omega^2(t)/\omega_0^2$", xlabel = raw"$t/t_f$",
                #xmin = -.05, xmax = 1.02, 
                #ylabel_style = "{ rotate = -90}", 
                xlabel_style = "{xshift  = 2.8cm, yshift = .2cm}",
                legend_style = "{at = {( 0.3,0.8 ) anchor = south}, draw = none}"
               }
              )
using SpecialPolynomials
function corrected_ω(control_parameter, correction_vector::Array)
    ys = vcat([ 0.0, correction_vector, 0.0 ]...)
    xs = range(0.0,1.0, length = length(ys) )|>collect
    l = Lagrange(xs,ys)
    return  t-> control_parameter(t) - l(t)
end
select_particle = n -> filter(row->row.N == n , data)
new_data = select_particle(20)
time = 0.0:.001:1.0
row = 55 
lambdas = new_data[row, r"λ"]|>Vector
ω_test = control_ω(new_data[row, "Λ0"], new_data[row, "Λf"], new_data[row, "tf"])
ω_esta = corrected_ω(ω_test, lambdas)
sta_control = @pgf Table([time, ω_test.(time) /ω_test(0)  ])
esta_control = @pgf Table([time, ω_esta.(time) /ω_test(0)])
@pgf push!(ax, Plot({color = sta_color, thick}, sta_control))
@pgf push!(ax, Plot({color = adiab_color, thick, dashed}, esta_control))
@pgf push!(ax, Legend("eSTA", "STA"))
row = 12
lambdas = new_data[row, r"λ"]|>Vector
ω_test = control_ω(new_data[row, "Λ0"], new_data[row, "Λf"], new_data[row, "tf"])
ω_esta = corrected_ω(ω_test, lambdas)
sta_control = @pgf Table([time, ω_test.(time)  /ω_test(0) ])
esta_control = @pgf Table([time, ω_esta.(time) /ω_test(0)])
@pgf push!(ax, Plot({color = sta_color, thick }, sta_control))
@pgf push!(ax, Plot({color = adiab_color, thick, dashed}, esta_control))
@pgf push!(ax, Legend("eSTA", "STA"))
pgfsave("enhanced_JJ_control.pdf", ax)
##}}} ##src
##{{{Josephson Junction example  ##src
x = -1:0.001:1
a = 1.1
b =1.3
pot = x-> -1/4*a*x^2 + 1/2b^2*x^4
ax = @pgf Axis({ticks = "none"})
ax = @pgf Axis({},
              )
xmin1 = findmin(pot.(x))[2]|>s->x[s]
xmin2 = -xmin1
point1 = @pgf Plot( {mark = "*",
                     color = "blue",
                     mark_options = "{scale = 4}",
                    },
                   Table([ [xmin1], [pot(xmin1) + .02] ]),
                  )
point2 = @pgf Plot( {mark = "*", color = "blue",mark_options = "{scale = 4}"}, Table([ [xmin2], [pot(xmin1) + .02] ]) )
point3 = @pgf Plot( {mark = "*",color = "red", mark_options = "{scale = 4}"}, Table([ [xmin1], [pot(xmin1) + .06] ]) )
@pgf push!(ax,Plot({color = "red", ultra_thick},Table([x, pot.(x)])))
@pgf push!(ax,point1, point2,point3)
pgfsave("jj_sketch.pdf", ax)
##}}} ##src
