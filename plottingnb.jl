### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a50603dc-e11c-11ed-0121-2f6b8cb66545
using CSV, DataFrames, PGFPlotsX, Colors, PlutoUI

# ╔═╡ 9fe7773c-91ff-4ead-9674-e01482daa8cf
html"""<style>
main {
    max-width: 1000px;
}
"""

# ╔═╡ 7068afff-532d-4e1f-9185-fdb00444abb8
md"""
### Select the number of particles $(@bind np Select([10, 50]))
"""

# ╔═╡ ea1a63dc-3c5c-41dc-847e-8ab412b2e163
md"""
First I will load the file and store the content into a DataFrame
"""

# ╔═╡ 7cb08aba-7ee3-42bc-97cc-01a92a483c88
md"""
Now I need to find a way to filter the data based on the number of particles
"""

# ╔═╡ ce978666-f44d-4603-beae-059f991ac520
filter_val(df::DataFrame, col_name::Symbol, value) = filter(row -> row[col_name]==value, df)

# ╔═╡ 7a6adb36-4a60-4a8f-9fe3-be95ee237d7b
filter_val(df::DataFrame, col_name::Symbol, min, max) = filter(row -> row[col_name] > min && row[col_name] < max, df)

# ╔═╡ 67cdf521-c5f6-4cc4-9f12-1c04a6a5992d
function filter_val(df::DataFrame, N::Int64, λ::Int64; Λ0::Float64=0.0, Λf::Float64=50.0)
	filter(row -> row.N == N && row.λ == λ , df)
end

# ╔═╡ 77bd3260-d0dc-43b0-aedc-1bdd1196b4bd
begin 
	df = CSV.read("./data/whole_data.dat", DataFrame)	
	df1 = filter_val(df, np, 1)
	df = filter_val(df, np, 5)
end

# ╔═╡ 2ae5ca0c-e2b6-4930-9626-0754b032a011
md"""
### Plotting Options
#### Colors
"""

# ╔═╡ 586f1eac-1671-4da6-a03a-86b020f51378
colors = (
    black=colorant"#000000", # STA	
    red=colorant"#FF0000", # eSTA Full Hamiltonian with Hessian 
    blue=colorant"#0000FF", # eSTA Intermediate Hamiltonian with Hessian
    yellow=colorant"#FFa000", # eSTA Full Hamiltonian with original version
    green=colorant"#00FF00"# eSTA Intermediate Hamiltonian original version
)

# ╔═╡ f02d195c-ea8e-4633-b0b6-cd1a1fde7d53
md"""
#### Line Style
"""

# ╔═╡ b43b91bf-81a8-49cb-b676-cc985ed456db
styles = (
    solid="solid",                         # eSTA Intermediate Hamiltonian with Hessian
    dot_dash="dash pattern={on 4pt off 1pt on 1pt off 1pt}",#STA
    ldash="dash pattern={on 2pt off 2pt}",# eSTA Full Hamiltonian with Hessian 
    dash=" dashed",                        # eSTA Full Hamiltonian with original version
    dot="dash pattern={on 1pt off 1pt}",  # eSTA Intermediate Hamiltonian original version
)

# ╔═╡ bfa45357-f1bc-4f47-804d-f07ff9a62eea
md"""
### Setting up the options for the lines
"""

# ╔═╡ 9e64b389-3b31-45ae-baa8-ac72ae0d0885
begin
	esta_opt = @pgf {color=colors.red, line_width=1, style=styles.solid}
	sta_opt = @pgf{color=colors.black, line_width=1, style=styles.dash}
	ad_opt = @pgf {color=colors.green, line_width=1, style=styles.dot_dash}
	extra_opt = @pgf {color=colors.yellow, line_width=1, style=styles.ldash}
end

# ╔═╡ 7520d9ee-2cf4-41a6-b5db-65dc27581512
md"""
### Here I will plot the fidelity
"""

# ╔═╡ 930454f1-c558-42db-9a6c-51c60359ae21
md"""
Min tf $(@bind Fmin NumberField(0.0:0.01:1.0, default=0.0))

Min tf $(@bind Fmax NumberField(0.0:0.01:1.0, default=1.0))
"""

# ╔═╡ 81cc7b20-11c1-4281-af19-428a987017d3
begin
	f_df = filter_val(df, :tf, Fmin, Fmax)
	f_df1 = filter_val(df1, :tf, Fmin, Fmax)
	@pgf TikzPicture({ "scale" => 3 }, 
	Axis({xlabel = raw"t/\tau", ylabel = "F"},
		Plot(esta_opt, Table(f_df.tf, f_df.F_eSTA)),
		Plot(sta_opt, Table(f_df.tf, f_df.F_STA)),
		Plot(ad_opt, Table(f_df.tf, f_df.F_ad)),
		Plot(extra_opt, Table(f_df1.tf, f_df1.F_eSTA))
)
	)
end

# ╔═╡ f220b4fc-9a45-4923-b303-6358079e142b
md"""
### Plotting the robustness with respect to the modulation error
"""

# ╔═╡ d75af0c8-5c53-4e91-b5a5-938df476ea80
md"""
Min $t_f$ $(@bind mnmin NumberField(0.0:0.01:1.0, default = 0.0))

Max $t_f$ $(@bind mnmax NumberField(0.0:0.01:1.0, default = 1.0))
"""

# ╔═╡ a96a2d64-1a2b-4aaf-a5c1-c86bd50f8f80
begin
	mn_df = filter_val(df, :tf, mnmin, mnmax)
	@pgf TikzPicture({ "scale" => 3 }, 
	Axis({xlabel=raw"t/\tau", ylabel=raw"S"},
		Plot(esta_opt, Table(mn_df.tf, mn_df.Mn_eSTA)),
		Plot(sta_opt, Table(mn_df.tf, mn_df.Mn_STA)),
)
)
end

# ╔═╡ a172948e-f0b5-4285-8a75-0dc1f58cbab2
md"""
### Plotting the time noise error
"""

# ╔═╡ 31a094cf-fe38-4691-8ea2-2f0feda83e37
md"""
Min tf $(@bind tnmin NumberField(0.0:0.01:1.0, default=0.0))

Min tf $(@bind tnmax NumberField(0.0:0.01:1.0, default=1.0))
"""

# ╔═╡ 9b947584-438d-4b3c-9b48-e7f50bb0fad3
begin
	tn_df = filter_val(df, :tf, tnmin, tnmax)
	@pgf TikzPicture({ "scale" => 3 }, 
	Axis({xlabel=raw"t/\tau", ylabel=raw"S"},
		Plot(esta_opt, Table(tn_df.tf, tn_df.Tn_eSTA)),
		Plot(sta_opt, Table(tn_df.tf, tn_df.Tn_STA)),
)
)
end

# ╔═╡ 180d9fc7-4bf4-450b-8d1b-d7b8e391e50c
md"""
### Plotting the error with respect to the white noise
"""

# ╔═╡ 6b76b553-0a5a-472a-8e69-a7e92faf11b3
md"""
Min $t_f$ $(@bind wnmin NumberField(0.0:0.01:1.0, default = 0.0))

Max $t_f$ $(@bind wnmax NumberField(0.0:0.01:1.0, default = 1.0))
"""

# ╔═╡ d21513f6-bd3b-41aa-a93c-744139c628ea
begin
	wn_df = filter_val(df, :tf, wnmin, wnmax)
	@pgf TikzPicture({ "scale" => 3 }, 
	Axis({xlabel=raw"t/\tau", ylabel=raw"S"},
		Plot(esta_opt, Table(wn_df.tf, wn_df.Wn_eSTA)),
		Plot(sta_opt, Table(wn_df.tf, wn_df.Wn_STA)),
)
)
end

# ╔═╡ bc97859e-f474-46b5-86ec-1567829b351d
md"""
### Plotting the squeezing 
"""

# ╔═╡ b2b8e791-c825-41c3-827d-c8b10b66f1d2
ξ_id = CSV.read("./data/squeezing_id.dat", DataFrame)[np,2]	 

# ╔═╡ 6e9fd844-ec54-4b54-ac96-db234241665c
md"""
Min $t_f$ $(@bind sqmin NumberField(0.0:0.01:1.0, default = 0.0))

Max $t_f$ $(@bind sqmax NumberField(0.0:0.01:1.0, default = 1.0))
"""

# ╔═╡ d642e1c5-04bb-4b73-bb3d-62967e948828
begin
sq_df = filter_val(df, :tf, sqmin, sqmax)
@pgf TikzPicture({ "scale" => 3 }, 
	Axis({xlabel=raw"t/\tau", ylabel=raw"S"},
		Plot(esta_opt, Table(sq_df.tf, sq_df.Sq_eSTA)),
		Plot(sta_opt, Table(sq_df.tf, sq_df.Sq_STA)),
		Plot(ad_opt, Table(sq_df.tf, ξ_id*ones(length(sq_df.tf))))
)
)
end

# ╔═╡ 156fb0c9-bea9-474c-bc03-31d60d4d6b06
md"""
## Now for the control parameters
"""

# ╔═╡ 7d76c66a-3929-46f4-b822-90ca6dab181c
evo_df = CSV.read("./data/evo_data.dat", DataFrame)

# ╔═╡ 1bd450ee-8f62-42f8-b7e7-96a53a014d65
md"""
### These are the control functions for $(evo_df.N |> unique) particles and final time $(evo_df.t[end])
"""

# ╔═╡ e939e438-7d4e-4a0b-8354-2d94451da076
@pgf TikzPicture({scale = "3"}, 
	Axis({xlabel = raw"$t$", ylabel = raw"$\Lambda$"},
			Plot(esta_opt, Table(evo_df.t, evo_df.Λ_eSTA)),
			Plot(sta_opt, Table(evo_df.t, evo_df.Λ_STA)),
			Plot(extra_opt, Table(evo_df.t, evo_df.Λ_ad))
	)
)

# ╔═╡ 991b344d-9d3b-44fc-878c-caf6bff8cce2
md"""
### And for the same parameters, we have the evolution of the squeezing
"""

# ╔═╡ 476d287a-a968-4400-b167-99bebf93f811
@pgf TikzPicture({scale = 3},
	Axis({xlabel = "t", ylabel = raw"$\xi^2_s$"},
		Plot(esta_opt, Table(evo_df.t, evo_df.ξs_eSTA)),
		Plot(sta_opt, Table(evo_df.t, evo_df.ξs_STA)),
	)
)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
PGFPlotsX = "8314cec4-20b6-5062-9cdb-752b83310925"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
CSV = "~0.10.9"
Colors = "~0.12.10"
DataFrames = "~1.5.0"
PGFPlotsX = "~1.5.3"
PlutoUI = "~0.7.50"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "d923c0a495166b9155c552df7a6b86c2c9b3c8c1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefaultApplication]]
deps = ["InteractiveUtils"]
git-tree-sha1 = "c0dfa5a35710a193d83f03124356eef3386688fc"
uuid = "3f0dd361-4fe0-5fc6-8523-80b14ec94d85"
version = "1.1.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PGFPlotsX]]
deps = ["ArgCheck", "Dates", "DefaultApplication", "DocStringExtensions", "MacroTools", "OrderedCollections", "Parameters", "Requires", "Tables"]
git-tree-sha1 = "e98a6743775e107062be357560977c06850a79be"
uuid = "8314cec4-20b6-5062-9cdb-752b83310925"
version = "1.5.3"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "548793c7859e28ef026dba514752275ee871169f"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "0b829474fed270a4b0ab07117dce9b9a2fa7581a"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.12"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─9fe7773c-91ff-4ead-9674-e01482daa8cf
# ╠═a50603dc-e11c-11ed-0121-2f6b8cb66545
# ╟─7068afff-532d-4e1f-9185-fdb00444abb8
# ╟─ea1a63dc-3c5c-41dc-847e-8ab412b2e163
# ╠═77bd3260-d0dc-43b0-aedc-1bdd1196b4bd
# ╟─7cb08aba-7ee3-42bc-97cc-01a92a483c88
# ╠═ce978666-f44d-4603-beae-059f991ac520
# ╠═7a6adb36-4a60-4a8f-9fe3-be95ee237d7b
# ╠═67cdf521-c5f6-4cc4-9f12-1c04a6a5992d
# ╟─2ae5ca0c-e2b6-4930-9626-0754b032a011
# ╟─586f1eac-1671-4da6-a03a-86b020f51378
# ╟─f02d195c-ea8e-4633-b0b6-cd1a1fde7d53
# ╟─b43b91bf-81a8-49cb-b676-cc985ed456db
# ╟─bfa45357-f1bc-4f47-804d-f07ff9a62eea
# ╟─9e64b389-3b31-45ae-baa8-ac72ae0d0885
# ╟─7520d9ee-2cf4-41a6-b5db-65dc27581512
# ╠═81cc7b20-11c1-4281-af19-428a987017d3
# ╟─930454f1-c558-42db-9a6c-51c60359ae21
# ╟─f220b4fc-9a45-4923-b303-6358079e142b
# ╟─a96a2d64-1a2b-4aaf-a5c1-c86bd50f8f80
# ╟─d75af0c8-5c53-4e91-b5a5-938df476ea80
# ╟─a172948e-f0b5-4285-8a75-0dc1f58cbab2
# ╟─9b947584-438d-4b3c-9b48-e7f50bb0fad3
# ╟─31a094cf-fe38-4691-8ea2-2f0feda83e37
# ╟─180d9fc7-4bf4-450b-8d1b-d7b8e391e50c
# ╟─d21513f6-bd3b-41aa-a93c-744139c628ea
# ╟─6b76b553-0a5a-472a-8e69-a7e92faf11b3
# ╟─bc97859e-f474-46b5-86ec-1567829b351d
# ╟─b2b8e791-c825-41c3-827d-c8b10b66f1d2
# ╟─d642e1c5-04bb-4b73-bb3d-62967e948828
# ╟─6e9fd844-ec54-4b54-ac96-db234241665c
# ╟─156fb0c9-bea9-474c-bc03-31d60d4d6b06
# ╠═7d76c66a-3929-46f4-b822-90ca6dab181c
# ╟─1bd450ee-8f62-42f8-b7e7-96a53a014d65
# ╟─e939e438-7d4e-4a0b-8354-2d94451da076
# ╟─991b344d-9d3b-44fc-878c-caf6bff8cce2
# ╟─476d287a-a968-4400-b167-99bebf93f811
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
