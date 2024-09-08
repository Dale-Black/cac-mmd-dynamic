### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 3934542d-9041-430f-999e-5fefa9a92791
# ╠═╡ show_logs = false
begin
	using Pkg;
	Pkg.activate(temp = true)
	Pkg.add(url = "https://github.com/Dale-Black/Attenuations.jl")
	Pkg.add("PlutoUI")
	Pkg.add("CairoMakie")
end

# ╔═╡ 4ed0fa55-662f-485a-95b4-b4bcb74d38d8
using Attenuations: Material, eV, keV, g, cm, μ, μᵨ, Elements, Materials, ustrip

# ╔═╡ a1e078bb-3420-4882-878c-3e5107ef5048
using PlutoUI: TableOfContents

# ╔═╡ 8d8295fc-79a9-4b48-8dba-29b7a9d6ad8a
using CairoMakie

# ╔═╡ 7f52e2ee-4765-4c81-ba14-cd148dcd3da6
TableOfContents()

# ╔═╡ 1d0169dd-7f4c-4854-badd-4db61e9ca31d
md"""
## Create omnipaque
"""

# ╔═╡ 2adabd6a-6d4c-11ef-2dfe-196231b9b8cd
begin
	# Number of mg of iohexol in omnipaque 300
	const mg_iohexol_in_omnipaque = 647 # mg

	# Number of mg of saline in omnipaque 300
	const mg_saline_in_omnipaque = 1000 - mg_iohexol_in_omnipaque # 1000 mg - mg_iohexol_in_omnipaque mg

	# Concentration of saline within omnipaque 300
	const concentration_saline_in_omnipaque = mg_saline_in_omnipaque / 1000
end

# ╔═╡ 7a2c80cf-8a60-43b4-b934-09e9fb143ec7
concentration_H2O_in_omnipaque = 0.9 * concentration_saline_in_omnipaque

# ╔═╡ 43f19c5f-dae8-422a-94d0-8f133c7235cf
concentration_NaCl_in_omnipaque = 0.1 * concentration_saline_in_omnipaque

# ╔═╡ 42f91628-87a2-4143-bae7-b5cc04a20071
concentration_iohexol_in_omnipaque = mg_iohexol_in_omnipaque / 1000

# ╔═╡ 356a1505-b072-4b20-900d-372cd7ea4d8d
md"""
**Plug this into [NIST XCOM](https://physics.nist.gov/cgi-bin/Xcom/xcom2?Method=Mix&Output2=Hand) and create material based on values**

H20 $(concentration_H2O_in_omnipaque) \
NaCl $(concentration_NaCl_in_omnipaque) \
C19H26I3N3O9 $(concentration_iohexol_in_omnipaque)

**Which should return something like:**
```
Constituents (Atomic Number : Fraction by Weight)
 Z=1	: 0.338349
 Z=6	: 0.179812
 Z=7	: 0.033109
 Z=8	: 0.113457
 Z=11	: 0.013886
 Z=17	: 0.021414
 Z=53	: 0.299974
```

**And can then be used to create the `Material`**
"""

# ╔═╡ 8b38d417-5b8f-48f9-99a0-c42875f40934
omnipaque_material = Material(
	"Omnipaque 300",
	0.0, # unnecessary for LAC
	0.0eV, # unnecessary for LAC
	
	#= density => specific gravity 
	(found from https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=442aed6e-6242-4a96-90aa-d988b62d55e8)
	=#
	
	1.0g/cm^3, # not accurate according to NIH specs ^^, but fits Mendonca better
	Dict(
		1 => 0.338349,
		6 => 0.179812,
		7 => 0.033109,
		8 => 0.113457,
		11 => 0.013886,
		17 => 0.021414,
		53 => 0.299974,
	)
)

# ╔═╡ 9dcfc538-10be-4db4-b6db-270d19f22a69
md"""
## Check against Mendoca et. al (2014) values

* bone ~ [0.49, 0.29] ✅
* blood ~ [0.20, 0.16] ✅
* fat ~ [0.18, 0.14] ✅
* air ~ [0.00, 0.00] ✅


* omnipaque ~ [1.70, 0.40] ✅


* omnipaque and blood ~ [0.41, 0.19]
* bone and blood ~ [0.28, 0.19]
"""

# ╔═╡ 53352af8-6c32-4f31-98d9-33bbdac4711b
energies = [70keV, 140keV]

# ╔═╡ 2d88de39-220c-4836-baed-895a98df183c
md"""
### Simple
"""

# ╔═╡ 92300301-a2c3-4ed9-a152-d342fd476562
bone_lacs = μ(Materials.corticalbone, energies) # matches mendoca paper values

# ╔═╡ 2df9bf1b-766e-4ea6-aa00-effa865b0bc6
blood_lacs = μ(Materials.wholeblood, energies) # matches mendoca paper values

# ╔═╡ 44df6f47-8179-4858-b59d-e79bb541d723
fat_lacs = μ(Materials.adipose, energies) # matches mendoca paper values

# ╔═╡ 459e195d-abda-427c-807a-a6335af0c744
air_lacs = μ(Materials.air, energies) # matches mendoca paper values

# ╔═╡ c8d48bcd-eda1-4122-8452-1e64ace83835
md"""
### Omnipaque
"""

# ╔═╡ 0da7aa66-4d91-40f1-abb5-4f0dbd482709
omnipaque_lacs = μ(omnipaque_material, energies)

# ╔═╡ 212e0fea-87c1-4fd0-92dc-79bb4f18f38c
md"""
### Mixtures
"""

# ╔═╡ 9dbeda01-95e6-434c-916e-e9d390e08dfa
begin
	alpha = 0.14
	omnipaque_and_blood_lacs = (alpha * omnipaque_lacs) + ((1 - alpha) * blood_lacs)
end

# ╔═╡ 72844bd1-183b-4823-804e-007a09054ee2
begin
	beta = 0.28
	bone_and_blood_lacs = (beta * bone_lacs) + ((1 - beta) * blood_lacs)
end

# ╔═╡ 248c9fc2-f2be-4487-b012-40ad79ba6f3a
md"""
## Plot triangles
"""

# ╔═╡ bb6c43d6-32cb-4246-8cba-d99d262bf18d
let
	f = Figure(size = (800, 400))
	
	ax = Axis(
		f[1, 1],
		title = """
			Linear attenuation coeffiecients (LAC)
			of multiple materials at nominal densities
		""",
		xlabel = "LAC (70 keV)",
		ylabel = "LAC (140 keV)",
		yticks = 0.0:0.05:0.5,
		xticks = 0.0:0.2:1.8
	)

	points = [
		Point2f.(ustrip.(air_lacs)[1], ustrip.(air_lacs)[2]),
		Point2f.(ustrip.(omnipaque_lacs)[1], ustrip.(omnipaque_lacs)[2]),
		Point2f.(ustrip.(bone_lacs)[1], ustrip.(bone_lacs)[2]),
		Point2f.(ustrip.(blood_lacs)[1], ustrip.(blood_lacs)[2]),
		Point2f.(ustrip.(fat_lacs)[1], ustrip.(fat_lacs)[2]),
		Point2f.(ustrip.(bone_and_blood_lacs)[1], ustrip.(bone_and_blood_lacs)[2]),
		Point2f.(
			ustrip.(omnipaque_and_blood_lacs)[1],
			ustrip.(omnipaque_and_blood_lacs)[2]
		),
	]
			

	triplot!(
		points; 
		show_points = true,
		triangle_color = :lightblue,
		show_convex_hull = true
	)

	markersize = 20
	
	scatter!(
		ustrip.(air_lacs)[1],
		ustrip.(air_lacs)[2],
		label = "Air",
		markersize = markersize
	)
	scatter!(
		ustrip.(bone_lacs)[1],
		ustrip.(bone_lacs)[2],
		label = "Bone",
		markersize = markersize
	)
	scatter!(
		ustrip.(omnipaque_lacs)[1],
		ustrip.(omnipaque_lacs)[2],
		label = "Omnipaque 300",
		markersize = markersize
	)
	scatter!(
		ustrip.(blood_lacs)[1],
		ustrip.(blood_lacs)[2],
		label = "Blood",
		markersize = markersize
	)
	scatter!(
		ustrip.(fat_lacs)[1],
		ustrip.(fat_lacs)[2],
		label = "Fat",
		markersize = markersize
	)
	scatter!(
		ustrip.(bone_and_blood_lacs)[1],
		ustrip.(bone_and_blood_lacs)[2],
		label = "Bone and blood",
		markersize = markersize
	)
	scatter!(
		ustrip.(omnipaque_and_blood_lacs)[1],
		ustrip.(omnipaque_and_blood_lacs)[2],
		label = "Omnipaque 300 and blood",
		markersize = markersize
	)

	axislegend(ax; position = :rb)

	# axislegend(ax; position = :rb)
	xlims!(low = -0.01, high = 1.8)
	ylims!(low = -0.01, high = 0.5)

	f
end	

# ╔═╡ Cell order:
# ╠═3934542d-9041-430f-999e-5fefa9a92791
# ╠═4ed0fa55-662f-485a-95b4-b4bcb74d38d8
# ╠═a1e078bb-3420-4882-878c-3e5107ef5048
# ╠═7f52e2ee-4765-4c81-ba14-cd148dcd3da6
# ╟─1d0169dd-7f4c-4854-badd-4db61e9ca31d
# ╠═2adabd6a-6d4c-11ef-2dfe-196231b9b8cd
# ╠═7a2c80cf-8a60-43b4-b934-09e9fb143ec7
# ╠═43f19c5f-dae8-422a-94d0-8f133c7235cf
# ╠═42f91628-87a2-4143-bae7-b5cc04a20071
# ╟─356a1505-b072-4b20-900d-372cd7ea4d8d
# ╠═8b38d417-5b8f-48f9-99a0-c42875f40934
# ╟─9dcfc538-10be-4db4-b6db-270d19f22a69
# ╠═53352af8-6c32-4f31-98d9-33bbdac4711b
# ╟─2d88de39-220c-4836-baed-895a98df183c
# ╠═92300301-a2c3-4ed9-a152-d342fd476562
# ╠═2df9bf1b-766e-4ea6-aa00-effa865b0bc6
# ╠═44df6f47-8179-4858-b59d-e79bb541d723
# ╠═459e195d-abda-427c-807a-a6335af0c744
# ╟─c8d48bcd-eda1-4122-8452-1e64ace83835
# ╠═0da7aa66-4d91-40f1-abb5-4f0dbd482709
# ╟─212e0fea-87c1-4fd0-92dc-79bb4f18f38c
# ╠═9dbeda01-95e6-434c-916e-e9d390e08dfa
# ╠═72844bd1-183b-4823-804e-007a09054ee2
# ╟─248c9fc2-f2be-4487-b012-40ad79ba6f3a
# ╠═8d8295fc-79a9-4b48-8dba-29b7a9d6ad8a
# ╟─bb6c43d6-32cb-4246-8cba-d99d262bf18d
