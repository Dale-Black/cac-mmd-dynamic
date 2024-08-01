### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ dbf189da-4eb1-11ef-00f4-eb12751279bc
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate(temp = true)
	Pkg.add(url = "https://github.com/kczimm/Attenuations.jl")
	Pkg.add("Unitful")
	Pkg.add("PlutoUI")
end

# ╔═╡ 1ebe8d80-705e-41e1-8854-52de02a96f42
using PlutoUI: TableOfContents

# ╔═╡ b28fcef2-1218-4171-8b18-af7d52d5c09b
using Attenuations

# ╔═╡ ecacc06c-21f3-46b8-9b36-d50ac9daf45f
using Attenuations: Material, Compound

# ╔═╡ 1ef09810-51c3-467b-a7c6-94a48ef348a8
using Unitful: g, cm, mL, mol, Quantity

# ╔═╡ 50271d9a-d0dc-4400-94f2-7dc172a3844e
TableOfContents()

# ╔═╡ cdf203de-9266-47e1-b932-3510f10012d0
md"""
# Main Functions
"""

# ╔═╡ 73b08e7c-8937-4a2c-9eb3-63b8f1b3c174
md"""
## `get_mass_fraction_dict`
"""

# ╔═╡ 42996dc4-bed7-469f-80bf-e88da0abe691
# function get_mass_fraction_dict(
# 	solute_elements,
# 	solvent_elements,
# 	concentration_solute,
# 	solute_of_interest)
	
# 	# Calculate molar masses of solute elements
# 	mm_solute_elements = Dict()
# 	for (element, _) in solute_elements
# 		mm_element = element.Z / element.ZAratio * g/mol
# 		mm_solute_elements[element.name] = mm_element
# 	end

# 	# Calculate molar masses of solvent elements
# 	mm_solvent_elements = Dict()
# 	for (element, _) in solvent_elements
# 		mm_element = element.Z / element.ZAratio * g/mol
# 		mm_solvent_elements[element.name] = mm_element
# 	end

# 	# Calculate total molar mass of solute
# 	mm_solute = 0g/mol
# 	for (element, count) in solute_elements
# 		mm_solute += mm_solute_elements[element.name] * count
# 	end

# 	# Calculate total molar mass of solvent
# 	mm_solvent = 0g/mol
# 	for (element, count) in solvent_elements
# 		mm_solvent += mm_solvent_elements[element.name] * count
# 	end

# 	# Find the molar mass of the solute of interest
# 	mm_solute_of_interest = mm_solute_elements[solute_of_interest]

# 	# Calculate the moles of the solute of interest
# 	moles_solute_of_interest = concentration_solute / mm_solute_of_interest * cm^3

# 	# Calculate mass of the solute compound
# 	mass_solute_compound = moles_solute_of_interest * mm_solute

# 	# Calculate masses of individual solute elements
# 	mass_solute_elements = Dict()
# 	for (element, count) in solute_elements
# 		mass_element = count * moles_solute_of_interest * mm_solute_elements[element.name]
# 		mass_solute_elements[element.name] = mass_element
# 	end

# 	# Calculate the mass of solvent needed
# 	mass_solvent = 1g - mass_solute_compound

# 	# Calculate moles of the solvent compound
# 	moles_solvent = mass_solvent / mm_solvent

# 	# Calculate masses of individual solvent elements
# 	mass_solvent_elements = Dict()
# 	for (element, count) in solvent_elements
# 		mass_element = count * moles_solvent * mm_solvent_elements[element.name]
# 		mass_solvent_elements[element.name] = mass_element
# 	end

# 	# Combine masses of solute and solvent elements
# 	combined_masses = Dict()
# 	for (symbol, mass) in mass_solute_elements
# 		combined_masses[symbol] = mass
# 	end
# 	for (symbol, mass) in mass_solvent_elements
# 		if haskey(combined_masses, symbol)
# 			combined_masses[symbol] += mass
# 		else
# 			combined_masses[symbol] = mass
# 		end
# 	end

# 	# Calculate mass fractions
# 	mass_fractions = Dict()
# 	total_mass_solution = 1g
# 	for (symbol, mass) in combined_masses
# 		mass_fractions[symbol] = mass / total_mass_solution
# 	end

# 	mass_fracs_updated = Dict{Int, Float64}(
# 	    Elements[Symbol(element)].Z => value
# 	    for (element, value) in mass_fractions
# 	)

# 	return mass_fracs_updated
# end

# ╔═╡ 5f26ca76-c776-4f58-85c2-abcb67f95aad
function get_mass_fraction_dict(
	solute_elements,
	solvent_elements,
	concentration_solute,
	solute_of_interest)
	
	# Calculate molar masses of solute elements
	mm_solute_elements = Dict()
	for (element, _) in solute_elements
		mm_element = element.Z / element.ZAratio * g/mol
		mm_solute_elements[element.name] = mm_element
	end

	# Calculate molar masses of solvent elements
	mm_solvent_elements = Dict()
	for (element, _) in solvent_elements
		mm_element = element.Z / element.ZAratio * g/mol
		mm_solvent_elements[element.name] = mm_element
	end

	# Calculate total molar mass of solute
	mm_solute = 0g/mol
	for (element, count) in solute_elements
		mm_solute += mm_solute_elements[element.name] * count
	end

	# Calculate total molar mass of solvent
	mm_solvent = 0g/mol
	for (element, count) in solvent_elements
		mm_solvent += mm_solvent_elements[element.name] * count
	end

	# Calculate the mass of the solute (assuming total volume is 1 cm^3)
	mass_solute = concentration_solute * cm^3

	# Calculate masses of individual solute elements
	mass_solute_elements = Dict()
	for (element, count) in solute_elements
		mass_fraction_element = count * mm_solute_elements[element.name] / mm_solute
		mass_solute_elements[element.name] = mass_fraction_element * mass_solute
	end

	# Calculate the mass of solvent needed (assuming the total mass of solution is 1g)
	mass_solvent = 1g - mass_solute

	# Calculate masses of individual solvent elements
	mass_solvent_elements = Dict()
	for (element, count) in solvent_elements
		mass_fraction_element = count * mm_solvent_elements[element.name] / mm_solvent
		mass_solvent_elements[element.name] = mass_fraction_element * mass_solvent
	end

	# Combine masses of solute and solvent elements
	combined_masses = Dict()
	for (symbol, mass) in mass_solute_elements
		combined_masses[symbol] = mass
	end
	for (symbol, mass) in mass_solvent_elements
		if haskey(combined_masses, symbol)
			combined_masses[symbol] += mass
		else
			combined_masses[symbol] = mass
		end
	end

	# Calculate mass fractions
	mass_fractions = Dict()
	total_mass_solution = 1g
	for (symbol, mass) in combined_masses
		mass_fractions[symbol] = mass / total_mass_solution
	end

	mass_fracs_updated = Dict{Int, Float64}(
	    Elements[Symbol(element)].Z => value
	    for (element, value) in mass_fractions
	)

	return mass_fracs_updated
end

# ╔═╡ 894cc8d7-7bfb-439e-b760-b73eab3b3b46
md"""
## `create_material`
"""

# ╔═╡ 876844a6-e0d0-4fd9-8953-f566371cee05
function create_material(name, mass_frac_dict; density = 1g/cm^3)
	return Material(
		name,
		NaN,
		NaN,
		density,
		mass_frac_dict
	)
end

# ╔═╡ 0ce91def-d568-4a21-b354-48411fa639f9
md"""
## `lac_to_hu`
"""

# ╔═╡ 223f8a56-5343-4670-b98d-ed1b19d44d63
function lac_to_hu(lac_material, energy)
	return 1000 * (lac_material - μ(Materials.water, energy)) / μ(Materials.water, energy)
end

# ╔═╡ d562648e-43d6-4a83-aa1f-e79a9e99c1c1
md"""
# Test
"""

# ╔═╡ d9a4f48c-7121-4ef2-b97b-69b4a152c947
md"""
!!! danger
	Currently there are a few errors:

	1. All of the high energy (150 keV) values are pretty far off compared to measured for calcium and iodine. Seems more like 130 keV for iodine.
	2. Anything above 300 mg/mL concentration (e.g. 400 mg/mL) causes `get_mass_fraction_dict` to error for calcium.
	3. Seems like CaCl2 in H2O is not the right mixture. Maybe hydroxyapatite or some other material.
"""

# ╔═╡ daf47b9b-a421-4906-a785-1f8d4472781c
begin
	const LOW_ENERGY = 70keV
	const HIGH_ENERGY = 130keV
end

# ╔═╡ 891c79ed-e27b-4bc6-968e-9bd0f1932ee0
md"""
## Calcium/Water Mixture
"""

# ╔═╡ 59f0f2eb-0162-4c1c-800a-11507bc15fd5
hydroxyapatite_solute = (
Elements[:Calcium] => 10,
Elements[:Phosphorus] => 6,
Elements[:Oxygen] => 26,
Elements[:Hydrogen] => 2
)

# ╔═╡ 6e79167f-bbd1-4d87-a3ca-542f1f18804b
water_solvent = (
Elements[:Hydrogen] => 2, 
Elements[:Oxygen] => 1
)

# ╔═╡ ddb41a39-52c1-45b4-817c-17b5eeb9b01c
mass_fractions_hydroxyapatite_rod = get_mass_fraction_dict(
	hydroxyapatite_solute,
	water_solvent,
	0.700g/cm^3,
	"Calcium"
)

# ╔═╡ 46d5c3fa-c8e9-49f8-bbf1-6544e10877c4
begin
	total_fraction = 0
	for (_, fraction) in mass_fractions_hydroxyapatite_rod
			total_fraction += fraction
	end
	total_fraction
end

# ╔═╡ 24c13049-8d24-4895-a744-8fd8c1599441
material_hydroxyapatite_rod = create_material(
	"Ca (as 2 * Ca5(PO4)3(OH)) in Water",
	mass_fractions_hydroxyapatite_rod
)

# ╔═╡ fadaeefb-0043-49b7-a85a-7d11e98be315
lac_to_hu(
	μ(material_hydroxyapatite_rod, LOW_ENERGY),
	LOW_ENERGY
)

# ╔═╡ deac5298-edec-4a2a-b428-496a3d939fb6
lac_to_hu(
	μ(material_hydroxyapatite_rod, HIGH_ENERGY),
	HIGH_ENERGY
)

# ╔═╡ 17155640-1c9f-4770-a3d7-3fd2590a6260
mass_fractions_calcium_rod = get_mass_fraction_dict(
	(Elements[:Calcium] => 1, Elements[:Chlorine] => 2), 
	(Elements[:Hydrogen] => 2, Elements[:Oxygen] => 1), 
	0.700g/cm^3,
	"Calcium"
)

# ╔═╡ 5eef21e1-9401-49d3-950c-00402313ce72
# material_calcium_rod = create_material(
# 	"Ca (as CaCl2) in Water",
# 	mass_fractions_calcium_rod
# )

# ╔═╡ 2995310d-5cc4-42c9-9b0a-b6ed1d223a23
# lac_to_hu(
# 	μ(material_calcium_rod, LOW_ENERGY),
# 	LOW_ENERGY
# )

# ╔═╡ d8c23b76-eda1-4ac5-8b46-d4b70326c8fa
# lac_to_hu(
# 	μ(material_calcium_rod, HIGH_ENERGY),
# 	HIGH_ENERGY
# )

# ╔═╡ f22a2b12-1b75-4491-91dd-125879c8a7a0
md"""
## Iodine/Water Mixture
"""

# ╔═╡ be33f88e-54b8-4217-af13-d62819796ea0
mass_fractions_iodine_rod = get_mass_fraction_dict(
	(Elements[:Sodium] => 1, Elements[:Iodine] => 1), 
	(Elements[:Hydrogen] => 2, Elements[:Oxygen] => 1),
	0.0025g/cm^3,
	"Iodine"
)

# ╔═╡ 5150f992-d88e-487d-815e-6d2b7b71ca9a
material_iodine_rod = create_material(
	"I (as NaI) in Water",
	mass_fractions_iodine_rod
)

# ╔═╡ f2d61494-8df4-4a7b-bde3-bd7e0cbdce4d
 lac_to_hu(
	 μ(material_iodine_rod, LOW_ENERGY),
	 LOW_ENERGY
 )

# ╔═╡ 66a901e2-09a0-4b52-aa1f-f828009841d8
lac_to_hu(
	μ(material_iodine_rod, HIGH_ENERGY),
	HIGH_ENERGY
)

# ╔═╡ Cell order:
# ╠═dbf189da-4eb1-11ef-00f4-eb12751279bc
# ╠═1ebe8d80-705e-41e1-8854-52de02a96f42
# ╠═b28fcef2-1218-4171-8b18-af7d52d5c09b
# ╠═ecacc06c-21f3-46b8-9b36-d50ac9daf45f
# ╠═1ef09810-51c3-467b-a7c6-94a48ef348a8
# ╠═50271d9a-d0dc-4400-94f2-7dc172a3844e
# ╟─cdf203de-9266-47e1-b932-3510f10012d0
# ╟─73b08e7c-8937-4a2c-9eb3-63b8f1b3c174
# ╠═42996dc4-bed7-469f-80bf-e88da0abe691
# ╠═5f26ca76-c776-4f58-85c2-abcb67f95aad
# ╟─894cc8d7-7bfb-439e-b760-b73eab3b3b46
# ╠═876844a6-e0d0-4fd9-8953-f566371cee05
# ╟─0ce91def-d568-4a21-b354-48411fa639f9
# ╠═223f8a56-5343-4670-b98d-ed1b19d44d63
# ╟─d562648e-43d6-4a83-aa1f-e79a9e99c1c1
# ╟─d9a4f48c-7121-4ef2-b97b-69b4a152c947
# ╠═daf47b9b-a421-4906-a785-1f8d4472781c
# ╟─891c79ed-e27b-4bc6-968e-9bd0f1932ee0
# ╠═59f0f2eb-0162-4c1c-800a-11507bc15fd5
# ╠═6e79167f-bbd1-4d87-a3ca-542f1f18804b
# ╠═ddb41a39-52c1-45b4-817c-17b5eeb9b01c
# ╠═46d5c3fa-c8e9-49f8-bbf1-6544e10877c4
# ╠═24c13049-8d24-4895-a744-8fd8c1599441
# ╠═fadaeefb-0043-49b7-a85a-7d11e98be315
# ╠═deac5298-edec-4a2a-b428-496a3d939fb6
# ╠═17155640-1c9f-4770-a3d7-3fd2590a6260
# ╠═5eef21e1-9401-49d3-950c-00402313ce72
# ╠═2995310d-5cc4-42c9-9b0a-b6ed1d223a23
# ╠═d8c23b76-eda1-4ac5-8b46-d4b70326c8fa
# ╟─f22a2b12-1b75-4491-91dd-125879c8a7a0
# ╠═be33f88e-54b8-4217-af13-d62819796ea0
# ╠═5150f992-d88e-487d-815e-6d2b7b71ca9a
# ╠═f2d61494-8df4-4a7b-bde3-bd7e0cbdce4d
# ╠═66a901e2-09a0-4b52-aa1f-f828009841d8
