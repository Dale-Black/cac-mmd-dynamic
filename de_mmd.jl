### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ ea478b7f-96d2-4802-9468-33f4bdf7e9d7
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate(temp = true)

	## Handle Python Stuff
	Pkg.add("CondaPkg")
	using CondaPkg
	CondaPkg.add("SimpleITK")
	CondaPkg.add("numpy")
	Pkg.add("PythonCall")
	
	Pkg.add("Statistics")
	Pkg.add("PlutoUI")
	Pkg.add("CairoMakie")
	Pkg.add("Random")
	Pkg.add("OrderedCollections")
	Pkg.add("LinearAlgebra")
	Pkg.add("Unitful")
	Pkg.add(url = "https://github.com/kczimm/Attenuations.jl")
	Pkg.add(url = "https://github.com/JuliaHealth/DICOM.jl")
end

# ╔═╡ 21d854cf-06a3-4214-a689-f57be77e84be
using Statistics: mean

# ╔═╡ 45e73788-4b7c-11ef-30bd-e161fe50a735
using DICOM: @tag_str, dcm_parse, dcmdir_parse, rescale!

# ╔═╡ 2a149318-4c73-4a57-a819-79f70adc1700
using PlutoUI: TableOfContents

# ╔═╡ e8af4ce9-72c7-4017-bd92-6fd3a419e973
using CairoMakie: Figure, Axis, lines!, scatter!, heatmap!

# ╔═╡ bed00701-b564-4a26-90eb-0f9132da0ce3
using Random: MersenneTwister, shuffle

# ╔═╡ d5dfce22-566f-4361-8cfc-237700840676
using OrderedCollections: OrderedDict

# ╔═╡ 7ec95459-8e70-4a06-9466-262b97fb871f
using LinearAlgebra: det

# ╔═╡ 1a04d22f-047b-42e7-84a3-8bd0fd222743
using Unitful: mg, mL, g, cm

# ╔═╡ d01cf5fa-a788-4efb-ad47-24e2a014f2f8
using Attenuations

# ╔═╡ 42b57d75-3b39-4c61-995d-8d0b20a6a7df
# ╠═╡ show_logs = false
using PythonCall

# ╔═╡ dd80ed67-d7c9-4551-a35f-c068f095d8a9
md"""
# Setup
"""

# ╔═╡ c9179fde-37af-4f7d-a87b-2301ea96b643
md"""
## Imports
"""

# ╔═╡ a88d38c6-c85a-414d-911d-146e20c7fd29
sitk = pyimport("SimpleITK")

# ╔═╡ 21415643-b7b9-420b-b7ef-a7c28f40e7e1
np = pyimport("numpy")

# ╔═╡ dab48e46-c732-44ff-9bed-391f94e4ca43
TableOfContents()

# ╔═╡ b06d9b63-67b3-413b-96cf-9035678c721c
md"""
## File Paths
"""

# ╔═╡ 07faae20-bd30-4412-bd57-27c95c8c8740
path_80kv = joinpath(pwd(), "data/80kV 40.dcm")

# ╔═╡ 3a6b40c4-a25a-4e28-9375-f06fe15f6cfc
path_135kv = joinpath(pwd(), "data/135kV 40.dcm")

# ╔═╡ a30419d5-2ca0-4af8-bc2b-676dcda84f60
path_70_keV = joinpath(pwd(), "data/70keV140mA 2.dcm")

# ╔═╡ ab512723-b5e5-496c-99c6-2277f70bb0a0
path_150_keV = joinpath(pwd(), "data/150keV.dcm")

# ╔═╡ 179fb57b-10c6-44c0-a47e-4d0d96f9f71b
md"""
## Dictionaries
"""

# ╔═╡ 9830994f-2f85-48fd-aae2-b70330019658
md"""
### Polyenergetic CT data
"""

# ╔═╡ 72c7e682-58db-4577-be77-10d5376604e8
begin
	Ca80550 = OrderedDict(
	    50 => (223, 156.55),
	    100 => (416, 189.29),
	    200 => (733, 246.30),
	    300 => (1096, 267.10),
	    400 => (1464, 282.52)
	)
	
	I80550 = OrderedDict(
	    2.0 => (52, 205.36),
	    2.5 => (99, 239.72),
	    5.0 => (163, 202.84),
	    7.5 => (289, 209.85),
	    10 => (357, 242.31),
	    15 => (533, 258.40),
	    20 => (712, 287.84)
	)
	
	BkgWater80550 = [-4, 181]
	
	Ca135200 = OrderedDict(
	    50 => (180, 64.80),
	    100 => (312, 73.64),
	    200 => (529, 76.86),
	    300 => (734, 74.64),
	    400 => (960, 82.30)
	)
	
	I135200 = OrderedDict(
	    2.0 => (33, 73.02),
	    2.5 => (47, 80.33),
	    5.0 => (92, 75.01),
	    7.5 => (118, 87.49),
	    10 => (180, 79.15),
	    15 => (269, 82.06),
	    20 => (329, 94.31)
	)
	
	BkgWater135200 = [-19, 61]
end

# ╔═╡ 0f7faad8-73cb-4417-8048-fca690c3e2db
md"""
### Monoenergetic (PCCT) data
"""

# ╔═╡ fd67cdc9-5e6f-4771-8d87-6a2183f3f93d
begin
	Ca70200 = OrderedDict(
	    50 => (189, 9.385),
	    100 => (329, 9.726),
	    200 => (605, 9.387),
	    300 => (867, 10.747),
	    400 => (1123, 11.241)
	)
	
	I70200 = OrderedDict(
	    2.0 => (40.350, 8.855),
	    2.5 => (54.613, 8.687),
	    5.0 => (117.650, 10.046),
	    7.5 => (180.218, 9.628),
	    10 => (242.005, 10.467),
	    15 => (352.105, 9.876),
	    20 => (476.491, 11.418)
	)
	
	BkgWater70200 = [0.077, 8.235]
	
	Ca150140 = OrderedDict(
	    50 => (153, 7.200),
	    100 => (226, 7.686),
	    200 => (375, 7.248),
	    300 => (523, 8.214),
	    400 => (676, 8.310)
	)
	
	I150140 = OrderedDict(
	    2.0 => (9.650, 8.062),
	    2.5 => (12.820, 7.345),
	    5.0 => (24.516, 9.425),
	    7.5 => (38.931, 8.479),
	    10 => (54.234, 7.859),
	    15 => (81.943, 9.025),
	    20 => (106.228, 8.688)
	)
	
	BkgWater150140 = [-0.431, 7.555]
end

# ╔═╡ 0c8c12fd-c439-44c7-a1e1-7f5f0674f2d3
md"""
# Function Declarations
"""

# ╔═╡ 935739af-dcde-49e4-bdc0-0875c7335587
md"""
## `DICOM.rescale!`
"""

# ╔═╡ 84ea30d1-f3b3-4a37-ae46-b991bd71072b
md"""
**Check that these values are correct**:

(verified in imageJ)

70 keV
- mean: -310.391
- min: -8192
- max: 1747

150 keV
- mean: -298.221
- min: -8192.0
- max: 723.0
"""

# ╔═╡ ba5ef8fe-4f03-4f01-8acf-27c990bed4d1
begin
	low_energy_mono = dcm_parse(path_70_keV)
	rescale!(low_energy_mono)
	low_energy_mono = low_energy_mono.PixelData
end;

# ╔═╡ b484d333-a9b8-4a8f-b325-6a360675f271
begin
	high_energy_mono = dcm_parse(path_150_keV)
	rescale!(high_energy_mono)
	high_energy_mono = high_energy_mono.PixelData
end;

# ╔═╡ b57e7f8e-0cb5-4591-bee4-6b755f8e0134
mean(low_energy_mono), minimum(low_energy_mono), maximum(low_energy_mono)

# ╔═╡ 3ac3b86e-20a3-48cb-bf17-8593353e9210
mean(high_energy_mono), minimum(high_energy_mono), maximum(high_energy_mono)

# ╔═╡ 2ba6333f-850b-47b7-b39c-60c2da627cb5
md"""
## `interpolate_edge`
"""

# ╔═╡ d59e48b9-7e94-4598-a590-2ead23389c81
function interpolate_edge(vertex1, vertex2; num_points=100)
    x1, y1 = vertex1
    x2, y2 = vertex2
    t = range(0, 1, length=num_points+1)
    x = x1 .+ t .* (x2 - x1)
    y = y1 .+ t .* (y2 - y1)
    return x, y
end

# ╔═╡ bc56c4ec-c232-4789-82d3-c290bbad77e1
let
	vertices = ([0,0], [1,0.5], [0.5,1.5])
	line1 = interpolate_edge(vertices[1], vertices[2])
	line2 = interpolate_edge(vertices[2], vertices[3])
	line3 = interpolate_edge(vertices[1], vertices[3])

	f = Figure()
		ax = Axis(f[1, 1])
		lines!(ax, line1[1], line1[2])
		lines!(ax, line2[1], line2[2])
		lines!(ax, line3[1], line3[2])
	f
end

# ╔═╡ fcb9cdfe-e516-4f34-afc9-bb5f66cf6b08
md"""
## `nearest_triangle_point`
"""

# ╔═╡ cba2ace6-90ff-4476-b576-0f87c0da0614
"""
	nearest_triangle_point(t_vertices, point)

Calculate the Hausdorff distance between a triangle and a point.

# Arguments
- `t_vertices`: List of 3 vertices defining the triangle
- `point`: The point to compare against the triangle

# Returns
- distance, projection_point
"""
function nearest_triangle_point(t_vertices, point)
	function find_closest_vertex(x1, y1, x2, y2, ptx, pty)
		dist1 = sqrt((x1 - ptx)^2 + (y1 - pty)^2)
		dist2 = sqrt((x2 - ptx)^2 + (y2 - pty)^2)
		return dist1 < dist2 ? (dist1, [x1, y1]) : (dist2, [x2, y2])
	end

	tol = 1e-9  # tolerance for floating-point comparisons
	dist = Inf
	proj_point = nothing
	ptx, pty = point

	# Loop through the three triangle edges to find min distance to the outside point
	# This can result in (a) projection point is on the triangle edge
	# or (b) projection point is on the extension line from the endpoints of the edge
	for i in 1:length(t_vertices)
		for j in (i+1):length(t_vertices)
			x1, y1 = t_vertices[i]
			x2, y2 = t_vertices[j]
			
			# Check if the triangle edge is vertical or horizontal
			if abs(x1 - x2) < tol  # edge is vertical
				x0 = x1
				y0 = pty
				# Check if the projection point [x0, y0] is on the triangle edge
				if min(y1, y2) < y0 < max(y1, y2)
					dist_from_line = abs(x0 - ptx)
					proj_point_from_line = [x0, y0]
				else
					# The projection is outside the edge but on the extension line, so take the closest endpoint
					dist_from_line, proj_point_from_line = find_closest_vertex(x1, y1, x2, y2, ptx, pty)
				end
			elseif abs(y1 - y2) < tol  # edge is horizontal
				x0 = ptx
				y0 = y1
				if min(x1, x2) < x0 < max(x1, x2)
					dist_from_line = abs(y0 - pty)
					proj_point_from_line = [x0, y0]
				else
					dist_from_line, proj_point_from_line = find_closest_vertex(x1, y1, x2, y2, ptx, pty)
				end
			else  # edge is neither vertical nor horizontal
				m1 = (y2 - y1) / (x2 - x1)
				m2 = -1 / m1
				# Point of projection on the extended line:
				x0 = (m1*x1 - y1 - m2*ptx + pty) / (m1 - m2)
				y0 = m2*(x0 - ptx) + pty
				if min(x1, x2) < x0 < max(x1, x2) && min(y1, y2) < y0 < max(y1, y2)
					dist_from_line = sqrt((x0 - ptx)^2 + (y0 - pty)^2)
					proj_point_from_line = [x0, y0]
				else
					dist_from_line, proj_point_from_line = find_closest_vertex(x1, y1, x2, y2, ptx, pty)
				end
			end

			if dist_from_line < dist
				dist = dist_from_line
				proj_point = proj_point_from_line
			end
		end
	end
	
	return dist, proj_point
end

# ╔═╡ 6082694f-5c6c-46f6-ba0a-8a70b19472a2
let
	# P_outside = [1.5, 0.3]
	point_outside = [-0.2, 0.1]
	vertex1 = [0, 0]
	vertex2 = [1, 1]
	vertex3 = [1, 0]
	t_vertices = [vertex1, vertex2, vertex3]
	t_sides = [
		interpolate_edge(vertex1, vertex2),
		interpolate_edge(vertex2, vertex3),
		interpolate_edge(vertex1, vertex3)
	]
	
	d_hausdorff, proj_point = nearest_triangle_point(t_vertices, point_outside)

	f = Figure()
	
	ax = Axis(f[1, 1])
	lines!(ax, t_sides[1][1], t_sides[1][2])
	lines!(ax, t_sides[2][1], t_sides[2][2])
	lines!(ax, t_sides[3][1], t_sides[3][2])
	
	scatter!(point_outside[1], point_outside[2]; color = :black)
	scatter!(proj_point[1], proj_point[2]; color = :red)
	lines!(
		[point_outside[1], proj_point[1]],
		[point_outside[2], proj_point[2]];
		linestyle = :dash
	)

	f
end

# ╔═╡ dd10f461-1b59-44d9-b669-16f5a0ff2e40
md"""
## `make_triangle`
"""

# ╔═╡ adef7ca0-60e0-401b-a753-83df3d51c8fc
function make_triangle(
	material1_density, material2_density;
	monoenergetic::Bool = true)
	
	if monoenergetic
		water_coord = [
			BkgWater70200[1],
			BkgWater150140[1]
		]
		material1_coord = [
			Ca70200[material1_density][1], 
			Ca150140[material1_density][1]
		]
		material2_coord = [
			I70200[material2_density][1],
			I150140[material2_density][1]
		]
		d_cutoff = maximum([
			BkgWater70200[2],
			BkgWater150140[2],
			Ca70200[material1_density][2],
			Ca150140[material1_density][2],
			I70200[material2_density][2],
			I150140[material2_density][2]
		])
	else
		water_coord = [
			BkgWater80550[1],
			BkgWater135200[1]
		]
		material1_coord = [
			Ca80550[material1_density][1],
			Ca135200[material1_density][1]
		]
		material2_coord = [
			I80550[material2_density][1],
			I135200[material2_density][1]
		]
		d_cutoff = maximum([
			BkgWater80550[2],
			BkgWater135200[2],
			Ca80550[material1_density][2],
			Ca135200[material1_density][2],
			I80550[material2_density][2],
			I135200[material2_density][2]
		])
	end

	line1 = interpolate_edge(water_coord, material1_coord)
	line2 = interpolate_edge(material1_coord, material2_coord)
	line3 = interpolate_edge(water_coord, material2_coord)

	return [line1, line2, line3], [water_coord, material1_coord, material2_coord], d_cutoff
end

# ╔═╡ e3c1cb7b-6dbd-4585-b76a-859a5cdeee90
begin
    Ca70200_keys = keys(Ca70200)
    I70200_keys = keys(I70200)

	vertices = [(c, i) for c in Ca70200_keys for i in I70200_keys]
end

# ╔═╡ 301ae9b0-fbda-4310-84d9-edd036dc881e
list_of_triangles = [make_triangle(vertex[1], vertex[2]) for vertex in vertices]

# ╔═╡ 4f8b68cb-9df5-4777-90ae-f652d7bed523
# (Ca300, I20): [[0.077, -0.431], [867, 523], [476.491, 106.228]]
t1_vertices = list_of_triangles[28][2]

# ╔═╡ 231c286c-2e6c-476b-9cf0-15200c12d59a
# 5 (Ca densities) * 7 (I densities) = 35 all possible triangles
length(list_of_triangles)

# ╔═╡ 7ac327f5-5602-4f8f-8b38-12e5ace84d9a
md"""
## `matrix_calculations`
"""

# ╔═╡ 93a76105-4cda-4b4b-ace7-c6c03932801d
function matrix_calculations(triangle_vertices, point)
    water_coord, Ca_coord, I_coord = triangle_vertices
    ptx, pty = point

    A = [water_coord[1] Ca_coord[1] I_coord[1];
         water_coord[2] Ca_coord[2] I_coord[2];
         1 1 1]
    b = [ptx, pty, 1]
    det_A = det(A)

    if abs(det_A) > 1e-6
        alpha = A \ b
        alpha[abs.(alpha) .< 1e-4] .= 0
        return (A, det_A, alpha)
    else
        return (A, 0, zeros(3))
    end
end

# ╔═╡ 506fdf0f-0533-490c-a226-33859518eac3
t1_sides = list_of_triangles[28][1]

# ╔═╡ 0101f4cd-2d46-4841-8f4f-cdcdb19213e0
begin
	line1 = t1_sides[1]
	line2 = t1_sides[2]
	line3 = t1_sides[3]
end

# ╔═╡ c1c7299b-f05d-4656-bed3-32d5b9233808
begin
	point_inside = [400, 180]
	point_outside = [800, 200]
end

# ╔═╡ bc0f6aff-3d7a-436f-b50f-b2f16c08c806
matrix_calculations(t1_vertices, point_inside)

# ╔═╡ e2695055-76cd-4aae-85f7-89de67aaf7ae
matrix_calculations(t1_vertices, point_outside)

# ╔═╡ 90ef4667-ba73-4a43-b9ef-3921bdad693b
md"""
## `mask_image`
"""

# ╔═╡ d7da1045-d598-4f0e-9414-da282b71cb4b
function mask_image(image; center = [256, 256], radius = 200)
	x = repeat(1:512, 1, 512)
	y = repeat((1:512)', 512, 1)
	distance = sqrt.((x .- center[1]).^2 .+ (y .- center[2]).^2)
	mask = distance .< radius
	image_copy = copy(image)
	image_copy[.!mask] .= 0

	return image_copy
end

# ╔═╡ 5ca24cf2-e4fe-463e-8b09-7c4ffcf48ff3
let
	low_e_dcm_rescaled = dcm_parse(path_70_keV)
	rescale!(low_e_dcm_rescaled)
	low_e_dcm_rescaled = low_e_dcm_rescaled.PixelData

	high_e_dcm_rescaled = dcm_parse(path_150_keV)
	rescale!(high_e_dcm_rescaled)
	high_e_dcm_rescaled = high_e_dcm_rescaled.PixelData
	
	low_e_dcm_masked = mask_image(low_e_dcm_rescaled)
	high_e_dcm_masked = mask_image(high_e_dcm_rescaled)

	f = Figure(size = (1000, 1000))
	ax = Axis(
		f[1, 1],
		title = "Low E Normal"
	)
	heatmap!(low_e_dcm_rescaled; colormap = :grays)

	ax = Axis(
		f[1, 2],
		title = "Low E Masked"
	)
	heatmap!(low_e_dcm_masked; colormap = :grays)

	ax = Axis(
		f[2, 1],
		title = "High E Normal"
	)
	heatmap!(high_e_dcm_rescaled; colormap = :grays)

	ax = Axis(
		f[2, 2],
		title = "High E Masked"
	)
	heatmap!(high_e_dcm_masked; colormap = :grays)
	f
end

# ╔═╡ b5f74b05-a95d-429a-9de1-1c8562f3eff2
md"""
## `register`
"""

# ╔═╡ eceec935-01a9-4f14-a51b-b1d7db57d53a
function register(
	fixed_image::AbstractArray, moving_image::AbstractArray;
	num_iterations = 10)
	
	fixed_image_pyarr = PyArray(fixed_image)
	fixed_image_sitk= sitk.GetImageFromArray(fixed_image_pyarr)
	fixed_image_sitk = sitk.Cast(fixed_image_sitk, sitk.sitkFloat32)

	moving_image_pyarr = PyArray(moving_image)
	moving_image_sitk = sitk.GetImageFromArray(moving_image_pyarr)
	moving_image_sitk = sitk.Cast(moving_image_sitk, sitk.sitkFloat32)

	# Perform Demons registration
	demons_filter = sitk.DemonsRegistrationFilter()
	demons_filter.SetNumberOfIterations(num_iterations)
	displacement_field = demons_filter.Execute(fixed_image_sitk, moving_image_sitk)

	# Create a displacement field transform and apply it to the entire moving image
	displacement_transform = sitk.DisplacementFieldTransform(displacement_field)
	moving_registered_sitk = sitk.Resample(moving_image_sitk, fixed_image_sitk, displacement_transform, sitk.sitkLinear, 0.0, moving_image_sitk.GetPixelID())
	moving_registered_pyarr = sitk.GetArrayFromImage(moving_registered_sitk)
	return pyconvert(Array, moving_registered_pyarr)
end

# ╔═╡ e21b0851-ed5b-444a-a322-a8636e7f59a8
begin
	low_e_dcm_rescaled = dcm_parse(path_70_keV)
	rescale!(low_e_dcm_rescaled)
	low_e_dcm_rescaled = low_e_dcm_rescaled.PixelData

	high_e_dcm_rescaled = dcm_parse(path_150_keV)
	rescale!(high_e_dcm_rescaled)
	high_e_dcm_rescaled = high_e_dcm_rescaled.PixelData
	
	low_e_dcm_masked = mask_image(low_e_dcm_rescaled)
	high_e_dcm_masked = mask_image(high_e_dcm_rescaled)
end;

# ╔═╡ b9fbea32-5500-4d37-b26b-c55e56d20f8f
high_e_dcm_masked_registered = register(low_e_dcm_masked, high_e_dcm_masked; num_iterations = 1000);

# ╔═╡ 6ad52280-b8dc-4b17-97ac-2ed68e659d49
let
	f = Figure()
	ax = Axis(
		f[1, 1],
		title = "Low E Fixed"
	)
	heatmap!(low_e_dcm_masked)

	ax = Axis(
		f[1, 2],
		title = "High E Registered"
	)
	heatmap!(high_e_dcm_masked_registered)

	ax = Axis(
		f[2, 1],
		title = "Unregistered Error"
	)
	heatmap!(low_e_dcm_masked .- high_e_dcm_masked)

	ax = Axis(
		f[2, 2],
		title = "Registered Error"
	)
	heatmap!(low_e_dcm_masked .- high_e_dcm_masked_registered)
	f
end

# ╔═╡ d9218299-dc81-402c-887e-a3593a660cd4
md"""
# MMD Algorithm
"""

# ╔═╡ 619f81c2-2cec-4ac5-b7f9-def57106a176
function mmd_algorithm(low_e_image, high_e_image, list_of_triangles, vertices)
	total_pixels = size(low_e_image, 1) * size(low_e_image, 2)

	# Initialize concentration maps
	ca_maps = OrderedDict(
		vertex[1] => zeros(size(low_e_image)) for vertex in vertices
	)
	i_maps = OrderedDict(
		vertex[2] => zeros(size(low_e_image)) for vertex in vertices
	)
	water_map = zeros(size(low_e_image))

	for px in 1:total_pixels
		row = div(px - 1, size(low_e_image, 1)) + 1
		col = mod(px - 1, size(low_e_image, 2)) + 1

		# Skip the region outside the phantom
		if low_e_image[row, col] == 0 && high_e_image[row, col] == 0
			continue
		end

		pt = [low_e_image[row, col], high_e_image[row, col]] # measured HU point
		dist = Inf
		soln = [0.0, 0.0, 0.0]

		local vertex_idx
		for (idx, triangle) in enumerate(list_of_triangles)
			t_sides, t_vertices, d_cutoff = triangle
			a, det_a, alpha = matrix_calculations(t_vertices, pt)

			if any(alpha .< 0)
				d_hausdorff, proj_point = nearest_triangle_point(t_vertices, pt)
				if d_hausdorff < dist
					dist = d_hausdorff
					a, det_a, alpha = matrix_calculations(t_vertices, proj_point)
					soln = alpha
					vertex_idx = idx
				else
					continue
				end
			else
				soln = alpha
				vertex_idx = idx
				break
			end
		end

		# Update maps
		if all(soln .== [0.0, 0.0, 0.0])
			water_map[row, col] = 0
			for ca_density in keys(ca_maps)
				ca_maps[ca_density][row, col] = 0
			end
			for i_density in keys(i_maps)
				i_maps[i_density][row, col] = 0
			end
		else
			water_map[row, col] = soln[1]
			ca_density, i_density = vertices[vertex_idx]
			ca_maps[ca_density][row, col] = soln[2]
			i_maps[i_density][row, col] = soln[3]
		end

	end

	return water_map, ca_maps, i_maps
end

# ╔═╡ 47785370-5b32-464e-bde5-f359215a3395
function calculate_final_density_maps(water_map, ca_maps, i_maps)
	# Initialize the final density maps
	final_ca_image = zeros(size(water_map))
	final_i_image = zeros(size(water_map))

	# Calculate total densities
	total_ca_density = sum(keys(ca_maps))
	total_i_density = sum(keys(i_maps))

	# Combine Ca maps
	for (ca_density, ca_map) in ca_maps
		final_ca_image .+= ca_density .* ca_map
	end

	# Combine I maps
	for (i_density, i_map) in i_maps
		final_i_image .+= i_density .* i_map
	end

	# Normalize the final images by the total densities
	final_ca_image ./= total_ca_density
	final_i_image ./= total_i_density

	return final_ca_image, final_i_image
end

# ╔═╡ ab8674df-58bc-40e8-871c-772723bc30cb
md"""
## Unregistered
"""

# ╔═╡ e3c007c5-e536-4b8f-ac05-f80396b630aa
water_map, ca_maps, i_maps = mmd_algorithm(low_e_dcm_masked, high_e_dcm_masked, list_of_triangles, vertices);

# ╔═╡ 84ce2c06-5222-45cb-804a-a08229887c94
final_ca_image, final_i_image = calculate_final_density_maps(water_map, ca_maps, i_maps);

# ╔═╡ e77da98a-21db-41e0-bf0a-fcee5bd9dad1
md"""
## Registered
"""

# ╔═╡ 72610fd6-770e-4129-b548-a06283b8d2e8
water_map_reg, ca_maps_reg, i_maps_reg = mmd_algorithm(low_e_dcm_masked, high_e_dcm_masked_registered, list_of_triangles, vertices);

# ╔═╡ 5281ab62-bc95-422d-8513-8c72882035d7
final_ca_image_reg, final_i_image_reg = calculate_final_density_maps(water_map_reg, ca_maps_reg, i_maps_reg);

# ╔═╡ da375fb5-7d6e-4339-84ba-e803c3a62680
let
	f = Figure(size = (800, 800))
	ax = Axis(
		f[1, 1],
		title = "Water Map"
	)
	heatmap!(water_map; colormap = :jet)

	ax = Axis(
		f[2, 1],
		title = "Calcium Map"
	)
	heatmap!(final_ca_image; colormap = :jet)

	ax = Axis(
		f[3, 1],
		title = "Iodine Map"
	)
	heatmap!(final_i_image; colormap = :jet)

	ax = Axis(
		f[1, 2],
		title = "Water Map Registered"
	)
	heatmap!(water_map_reg; colormap = :jet)

	ax = Axis(
		f[2, 2],
		title = "Calcium Map Registered"
	)
	heatmap!(final_ca_image_reg; colormap = :jet)

	ax = Axis(
		f[3, 2],
		title = "Iodine Map Registered"
	)
	heatmap!(final_i_image_reg; colormap = :jet)
	f
end

# ╔═╡ Cell order:
# ╟─dd80ed67-d7c9-4551-a35f-c068f095d8a9
# ╟─c9179fde-37af-4f7d-a87b-2301ea96b643
# ╠═ea478b7f-96d2-4802-9468-33f4bdf7e9d7
# ╠═21d854cf-06a3-4214-a689-f57be77e84be
# ╠═45e73788-4b7c-11ef-30bd-e161fe50a735
# ╠═2a149318-4c73-4a57-a819-79f70adc1700
# ╠═e8af4ce9-72c7-4017-bd92-6fd3a419e973
# ╠═bed00701-b564-4a26-90eb-0f9132da0ce3
# ╠═d5dfce22-566f-4361-8cfc-237700840676
# ╠═7ec95459-8e70-4a06-9466-262b97fb871f
# ╠═1a04d22f-047b-42e7-84a3-8bd0fd222743
# ╠═d01cf5fa-a788-4efb-ad47-24e2a014f2f8
# ╠═42b57d75-3b39-4c61-995d-8d0b20a6a7df
# ╠═a88d38c6-c85a-414d-911d-146e20c7fd29
# ╠═21415643-b7b9-420b-b7ef-a7c28f40e7e1
# ╠═dab48e46-c732-44ff-9bed-391f94e4ca43
# ╟─b06d9b63-67b3-413b-96cf-9035678c721c
# ╠═07faae20-bd30-4412-bd57-27c95c8c8740
# ╠═3a6b40c4-a25a-4e28-9375-f06fe15f6cfc
# ╠═a30419d5-2ca0-4af8-bc2b-676dcda84f60
# ╠═ab512723-b5e5-496c-99c6-2277f70bb0a0
# ╟─179fb57b-10c6-44c0-a47e-4d0d96f9f71b
# ╟─9830994f-2f85-48fd-aae2-b70330019658
# ╠═72c7e682-58db-4577-be77-10d5376604e8
# ╟─0f7faad8-73cb-4417-8048-fca690c3e2db
# ╠═fd67cdc9-5e6f-4771-8d87-6a2183f3f93d
# ╟─0c8c12fd-c439-44c7-a1e1-7f5f0674f2d3
# ╟─935739af-dcde-49e4-bdc0-0875c7335587
# ╟─84ea30d1-f3b3-4a37-ae46-b991bd71072b
# ╠═ba5ef8fe-4f03-4f01-8acf-27c990bed4d1
# ╠═b484d333-a9b8-4a8f-b325-6a360675f271
# ╠═b57e7f8e-0cb5-4591-bee4-6b755f8e0134
# ╠═3ac3b86e-20a3-48cb-bf17-8593353e9210
# ╟─2ba6333f-850b-47b7-b39c-60c2da627cb5
# ╠═d59e48b9-7e94-4598-a590-2ead23389c81
# ╠═bc56c4ec-c232-4789-82d3-c290bbad77e1
# ╟─fcb9cdfe-e516-4f34-afc9-bb5f66cf6b08
# ╠═cba2ace6-90ff-4476-b576-0f87c0da0614
# ╠═6082694f-5c6c-46f6-ba0a-8a70b19472a2
# ╟─dd10f461-1b59-44d9-b669-16f5a0ff2e40
# ╠═adef7ca0-60e0-401b-a753-83df3d51c8fc
# ╠═e3c1cb7b-6dbd-4585-b76a-859a5cdeee90
# ╠═4f8b68cb-9df5-4777-90ae-f652d7bed523
# ╠═301ae9b0-fbda-4310-84d9-edd036dc881e
# ╠═231c286c-2e6c-476b-9cf0-15200c12d59a
# ╟─7ac327f5-5602-4f8f-8b38-12e5ace84d9a
# ╠═93a76105-4cda-4b4b-ace7-c6c03932801d
# ╠═506fdf0f-0533-490c-a226-33859518eac3
# ╠═0101f4cd-2d46-4841-8f4f-cdcdb19213e0
# ╠═c1c7299b-f05d-4656-bed3-32d5b9233808
# ╠═bc0f6aff-3d7a-436f-b50f-b2f16c08c806
# ╠═e2695055-76cd-4aae-85f7-89de67aaf7ae
# ╟─90ef4667-ba73-4a43-b9ef-3921bdad693b
# ╠═d7da1045-d598-4f0e-9414-da282b71cb4b
# ╟─5ca24cf2-e4fe-463e-8b09-7c4ffcf48ff3
# ╟─b5f74b05-a95d-429a-9de1-1c8562f3eff2
# ╠═eceec935-01a9-4f14-a51b-b1d7db57d53a
# ╠═e21b0851-ed5b-444a-a322-a8636e7f59a8
# ╠═b9fbea32-5500-4d37-b26b-c55e56d20f8f
# ╟─6ad52280-b8dc-4b17-97ac-2ed68e659d49
# ╟─d9218299-dc81-402c-887e-a3593a660cd4
# ╠═619f81c2-2cec-4ac5-b7f9-def57106a176
# ╠═47785370-5b32-464e-bde5-f359215a3395
# ╟─ab8674df-58bc-40e8-871c-772723bc30cb
# ╠═e3c007c5-e536-4b8f-ac05-f80396b630aa
# ╠═84ce2c06-5222-45cb-804a-a08229887c94
# ╟─e77da98a-21db-41e0-bf0a-fcee5bd9dad1
# ╠═72610fd6-770e-4129-b548-a06283b8d2e8
# ╠═5281ab62-bc95-422d-8513-8c72882035d7
# ╟─da375fb5-7d6e-4339-84ba-e803c3a62680
