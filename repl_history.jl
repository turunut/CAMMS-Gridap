# time: 2023-03-08 23:37:06 CET
# mode: julia
	using Pkg
# time: 2023-03-08 23:37:16 CET
# mode: julia
	Pkg.instantiate()
# time: 2023-03-08 23:38:02 CET
# mode: julia
	import Pkg; Pkg.add("Gridap")
# time: 2023-03-08 23:42:04 CET
# mode: julia
	import Pkg; Pkg.add("GLMakie")
# time: 2023-03-08 23:53:38 CET
# mode: julia
	import Pkg; Pkg.add("Makie")
# time: 2023-03-08 23:57:29 CET
# mode: pkg
	add Revise
# time: 2023-03-08 23:57:56 CET
# mode: pkg
	add OhMyREPL
# time: 2023-03-08 23:58:08 CET
# mode: julia
	exit()
# time: 2023-03-09 00:01:48 CET
# mode: julia
	x = [1,2,3]
# time: 2023-03-09 05:45:42 CET
# mode: julia
	get_face_labeling(model)
# time: 2023-03-09 05:59:46 CET
# mode: julia
	import Pkg; Pkg.add("Metis")
# time: 2023-03-09 06:02:10 CET
# mode: julia
	import Pkg; Pkg.add("GridapDistributed")
# time: 2023-03-09 06:10:16 CET
# mode: julia
	status --outdated -m
# time: 2023-03-09 06:10:28 CET
# mode: julia
	import Pkg; Pkg.add("PartitionedArrays")
# time: 2023-03-09 06:13:22 CET
# mode: julia
	Pkg.build("GridapGmsh")
# time: 2023-03-09 06:13:45 CET
# mode: julia
	exit()
# time: 2023-03-09 06:15:56 CET
# mode: julia
	exit()
# time: 2023-03-09 06:16:03 CET
# mode: julia
	Pkg.build("GridapGmsh")
# time: 2023-03-09 06:16:10 CET
# mode: julia
	using Pkg
# time: 2023-03-09 06:16:12 CET
# mode: julia
	Pkg.build("GridapGmsh")
# time: 2023-03-09 06:16:48 CET
# mode: julia
	Pkg.instatiate()
# time: 2023-03-09 06:17:10 CET
# mode: julia
	Pkg.add("GridapGmsh")
# time: 2023-03-09 06:38:14 CET
# mode: julia
	paraview
# time: 2023-03-10 03:07:03 CET
# mode: julia
	exit()
# time: 2023-03-10 03:07:12 CET
# mode: julia
	using Pkg
# time: 2023-03-10 03:07:26 CET
# mode: julia
	pkg> add GridapGmsh
# time: 2023-03-10 03:07:38 CET
# mode: julia
	pkg
# time: 2023-03-10 03:07:43 CET
# mode: julia
	Pkg
# time: 2023-03-10 03:08:54 CET
# mode: julia
	GridapGmsh
# time: 2023-03-10 03:09:00 CET
# mode: julia
	using GridapGmsh
# time: 2023-03-10 03:09:39 CET
# mode: julia
	GridapGmsh
# time: 2023-03-10 03:11:46 CET
# mode: julia
	pkg> instantiate
# time: 2023-03-10 04:01:29 CET
# mode: julia
	gid
# time: 2023-03-10 06:59:52 CET
# mode: julia
	Γ
# time: 2023-03-10 06:59:58 CET
# mode: julia
	dΓ
# time: 2023-03-10 07:00:08 CET
# mode: julia
	measure(dΓ\)
# time: 2023-03-10 07:00:10 CET
# mode: julia
	measure(dΓ)
# time: 2023-03-10 07:09:49 CET
# mode: julia
	test = get_fe_basis(U)
# time: 2023-03-10 07:10:43 CET
# mode: julia
	trial = get_trial_fe_basis(V)
# time: 2023-03-10 07:10:48 CET
# mode: julia
	trial = get_trial_fe_basis(U)
# time: 2023-03-10 07:10:52 CET
# mode: julia
	U
# time: 2023-03-10 07:10:58 CET
# mode: julia
	U.spaces
# time: 2023-03-10 07:11:05 CET
# mode: julia
	V.spaces
# time: 2023-03-10 07:11:30 CET
# mode: julia
	test = get_fe_basis(V)
# time: 2023-03-10 07:11:36 CET
# mode: julia
	test = get_fe_basis(U)
# time: 2023-03-10 07:12:26 CET
# mode: julia
	test
# time: 2023-03-10 07:12:31 CET
# mode: julia
	test.single_fields
# time: 2023-03-10 07:12:38 CET
# mode: julia
	test.single_fields[1]
# time: 2023-03-10 07:12:42 CET
# mode: julia
	test.single_fields[2]
# time: 2023-03-10 07:14:20 CET
# mode: julia
	uh
# time: 2023-03-10 07:14:30 CET
# mode: julia
	uh.single_fe_functions
# time: 2023-03-10 07:14:33 CET
# mode: julia
	uh.single_fe_functions[1]
# time: 2023-03-13 07:18:59 CET
# mode: julia
	v = get_fe_basis(Vu)
# time: 2023-03-13 07:19:12 CET
# mode: julia
	u = get_trial_fe_basis(Uu)
# time: 2023-03-13 07:20:10 CET
# mode: julia
	Gridap.CellData.get_data(u)
# time: 2023-03-13 07:20:30 CET
# mode: julia
	u_data = Gridap.CellData.get_data(u)
# time: 2023-03-13 07:20:34 CET
# mode: julia
	u_data[1]
# time: 2023-03-13 07:21:14 CET
# mode: julia
	v_data = Gridap.CellData.get_data(v)
# time: 2023-03-13 07:21:21 CET
# mode: julia
	v_data[1]
# time: 2023-03-13 07:23:06 CET
# mode: julia
	σ(u)
# time: 2023-03-13 07:23:37 CET
# mode: julia
	du = ε(u)
# time: 2023-03-13 07:24:05 CET
# mode: julia
	du_data = Gridap.CellData.get_data(du)
# time: 2023-03-13 07:24:16 CET
# mode: julia
	du_data[1]
# time: 2023-03-13 07:25:32 CET
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-03-13 07:25:55 CET
# mode: julia
	pts = get_cell_points(dΩ)
# time: 2023-03-13 07:26:06 CET
# mode: julia
	pts.cell_ref_point[1]
# time: 2023-03-13 07:26:32 CET
# mode: julia
	quad = get_cell_quadrature(dΩ)
# time: 2023-03-13 07:27:50 CET
# mode: julia
	pts = get_cell_quadrature(dΩ)
# time: 2023-03-13 07:28:20 CET
# mode: julia
	dΩ.quad
# time: 2023-03-13 07:28:26 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-13 07:28:37 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-13 07:28:56 CET
# mode: julia
	w = quad.cell_weight
# time: 2023-03-13 07:29:01 CET
# mode: julia
	w[1]
# time: 2023-03-13 07:29:34 CET
# mode: julia
	du(pts)
# time: 2023-03-13 07:29:42 CET
# mode: julia
	dux = du(pts)
# time: 2023-03-13 07:29:45 CET
# mode: julia
	dux[1]
# time: 2023-03-13 07:30:34 CET
# mode: julia
	dux[1][1,1,1]
# time: 2023-03-13 07:33:12 CET
# mode: julia
	du = σ∘ε(u)
# time: 2023-03-13 07:33:22 CET
# mode: julia
	dux = du(pts)
# time: 2023-03-13 07:33:28 CET
# mode: julia
	dux[1]
# time: 2023-03-13 07:33:44 CET
# mode: julia
	dux[1][1,1,1]
# time: 2023-03-13 07:40:58 CET
# mode: julia
	du = get_x(σ∘ε(u))
# time: 2023-03-13 07:41:23 CET
# mode: julia
	du = get_x∘σ∘ε(u)
# time: 2023-03-13 07:41:28 CET
# mode: julia
	dux = du(pts)
# time: 2023-03-13 07:41:31 CET
# mode: julia
	dux[1][1,1,1]
# time: 2023-03-13 08:02:21 CET
# mode: julia
	tr(ε)
# time: 2023-03-13 08:02:31 CET
# mode: julia
	tr(ε0)
# time: 2023-03-13 08:03:16 CET
# mode: julia
	tr(ε0(1.0))
# time: 2023-03-13 08:03:50 CET
# mode: julia
	one(ε0(1.0))
# time: 2023-03-13 08:04:51 CET
# mode: julia
	u.cell_basis[1]
# time: 2023-03-13 08:04:56 CET
# mode: julia
	u.cell_basis[1][1]
# time: 2023-03-13 08:05:06 CET
# mode: julia
	ε(u.cell_basis[1][1])
# time: 2023-03-13 08:05:23 CET
# mode: julia
	ε(u.cell_basis[1][1])(VectorValue(0.0,0.0))
# time: 2023-03-13 08:06:49 CET
# mode: julia
	SymTensorValue(2.0,2.0,2.0)
# time: 2023-03-13 23:14:23 CET
# mode: julia
	test = get_fe_basis(U)
# time: 2023-03-13 23:14:30 CET
# mode: julia
	U
# time: 2023-03-13 23:14:44 CET
# mode: julia
	MultiFieldFESpace
# time: 2023-03-13 23:26:47 CET
# mode: julia
	du = σ∘ε(u)
# time: 2023-03-14 02:55:29 CET
# mode: julia
	model
# time: 2023-03-15 01:22:01 CET
# mode: julia
	\varepsilon
# time: 2023-03-15 01:22:04 CET
# mode: julia
	ε
# time: 2023-03-15 01:22:42 CET
# mode: julia
	symmetric_gradient
# time: 2023-03-15 01:23:00 CET
# mode: julia
	symmetric_gradient()
# time: 2023-03-15 02:26:24 CET
# mode: julia
	ε0
# time: 2023-03-15 02:32:23 CET
# mode: julia
	ε
# time: 2023-03-15 02:32:26 CET
# mode: julia
	ε0
# time: 2023-03-15 02:32:50 CET
# mode: julia
	symmetric_gradient
# time: 2023-03-15 02:34:58 CET
# mode: julia
	(ε(u),tags)
# time: 2023-03-15 02:35:04 CET
# mode: julia
	(ε,tags)
# time: 2023-03-15 02:35:15 CET
# mode: julia
	(ε0,tags)
# time: 2023-03-15 02:35:26 CET
# mode: julia
	(symTen,tags)
# time: 2023-03-15 02:35:40 CET
# mode: julia
	(symTen(ε0),tags)
# time: 2023-03-15 02:36:28 CET
# mode: julia
	(symmetric_gradient(ε0),tags)
# time: 2023-03-15 02:36:46 CET
# mode: julia
	(ε0,tags)
# time: 2023-03-15 02:37:28 CET
# mode: julia
	(ε, tags)
# time: 2023-03-16 02:52:04 CET
# mode: julia
	\varepsilon
# time: 2023-03-16 02:52:06 CET
# mode: julia
	ε
# time: 2023-03-16 05:31:59 CET
# mode: julia
	κ(u)
# time: 2023-03-16 05:32:05 CET
# mode: julia
	aux = κ(u)
# time: 2023-03-16 05:32:16 CET
# mode: julia
	get_data(aux)
# time: 2023-03-16 05:33:11 CET
# mode: julia
	κ(u)
# time: 2023-03-16 05:34:02 CET
# mode: julia
	pts = get_cell_quadrature(dΩ)
# time: 2023-03-16 05:34:08 CET
# mode: julia
	pts = get_cell_quadratures(dΩ)
# time: 2023-03-16 05:34:15 CET
# mode: julia
	using Gridap.CellData
# time: 2023-03-16 05:34:16 CET
# mode: julia
	pts = get_cell_quadratures(dΩ)
# time: 2023-03-16 05:34:19 CET
# mode: julia
	pts = get_cell_quadrature(dΩ)
# time: 2023-03-16 05:34:32 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 05:34:45 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 05:34:51 CET
# mode: julia
	cf(pts)
# time: 2023-03-16 05:35:02 CET
# mode: julia
	cf(pts)[1]
# time: 2023-03-16 05:43:47 CET
# mode: julia
	pts = get_cell_quadratures(dΩ)
# time: 2023-03-16 05:44:07 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 05:44:52 CET
# mode: julia
	cf(pts)
# time: 2023-03-16 05:44:57 CET
# mode: julia
	cf(pts)[1]
# time: 2023-03-16 05:45:04 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 05:45:09 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 05:45:16 CET
# mode: julia
	cf(pts)[1]
# time: 2023-03-16 06:36:10 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 06:36:21 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 06:36:28 CET
# mode: julia
	cf(pts)[1]
# time: 2023-03-16 06:36:58 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 06:37:01 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 06:37:03 CET
# mode: julia
	cf(pts)[1]
# time: 2023-03-16 06:38:19 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 06:38:21 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 06:38:26 CET
# mode: julia
	cf(pts)
# time: 2023-03-16 06:38:34 CET
# mode: julia
	cf(pts)[1]
# time: 2023-03-16 06:39:37 CET
# mode: julia
	cp(pts)
# time: 2023-03-16 06:39:41 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 06:39:44 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 06:39:47 CET
# mode: julia
	cp(pts)
# time: 2023-03-16 06:39:58 CET
# mode: julia
	cp(pts)[1]
# time: 2023-03-16 06:43:04 CET
# mode: julia
	cq(pts)
# time: 2023-03-16 06:43:08 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 06:43:12 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 06:43:15 CET
# mode: julia
	cq(pts)
# time: 2023-03-16 06:43:22 CET
# mode: julia
	cq(pts)[q]
# time: 2023-03-16 06:43:24 CET
# mode: julia
	cq(pts)[1]
# time: 2023-03-16 06:44:09 CET
# mode: julia
	cp(pts]
# time: 2023-03-16 06:44:14 CET
# mode: julia
	cp(pts)
# time: 2023-03-16 06:44:18 CET
# mode: julia
	cp(pts)[1]
# time: 2023-03-16 07:02:33 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 07:02:35 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 07:02:38 CET
# mode: julia
	cp(pts)
# time: 2023-03-16 07:02:41 CET
# mode: julia
	cp(pts)[1]
# time: 2023-03-16 07:03:14 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 07:03:15 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 07:03:17 CET
# mode: julia
	cp(pts)
# time: 2023-03-16 07:05:40 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 07:05:42 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 07:05:44 CET
# mode: julia
	cp(pts)
# time: 2023-03-16 07:05:48 CET
# mode: julia
	cp(pts)[1]
# time: 2023-03-16 07:06:08 CET
# mode: julia
	cp(pts)[1,1]
# time: 2023-03-16 07:06:23 CET
# mode: julia
	cp(pts)[1,1,1]
# time: 2023-03-16 07:06:26 CET
# mode: julia
	cp(pts)[1,1,1,1]
# time: 2023-03-16 07:06:32 CET
# mode: julia
	cp(pts)[1][1]
# time: 2023-03-16 07:06:37 CET
# mode: julia
	cp(pts)[1][2]
# time: 2023-03-16 07:06:41 CET
# mode: julia
	cp(pts)[2][1]
# time: 2023-03-16 07:06:44 CET
# mode: julia
	cp(pts)[2][2]
# time: 2023-03-16 23:05:04 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 23:05:08 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:05:12 CET
# mode: julia
	cp(pts)[1]
# time: 2023-03-16 23:05:17 CET
# mode: julia
	cp(pts)
# time: 2023-03-16 23:05:25 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:05:29 CET
# mode: julia
	cp(pts)
# time: 2023-03-16 23:05:57 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 23:05:59 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:06:10 CET
# mode: julia
	c(pts)
# time: 2023-03-16 23:06:13 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-16 23:08:20 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 23:08:22 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:08:25 CET
# mode: julia
	c(pts)
# time: 2023-03-16 23:08:30 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-16 23:15:26 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 23:15:27 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:15:32 CET
# mode: julia
	c(pts)
# time: 2023-03-16 23:15:35 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-16 23:19:49 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 23:19:50 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:19:53 CET
# mode: julia
	c(pts)
# time: 2023-03-16 23:19:56 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-16 23:22:02 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-16 23:22:03 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:22:05 CET
# mode: julia
	c(pts)
# time: 2023-03-16 23:22:09 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-16 23:23:21 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:23:22 CET
# mode: julia
	c(pts)
# time: 2023-03-16 23:23:26 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-16 23:23:28 CET
# mode: julia
	c(pts)
# time: 2023-03-16 23:23:59 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-17 00:10:57 CET
# mode: julia
	uw(uh)
# time: 2023-03-17 00:11:07 CET
# mode: julia
	uw(uh)[1]
# time: 2023-03-17 00:50:50 CET
# mode: julia
	γ(uh)
# time: 2023-03-17 00:51:20 CET
# mode: julia
	get_cell_points(γ(uh))
# time: 2023-03-17 00:52:13 CET
# mode: julia
	dv = γ(uh)
# time: 2023-03-17 00:52:24 CET
# mode: julia
	Gridap.CellData.get_data(dv)
# time: 2023-03-17 00:52:44 CET
# mode: julia
	dv_data = Gridap.CellData.get_data(dv)
# time: 2023-03-17 00:52:50 CET
# mode: julia
	dv_data[1]
# time: 2023-03-17 00:52:56 CET
# mode: julia
	dv_data[1][1]
# time: 2023-03-17 01:03:50 CET
# mode: julia
	import Pkg; Pkg.add("Interpolations")
# time: 2023-03-17 01:06:36 CET
# mode: julia
	import Pkg; Pkg.add("QuadGK")
# time: 2023-03-17 01:14:12 CET
# mode: julia
	import Pkg; Pkg.add("Pardiso")
# time: 2023-03-17 02:08:26 CET
# mode: julia
	u = get_trial_fe_basis(U)
# time: 2023-03-17 02:08:53 CET
# mode: julia
	v = get_fe_basis(U)
# time: 2023-03-17 02:20:52 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-17 02:21:01 CET
# mode: julia
	dΩ.quad
# time: 2023-03-17 02:21:06 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-17 02:21:13 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-17 02:21:19 CET
# mode: julia
	c(pts)
# time: 2023-03-17 02:21:22 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-17 02:52:24 CET
# mode: julia
	pts
# time: 2023-03-17 02:52:33 CET
# mode: julia
	pts[1,1]
# time: 2023-03-17 02:52:38 CET
# mode: julia
	pts[1][1]
# time: 2023-03-17 02:52:41 CET
# mode: julia
	pts[1]
# time: 2023-03-20 00:14:09 CET
# mode: julia
	ε_ps
# time: 2023-03-20 00:14:42 CET
# mode: julia
	C
# time: 2023-03-20 00:26:09 CET
# mode: julia
	u
# time: 2023-03-20 00:26:25 CET
# mode: julia
	ε(u)
# time: 2023-03-20 00:33:31 CET
# mode: julia
	u(x) = x
# time: 2023-03-20 00:33:37 CET
# mode: julia
	a(x) = x
# time: 2023-03-20 00:33:44 CET
# mode: julia
	ε(a)
# time: 2023-03-20 00:33:53 CET
# mode: julia
	eee = ε(a)
# time: 2023-03-20 00:34:04 CET
# mode: julia
	get_array(eee)
# time: 2023-03-20 00:34:43 CET
# mode: julia
	eee = ε(a)
# time: 2023-03-20 00:34:45 CET
# mode: julia
	get_array(eee)
# time: 2023-03-20 00:35:30 CET
# mode: julia
	get_array(eee{2,2,Float64})
# time: 2023-03-20 00:36:21 CET
# mode: julia
	e.data
# time: 2023-03-20 00:36:24 CET
# mode: julia
	eee.data
# time: 2023-03-20 00:59:04 CET
# mode: julia
	u = get_fe_basis(U)
# time: 2023-03-20 00:59:23 CET
# mode: julia
	ε_ps(u)
# time: 2023-03-20 01:03:42 CET
# mode: julia
	Field(x->x+1)
# time: 2023-03-20 01:03:49 CET
# mode: julia
	Gridap.Fields.Field(x->x+1)
# time: 2023-03-20 01:03:56 CET
# mode: julia
	Gridap.Fields.GenericField(x->x+1)
# time: 2023-03-20 01:04:08 CET
# mode: julia
	GenericField(x->x+1)
# time: 2023-03-20 01:04:48 CET
# mode: julia
	ε_ps(u)
# time: 2023-03-20 01:05:05 CET
# mode: julia
	ε_ps∘u
# time: 2023-03-20 01:17:45 CET
# mode: julia
	e
# time: 2023-03-20 01:17:51 CET
# mode: julia
	e = ε(u)
# time: 2023-03-20 01:18:37 CET
# mode: julia
	to_array(e)
# time: 2023-03-20 01:18:48 CET
# mode: julia
	to_array∘e
# time: 2023-03-20 01:19:54 CET
# mode: julia
	e(VectorValue(1.0,0.0))
# time: 2023-03-20 01:20:12 CET
# mode: julia
	e(VectorValue(1.0,0.0))[1]
# time: 2023-03-20 01:20:19 CET
# mode: julia
	res = e(VectorValue(1.0,0.0))[1]
# time: 2023-03-20 01:20:24 CET
# mode: julia
	res[1]
# time: 2023-03-20 01:20:26 CET
# mode: julia
	res[2]
# time: 2023-03-20 01:20:30 CET
# mode: julia
	res[3]
# time: 2023-03-20 01:20:32 CET
# mode: julia
	res[4]
# time: 2023-03-20 01:20:38 CET
# mode: julia
	res[5]
# time: 2023-03-20 01:20:40 CET
# mode: julia
	res[4]
# time: 2023-03-20 01:23:39 CET
# mode: julia
	ε_ps∘u
# time: 2023-03-20 01:24:32 CET
# mode: julia
	L
# time: 2023-03-20 01:24:34 CET
# mode: julia
	ε_ps∘u
# time: 2023-03-20 01:25:10 CET
# mode: julia
	L
# time: 2023-03-20 01:25:12 CET
# mode: julia
	ε_ps∘u
# time: 2023-03-20 01:29:53 CET
# mode: julia
	e
# time: 2023-03-20 01:29:56 CET
# mode: julia
	res
# time: 2023-03-20 01:31:28 CET
# mode: julia
	ε_ps(res)
# time: 2023-03-20 01:31:37 CET
# mode: julia
	res
# time: 2023-03-20 01:32:09 CET
# mode: julia
	L
# time: 2023-03-20 01:32:11 CET
# mode: julia
	res
# time: 2023-03-20 01:32:13 CET
# mode: julia
	ε_ps(res)
# time: 2023-03-20 01:32:38 CET
# mode: julia
	res
# time: 2023-03-20 01:33:53 CET
# mode: julia
	ε_ps(u)
# time: 2023-03-20 01:34:07 CET
# mode: julia
	ε_ps∘u
# time: 2023-03-20 01:37:19 CET
# mode: julia
	sol(x) = sum(x)
# time: 2023-03-20 01:37:27 CET
# mode: julia
	u = interpolate(sol,U)
# time: 2023-03-20 01:37:53 CET
# mode: julia
	sol(x) = x
# time: 2023-03-20 01:37:54 CET
# mode: julia
	u = interpolate(sol,U)
# time: 2023-03-20 01:38:12 CET
# mode: julia
	ε_ps∘u
# time: 2023-03-20 01:40:09 CET
# mode: julia
	ε_ps_op∘ε(u)
# time: 2023-03-20 01:41:02 CET
# mode: julia
	ε_ps(u)
# time: 2023-03-20 01:43:38 CET
# mode: julia
	σ(u)
# time: 2023-03-20 01:44:22 CET
# mode: julia
	G
# time: 2023-03-20 01:44:25 CET
# mode: julia
	S
# time: 2023-03-20 01:44:51 CET
# mode: julia
	S = Gridap.TensorValues.SymTensorValue(S)
# time: 2023-03-20 01:44:55 CET
# mode: julia
	S
# time: 2023-03-20 01:45:02 CET
# mode: julia
	get_array(S)
# time: 2023-03-20 01:45:15 CET
# mode: julia
	inv(S)
# time: 2023-03-20 01:45:42 CET
# mode: julia
	σ(u)
# time: 2023-03-20 01:50:21 CET
# mode: julia
	res
# time: 2023-03-20 01:50:39 CET
# mode: julia
	Gridap.TensorValues.SymTensorValue(zeros(2,2))
# time: 2023-03-20 02:51:59 CET
# mode: julia
	F
# time: 2023-03-20 02:52:37 CET
# mode: julia
	x
# time: 2023-03-20 02:53:29 CET
# mode: julia
	x = VectorValue(1.0, 2.0)
# time: 2023-03-20 02:53:46 CET
# mode: julia
	(F⋅x)
# time: 2023-03-20 02:54:17 CET
# mode: julia
	Ften = TensorValue(F)
# time: 2023-03-20 02:54:23 CET
# mode: julia
	(Ften⋅x)
# time: 2023-03-20 06:53:11 CET
# mode: julia
	Ds
# time: 2023-03-20 06:53:21 CET
# mode: julia
	Dd
# time: 2023-03-20 06:56:10 CET
# mode: julia
	c[1]
# time: 2023-03-20 06:56:21 CET
# mode: julia
	vτ
# time: 2023-03-20 06:56:38 CET
# mode: julia
	c
# time: 2023-03-20 06:56:44 CET
# mode: julia
	quad
# time: 2023-03-20 06:56:53 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-20 06:57:01 CET
# mode: julia
	c(pts)
# time: 2023-03-20 06:57:06 CET
# mode: julia
	c(pts)[1]
# time: 2023-03-20 23:15:05 CET
# mode: julia
	arrDd
# time: 2023-03-20 23:15:14 CET
# mode: julia
	arrDs
# time: 2023-03-21 01:34:10 CET
# mode: julia
	ε_p
# time: 2023-03-21 01:37:24 CET
# mode: julia
	u = get_trial_fe_basis(U)
# time: 2023-03-21 01:37:39 CET
# mode: julia
	u = get_trial_fe_basis(Uv)
# time: 2023-03-21 01:37:50 CET
# mode: julia
	w = get_trial_fe_basis(Uw)
# time: 2023-03-21 01:38:07 CET
# mode: julia
	t = get_trial_fe_basis(Ut)
# time: 2023-03-21 01:38:44 CET
# mode: julia
	ε_p(u,t)
# time: 2023-03-21 01:39:24 CET
# mode: julia
	ε_p(v,t) ⊙ (P(u,θ))
# time: 2023-03-21 01:40:02 CET
# mode: julia
	u = get_trial_fe_basis(Uv)
# time: 2023-03-21 01:40:19 CET
# mode: julia
	v = get_fe_basis(Vv)
# time: 2023-03-21 01:40:39 CET
# mode: julia
	θ = get_trial_fe_basis(Ut)
# time: 2023-03-21 01:40:52 CET
# mode: julia
	t = get_fe_basis(Vt)
# time: 2023-03-21 01:40:58 CET
# mode: julia
	ε_p(v,t) ⊙ (P(u,θ))
# time: 2023-03-21 01:41:45 CET
# mode: julia
	u(x) = 3*x
# time: 2023-03-21 01:41:51 CET
# mode: julia
	fu(x) = 3*x
# time: 2023-03-21 01:41:59 CET
# mode: julia
	ft(x) = 3
# time: 2023-03-21 01:42:14 CET
# mode: julia
	ε_p_op(fu,ft)
# time: 2023-03-21 01:42:33 CET
# mode: julia
	ε_p(fu,ft)
# time: 2023-03-21 02:04:24 CET
# mode: julia
	Float64(3)
# time: 2023-03-21 03:02:32 CET
# mode: julia
	vτ.cell_basis
# time: 2023-03-21 03:03:07 CET
# mode: julia
	vτ.single_fields
# time: 2023-03-21 03:03:19 CET
# mode: julia
	vτ.single_fields[1]
# time: 2023-03-21 03:04:57 CET
# mode: julia
	contr.dict
# time: 2023-03-21 03:05:15 CET
# mode: julia
	contr
# time: 2023-03-21 03:05:18 CET
# mode: julia
	contr.dict
# time: 2023-03-21 03:05:27 CET
# mode: julia
	first(contr.dict)
# time: 2023-03-21 03:05:36 CET
# mode: julia
	first(contr.dict)[2]
# time: 2023-03-21 03:05:58 CET
# mode: julia
	first(contr.dict)[2][1]
# time: 2023-03-21 03:07:57 CET
# mode: julia
	first(contr.dict)[2][1][1,1]
# time: 2023-03-21 03:10:29 CET
# mode: julia
	reffe1D
# time: 2023-03-21 03:10:59 CET
# mode: julia
	reffe = ReferenceFE(QUAD,lagrangian,Float64,2)
# time: 2023-03-21 03:11:28 CET
# mode: julia
	basis = get_basis(reffe)
# time: 2023-03-21 03:11:37 CET
# mode: julia
	basis = Gridap.ReferenceFEs.get_basis(reffe)
# time: 2023-03-21 03:11:42 CET
# mode: julia
	basis = Gridap.ReferenceFEs.get_fe_basis(reffe)
# time: 2023-03-21 03:11:46 CET
# mode: julia
	basis = Gridap.ReferenceFEs.get_dof_basis(reffe)
# time: 2023-03-21 03:11:57 CET
# mode: julia
	methods(reffe)
# time: 2023-03-21 03:12:19 CET
# mode: julia
	basis.nodes
# time: 2023-03-21 03:13:47 CET
# mode: julia
	Vv.fe_dof_basis
# time: 2023-03-21 03:13:52 CET
# mode: julia
	Vv.fe_dof_basis.cell_dof
# time: 2023-03-21 03:13:56 CET
# mode: julia
	Vv.fe_dof_basis.cell_dof[1]
# time: 2023-03-21 03:13:59 CET
# mode: julia
	Vv.fe_dof_basis.cell_dof[1].nodes
# time: 2023-03-21 03:14:29 CET
# mode: julia
	Vv.fe_dof_basis.cell_dof[1].dof_to_comp
# time: 2023-03-21 05:19:16 CET
# mode: julia
	first(contr.dict)
# time: 2023-03-21 05:19:26 CET
# mode: julia
	first(contr.dict)[1]
# time: 2023-03-21 05:19:29 CET
# mode: julia
	first(contr.dict)[2]
# time: 2023-03-21 05:19:35 CET
# mode: julia
	first(contr.dict)
# time: 2023-03-21 05:19:55 CET
# mode: julia
	first(contr.dict)[2][1]
# time: 2023-03-21 05:20:11 CET
# mode: julia
	first(contr.dict)[2][1][1,1]
# time: 2023-03-21 05:20:54 CET
# mode: julia
	first(contr.dict)[2][1]
# time: 2023-03-21 05:25:02 CET
# mode: julia
	first(contr.dict)[2][1][1,1]
# time: 2023-03-21 05:31:35 CET
# mode: julia
	first(contr.dict)[2][1][1,1][1]
# time: 2023-03-21 05:31:43 CET
# mode: julia
	first(contr.dict)[2][1][1,1]
# time: 2023-03-21 05:31:46 CET
# mode: julia
	first(contr.dict)[2][1][1,2]
# time: 2023-03-21 05:31:52 CET
# mode: julia
	first(contr.dict)[2][1][2,2]
# time: 2023-03-21 05:40:56 CET
# mode: julia
	uh
# time: 2023-03-21 05:43:09 CET
# mode: julia
	first(contr.dict)[2]
# time: 2023-03-21 05:43:20 CET
# mode: julia
	first(contr.dict)[2][1]
# time: 2023-03-21 05:43:30 CET
# mode: julia
	first(contr.dict)[2][1][1,3]
# time: 2023-03-21 05:44:19 CET
# mode: julia
	Vv
# time: 2023-03-21 05:44:26 CET
# mode: julia
	reffe1D
# time: 2023-03-21 05:44:52 CET
# mode: julia
	dΩ
# time: 2023-03-21 05:44:57 CET
# mode: julia
	dΩ.quad
# time: 2023-03-21 05:45:08 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-21 05:45:10 CET
# mode: julia
	quad
# time: 2023-03-21 05:45:32 CET
# mode: julia
	quad.cell_quad
# time: 2023-03-21 05:45:40 CET
# mode: julia
	quad.cell_quad.axes
# time: 2023-03-21 05:45:45 CET
# mode: julia
	quad.cell_quad.axes[1]
# time: 2023-03-21 05:45:49 CET
# mode: julia
	quad.cell_quad.axes[1][1]
# time: 2023-03-21 05:45:52 CET
# mode: julia
	quad.cell_quad.axes[1][2]
# time: 2023-03-21 05:45:54 CET
# mode: julia
	quad.cell_quad.axes[1][3]
# time: 2023-03-21 05:46:08 CET
# mode: julia
	quad.cell_quad.axes[1][100]
# time: 2023-03-21 05:46:11 CET
# mode: julia
	quad.cell_quad.axes
# time: 2023-03-21 05:46:23 CET
# mode: julia
	OneTo.(10)
# time: 2023-03-21 05:46:30 CET
# mode: julia
	Base.OneTo.(10)
# time: 2023-03-21 05:46:38 CET
# mode: julia
	print(Base.OneTo.(10))
# time: 2023-03-21 05:46:45 CET
# mode: julia
	println(Base.OneTo.(10))
# time: 2023-03-21 05:46:57 CET
# mode: julia
	println(Base.OneTo(10))
# time: 2023-03-21 05:47:06 CET
# mode: julia
	Base.OneTo(10)[1]
# time: 2023-03-21 05:47:08 CET
# mode: julia
	Base.OneTo(10)[2]
# time: 2023-03-21 05:52:37 CET
# mode: julia
	Ω
# time: 2023-03-21 05:52:45 CET
# mode: julia
	Ω.grid
# time: 2023-03-21 05:52:58 CET
# mode: julia
	Ω.grid.CartesianGrid()
# time: 2023-03-21 05:53:14 CET
# mode: julia
	Ω.grid.node_coords
# time: 2023-03-21 05:53:19 CET
# mode: julia
	Ω.grid.node_coords[1]
# time: 2023-03-21 05:53:23 CET
# mode: julia
	Ω.grid.node_coords[2]
# time: 2023-03-21 05:53:48 CET
# mode: julia
	Ω.grid.cell_node_ids
# time: 2023-03-21 05:54:16 CET
# mode: julia
	Ω.grid.cell_type
# time: 2023-03-22 00:08:16 CET
# mode: julia
	uuu = get_trial_fe_basis(U)
# time: 2023-03-22 00:08:46 CET
# mode: julia
	eee = ε0(u)
# time: 2023-03-22 00:08:50 CET
# mode: julia
	eee = ε0(uuu)
# time: 2023-03-22 00:09:08 CET
# mode: julia
	eeetr = tr(ε0(uuu))
# time: 2023-03-22 02:29:05 CET
# mode: julia
	aux = σ_bimat∘(ε0,tags)
# time: 2023-03-22 02:40:38 CET
# mode: julia
	∫(1.0)dΓ
# time: 2023-03-22 02:40:50 CET
# mode: julia
	aaa = ∫(1.0)dΓ
# time: 2023-03-22 02:40:52 CET
# mode: julia
	aaa
# time: 2023-03-22 06:42:34 CET
# mode: julia
	tags
# time: 2023-03-22 06:44:25 CET
# mode: julia
	x = [1,2,3]
# time: 2023-03-22 06:44:30 CET
# mode: julia
	x+1
# time: 2023-03-22 06:44:35 CET
# mode: julia
	x.+1
# time: 2023-03-22 06:44:53 CET
# mode: julia
	E_layer(tags)
# time: 2023-03-22 06:45:01 CET
# mode: julia
	E_layer.(tags)
# time: 2023-03-22 06:48:14 CET
# mode: julia
	pts = get_cell_coords(Ω)
# time: 2023-03-22 06:48:21 CET
# mode: julia
	pts = get_cell_coordinatess(Ω)
# time: 2023-03-22 06:48:23 CET
# mode: julia
	pts = get_cell_coordinates(Ω)
# time: 2023-03-22 06:48:34 CET
# mode: julia
	z_cf(pts)
# time: 2023-03-22 06:49:51 CET
# mode: julia
	z_cf
# time: 2023-03-22 06:49:58 CET
# mode: julia
	z_cf.cell_field
# time: 2023-03-22 06:50:10 CET
# mode: julia
	z_cf.cell_field[1]
# time: 2023-03-22 06:50:25 CET
# mode: julia
	z_cf.cell_field[1](VectorValue(1.0,2.0))
# time: 2023-03-22 06:51:35 CET
# mode: julia
	lazy_map(evaluate,z_cf.cell_field,)
# time: 2023-03-22 06:51:39 CET
# mode: julia
	pts.args
# time: 2023-03-22 06:51:49 CET
# mode: julia
	get_data(pts)
# time: 2023-03-22 06:51:59 CET
# mode: julia
	Gridap.CellData.get_data(pts)
# time: 2023-03-22 06:52:13 CET
# mode: julia
	pts
# time: 2023-03-22 06:52:23 CET
# mode: julia
	lazy_map(evaluate,z_cf.cell_field,pts)
# time: 2023-03-22 06:53:31 CET
# mode: julia
	cpts = CellPoint(pts)
# time: 2023-03-22 06:53:59 CET
# mode: julia
	Gridap.CellData.CellPoint(pts,Ω)
# time: 2023-03-22 07:49:10 CET
# mode: julia
	4/6
# time: 2023-03-22 22:47:59 CET
# mode: julia
	E_field
# time: 2023-03-22 22:48:26 CET
# mode: julia
	Da(E) = ∫( E )*dΓ
# time: 2023-03-22 22:48:38 CET
# mode: julia
	Da∘E_field
# time: 2023-03-22 22:48:57 CET
# mode: julia
	Da = ∫( E )*dΓ
# time: 2023-03-22 22:49:07 CET
# mode: julia
	Da = ∫( E_field )*dΓ
# time: 2023-03-22 22:49:15 CET
# mode: julia
	Da2 = ∫( E_field )*dΓ
# time: 2023-03-22 22:49:33 CET
# mode: julia
	Da2
# time: 2023-03-22 22:50:01 CET
# mode: julia
	quad = dΩ.quad
# time: 2023-03-22 22:50:13 CET
# mode: julia
	pts = get_cell_points(quad)
# time: 2023-03-22 22:50:29 CET
# mode: julia
	Da2(pts)
# time: 2023-03-22 22:50:39 CET
# mode: julia
	Da(pts)
# time: 2023-03-22 22:51:24 CET
# mode: julia
	Da2 = ∫( E_field )*dΓ
# time: 2023-03-22 22:51:34 CET
# mode: julia
	Da2(E_field) = ∫( E_field )*dΓ
# time: 2023-03-22 22:51:38 CET
# mode: julia
	Da3(E_field) = ∫( E_field )*dΓ
# time: 2023-03-23 00:54:52 CET
# mode: julia
	Da(E_field)
# time: 2023-03-23 00:55:03 CET
# mode: julia
	contr = Da(E_field)
# time: 2023-03-23 00:55:15 CET
# mode: julia
	first(contr.dict)
# time: 2023-03-23 00:55:35 CET
# mode: julia
	first(contr.dict)[2]
# time: 2023-03-23 00:56:03 CET
# mode: julia
	Da(E_field)
# time: 2023-03-27 03:04:06 CEST
# mode: julia
	step(3,5)
# time: 2023-03-27 03:04:11 CEST
# mode: julia
	step(8,5)
# time: 2023-03-27 03:10:15 CEST
# mode: julia
	da_fun(3.0)
# time: 2023-03-27 03:19:36 CEST
# mode: julia
	Da_fun(E_field)
# time: 2023-03-27 03:19:44 CEST
# mode: julia
	Dd_fun(E_field,z_cf)
# time: 2023-03-27 03:20:17 CEST
# mode: julia
	z_cf
# time: 2023-03-27 03:29:06 CEST
# mode: julia
	E_layer
# time: 2023-03-27 03:29:13 CEST
# mode: julia
	E_layer.1
# time: 2023-03-27 03:30:08 CEST
# mode: julia
	E_layer∘z_cf
# time: 2023-03-27 03:30:23 CEST
# mode: julia
	E_layer∘z_cf[1]
# time: 2023-03-27 03:30:34 CEST
# mode: julia
	aaa = E_layer∘z_cf
# time: 2023-03-27 03:30:35 CEST
# mode: julia
	aaa
# time: 2023-03-27 03:30:42 CEST
# mode: julia
	first(aaa)
# time: 2023-03-27 03:31:34 CEST
# mode: julia
	aaa
# time: 2023-03-27 05:20:30 CEST
# mode: julia
	dd(3.0)
# time: 2023-03-27 05:20:33 CEST
# mode: julia
	daa(3.0)
# time: 2023-03-27 06:23:51 CEST
# mode: julia
	line quad = dΩ.quad
# time: 2023-03-27 06:24:03 CEST
# mode: julia
	line = dΓ.quad
# time: 2023-03-27 06:24:24 CEST
# mode: julia
	pts = get_cell_points(line)
# time: 2023-03-27 06:25:31 CEST
# mode: julia
	db∘z_cf
# time: 2023-03-27 06:26:12 CEST
# mode: julia
	db_fun(E_field,z_cf,z_val)
# time: 2023-03-27 06:46:59 CEST
# mode: julia
	aux
# time: 2023-03-27 06:47:14 CEST
# mode: julia
	aux[10]
# time: 2023-03-27 06:47:18 CEST
# mode: julia
	aux[9]
# time: 2023-03-27 06:47:22 CEST
# mode: julia
	aux[8]
# time: 2023-03-29 04:13:29 CEST
# mode: julia
	3.0 ⋅ VectorValue(1.0)
# time: 2023-03-29 06:02:19 CEST
# mode: julia
	ca
# time: 2023-03-29 06:02:29 CEST
# mode: julia
	ca(pts)
# time: 2023-03-29 06:02:36 CEST
# mode: julia
	ca(pts)[2]
# time: 2023-03-29 06:02:39 CEST
# mode: julia
	ca(pts)[1]
# time: 2023-03-29 06:03:24 CEST
# mode: julia
	cc1(pts)[1]
# time: 2023-03-29 06:15:49 CEST
# mode: julia
	first(contr.dict)
# time: 2023-03-29 06:15:53 CEST
# mode: julia
	first(contr.dict)[1]
# time: 2023-03-29 06:16:04 CEST
# mode: julia
	first(contr.dict)[1][1]
# time: 2023-03-29 06:16:10 CEST
# mode: julia
	first(contr.dict)[2][1]
# time: 2023-03-29 06:16:57 CEST
# mode: julia
	first(contr.dict)[2]
# time: 2023-03-29 06:17:03 CEST
# mode: julia
	first(contr.dict)[2][1][1,1]
# time: 2023-03-29 06:17:07 CEST
# mode: julia
	first(contr.dict)[2][1][3,1]
# time: 2023-03-29 06:18:15 CEST
# mode: julia
	first(contr.dict)[2][1]
# time: 2023-03-29 06:18:31 CEST
# mode: julia
	first(contr.dict)[2][1][3,1]
# time: 2023-03-29 06:18:36 CEST
# mode: julia
	first(contr.dict)[2][1][2,2]
# time: 2023-03-29 06:47:19 CEST
# mode: julia
	AAA = ε(u)
# time: 2023-03-29 06:47:22 CEST
# mode: julia
	AAA
# time: 2023-03-29 06:47:40 CEST
# mode: julia
	ε∘u
# time: 2023-03-29 06:48:06 CEST
# mode: julia
	ddd = ε∘u
# time: 2023-03-29 06:48:08 CEST
# mode: julia
	ddd
# time: 2023-03-29 06:53:39 CEST
# mode: julia
	ε∘u
# time: 2023-03-29 06:53:45 CEST
# mode: julia
	ddd
# time: 2023-03-30 00:02:15 CEST
# mode: julia
	computeCT
# time: 2023-03-30 00:10:25 CEST
# mode: julia
	u(x) = 1.0
# time: 2023-03-30 00:10:47 CEST
# mode: julia
	var_a = ε(u)
# time: 2023-03-30 00:11:11 CEST
# mode: julia
	var_a.f
# time: 2023-03-30 00:11:17 CEST
# mode: julia
	var_a.data
# time: 2023-03-30 02:01:35 CEST
# mode: julia
	aaaa = ε_data_op(ε)
# time: 2023-03-30 02:01:49 CEST
# mode: julia
	\varepsilon
# time: 2023-03-30 02:01:52 CEST
# mode: julia
	ε
# time: 2023-03-30 02:02:18 CEST
# mode: julia
	rr = ε(u)
# time: 2023-03-30 02:02:23 CEST
# mode: julia
	rr
# time: 2023-03-30 02:02:29 CEST
# mode: julia
	rr.f
# time: 2023-03-30 02:18:45 CEST
# mode: julia
	import Pkg; Pkg.add("PkgTemplates")
# time: 2023-03-30 02:19:08 CEST
# mode: julia
	using PkgTemplates
# time: 2023-03-30 02:21:16 CEST
# mode: julia
	t = Template(;
	 user="gridap",
	 authors=["Francesc Turon <fturon@cimne.upc.edu>"],
	 dir=pwd(),
	 julia=v"1.8.3",
	 plugins=[
	     License(; name="MIT", path=nothing, destination="LICENSE.md"),
	     CompatHelper(),
	     Codecov(),    
	     GitHubActions(;
	     osx=false,
	     windows=false,
	     ),
	     Documenter{GitHubActions}(),
	     Git(;
	         ignore=["*.jl.*.cov",
	                 "*.jl.cov",
	                 "*.jl.mem",
	                 "*.code-workspace",
	                 ".DS_Store",
	                 "docs/build/",
	                 "Manifest.toml",
	                 "tmp/"],
	         ssh=true
	     ),
	 ],
	)
# time: 2023-03-30 02:21:53 CEST
# mode: julia
	t("CAMMS")
# time: 2023-03-30 03:02:47 CEST
# mode: julia
	paraview
# time: 2023-03-30 03:02:52 CEST
# mode: julia
	exit()
# time: 2023-03-30 03:39:33 CEST
# mode: julia
	x = [1,2,3]
# time: 2023-03-30 03:39:44 CEST
# mode: julia
	x = (1,2,3)
# time: 2023-03-30 03:40:04 CEST
# mode: julia
	y = (4,5,6)
# time: 2023-03-30 03:40:09 CEST
# mode: julia
	VectorValue
# time: 2023-03-30 03:40:13 CEST
# mode: julia
	VectorValue(x,y)
# time: 2023-03-30 03:40:21 CEST
# mode: julia
	VectorValue(x...,y...)
# time: 2023-03-30 03:40:58 CEST
# mode: julia
	x = SymTensorValue(2,3,4)
# time: 2023-03-30 03:41:06 CEST
# mode: julia
	x.data
# time: 2023-03-30 03:41:18 CEST
# mode: julia
	VectorValue(x.data...)
# time: 2023-03-30 03:41:23 CEST
# mode: julia
	VectorValue(x.data...,x.data...)
# time: 2023-03-30 03:43:15 CEST
# mode: julia
	x = zero(U)
# time: 2023-03-30 03:43:51 CEST
# mode: julia
	x
# time: 2023-03-30 03:44:07 CEST
# mode: julia
	y = ε(x)
# time: 2023-03-30 03:44:33 CEST
# mode: julia
	εₘ_op∘(ct,y,y,y)
# time: 2023-03-30 03:44:40 CEST
# mode: julia
	modModel.εₘ_op∘(ct,y,y,y)
# time: 2023-03-30 03:44:51 CEST
# mode: julia
	z = modModel.εₘ_op∘(ct,y,y,y)
# time: 2023-03-30 03:44:54 CEST
# mode: julia
	pts
# time: 2023-03-30 03:45:14 CEST
# mode: julia
	pts = get_cell_points(dΩ.quad)
# time: 2023-03-30 03:45:18 CEST
# mode: julia
	z(pts)
# time: 2023-03-30 03:45:38 CEST
# mode: julia
	pts.cell_ref_point
# time: 2023-03-30 03:46:02 CEST
# mode: julia
	ct
# time: 2023-03-30 03:46:18 CEST
# mode: julia
	ctcf = CellField(ct,Ω)
# time: 2023-03-30 03:46:28 CEST
# mode: julia
	z = modModel.εₘ_op∘(ctcf,y,y,y)
# time: 2023-03-30 03:46:44 CEST
# mode: julia
	z(pts)
# time: 2023-03-30 03:46:58 CEST
# mode: julia
	zz = z(pts)
# time: 2023-03-30 03:47:06 CEST
# mode: julia
	zz[1]
# time: 2023-03-30 03:47:21 CEST
# mode: julia
	x[1]
# time: 2023-03-30 03:52:49 CEST
# mode: julia
	using Gridap.Arrays
# time: 2023-03-30 03:53:05 CEST
# mode: julia
	vals = [1.0,0.0]
# time: 2023-03-30 03:53:57 CEST
# mode: julia
	ptrs = [1,1,1,1,1,1,1,2,2,2,2,2,2,2,2]
# time: 2023-03-30 03:54:18 CEST
# mode: julia
	x = CompressedArray(vals,ptrs)
# time: 2023-03-30 03:56:33 CEST
# mode: julia
	Ω
# time: 2023-03-30 03:56:40 CEST
# mode: julia
	num_cells(Ω)
# time: 2023-03-30 03:56:55 CEST
# mode: julia
	cts = Fill(1.0,100)
# time: 2023-03-30 03:57:05 CEST
# mode: julia
	cts = fill(1.0,100)
# time: 2023-03-30 03:57:12 CEST
# mode: julia
	CellField(cts,Ω)
# time: 2023-03-30 03:57:35 CEST
# mode: julia
	get_model(Ω)
# time: 2023-03-30 03:57:42 CEST
# mode: julia
	get_background_model(Ω)
# time: 2023-03-30 03:58:12 CEST
# mode: julia
	labels = get_face_labeling(Ω)
# time: 2023-03-30 03:58:21 CEST
# mode: julia
	labels = get_face_labeling(model)
# time: 2023-03-30 03:58:38 CEST
# mode: julia
	labels.tag_to_entities
# time: 2023-03-30 05:19:54 CEST
# mode: julia
	JK = SymTensorValue( (11) )
# time: 2023-03-30 06:13:24 CEST
# mode: julia
	tags
# time: 2023-03-30 06:16:07 CEST
# mode: julia
	CT_CF
# time: 2023-03-30 06:16:12 CEST
# mode: julia
	u
# time: 2023-03-30 06:19:39 CEST
# mode: julia
	CT_CF∘tags
# time: 2023-03-30 06:27:49 CEST
# mode: julia
	CT_CF
# time: 2023-03-30 06:28:05 CEST
# mode: julia
	pts = \Omega
# time: 2023-03-30 06:28:11 CEST
# mode: julia
	pts = Ω.quad
# time: 2023-03-30 06:28:41 CEST
# mode: julia
	pts = get_cell_points(dΩ.quad)
# time: 2023-03-30 06:28:46 CEST
# mode: julia
	z(pts)
# time: 2023-03-30 06:28:55 CEST
# mode: julia
	CT_CF(pts)
# time: 2023-03-30 06:29:13 CEST
# mode: julia
	CT_CF(pts)[1]
# time: 2023-03-30 06:29:44 CEST
# mode: julia
	u = 3.0
# time: 2023-03-30 06:29:51 CEST
# mode: julia
	ε()
# time: 2023-03-30 06:29:54 CEST
# mode: julia
	ε(u)
# time: 2023-03-30 06:30:51 CEST
# mode: julia
	CT_CF(pts)[4]
# time: 2023-03-30 06:31:33 CEST
# mode: julia
	CT_CF∘(tags)
# time: 2023-03-30 06:31:49 CEST
# mode: julia
	CT_CF
# time: 2023-03-30 06:31:53 CEST
# mode: julia
	CT_CF∘(tags)
# time: 2023-03-30 06:32:07 CEST
# mode: julia
	cts
# time: 2023-03-30 06:33:20 CEST
# mode: julia
	tags
# time: 2023-03-30 06:34:20 CEST
# mode: julia
	tags[1]
# time: 2023-03-30 06:34:24 CEST
# mode: julia
	tags[1] = 4
# time: 2023-03-30 06:34:26 CEST
# mode: julia
	tags
# time: 2023-03-30 06:43:21 CEST
# mode: julia
	replace(tags, 4 => 7)
# time: 2023-03-30 06:51:56 CEST
# mode: julia
	tags
# time: 2023-03-30 06:52:26 CEST
# mode: julia
	iflag
# time: 2023-03-30 06:52:34 CEST
# mode: julia
	matFlag
# time: 2023-03-30 06:52:49 CEST
# mode: julia
	enumerate(matFlag)
# time: 2023-03-30 06:55:51 CEST
# mode: julia
	replace(tags, 7 => 3)
# time: 2023-03-30 06:55:58 CEST
# mode: julia
	tags
# time: 2023-03-30 07:37:56 CEST
# mode: julia
	ct1
# time: 2023-03-30 07:38:01 CEST
# mode: julia
	ct1[1]
# time: 2023-03-30 07:38:07 CEST
# mode: julia
	ct1[2]
# time: 2023-03-30 07:38:09 CEST
# mode: julia
	ct1[3]
# time: 2023-03-30 07:38:11 CEST
# mode: julia
	ct1[4]
# time: 2023-03-30 07:38:14 CEST
# mode: julia
	ct1[5]
# time: 2023-03-30 07:38:16 CEST
# mode: julia
	ct1[6]
# time: 2023-03-30 07:38:18 CEST
# mode: julia
	ct1[7]
# time: 2023-03-30 07:38:20 CEST
# mode: julia
	ct1[8]
# time: 2023-03-30 07:40:08 CEST
# mode: julia
	ct1[7]
# time: 2023-03-30 07:40:11 CEST
# mode: julia
	ct1[6]
# time: 2023-03-30 07:40:18 CEST
# mode: julia
	ct1
# time: 2023-03-30 07:41:52 CEST
# mode: julia
	SymTensorValue{3,Int}(11,21.0,22,11,21.0,22)
# time: 2023-03-30 07:44:52 CEST
# mode: julia
	ct1
# time: 2023-03-31 07:22:52 CEST
# mode: julia
	pwd()*"/models/segonModel.msh"
# time: 2023-03-31 07:25:53 CEST
# mode: julia
	pwd()
# time: 2023-03-31 08:18:22 CEST
# mode: julia
	tagsPS
# time: 2023-03-31 08:20:45 CEST
# mode: julia
	typeof(tagsPS)
# time: 2023-03-31 08:21:21 CEST
# mode: julia
	tagsPS2 = Vector{Int8}
# time: 2023-03-31 08:21:26 CEST
# mode: julia
	tagsPS2
# time: 2023-03-31 08:24:11 CEST
# mode: julia
	tagsPS
# time: 2023-03-31 08:26:18 CEST
# mode: julia
	replace!
# time: 2023-04-01 04:29:12 CEST
# mode: julia
	b*(1/1)*(z2^1-z1^1)*G
# time: 2023-04-01 04:29:22 CEST
# mode: julia
	b*(1/3)*(z2^3-z1^3)*E11
# time: 2023-04-01 04:38:07 CEST
# mode: julia
	G
# time: 2023-04-01 04:38:25 CEST
# mode: julia
	50*100*G
# time: 2023-04-01 05:09:20 CEST
# mode: julia
	G
# time: 2023-04-03 07:25:51 CEST
# mode: julia
	LOAD_PATH
# time: 2023-04-03 08:11:40 CEST
# mode: julia
	MT
# time: 2023-04-03 08:12:15 CEST
# mode: julia
	CTf
# time: 2023-04-03 08:42:42 CEST
# mode: julia
	CTf[1]
# time: 2023-04-03 08:51:24 CEST
# mode: julia
	mt
# time: 2023-04-03 08:51:27 CEST
# mode: julia
	MT
# time: 2023-04-03 08:52:24 CEST
# mode: julia
	SymTensorValue((3.0))
# time: 2023-04-03 08:52:42 CEST
# mode: julia
	clear
# time: 2023-04-04 00:49:29 CEST
# mode: julia
	SymTensorValue( (4.0) )
# time: 2023-04-04 00:49:49 CEST
# mode: julia
	AAA = SymTensorValue((3.0))
# time: 2023-04-04 00:49:53 CEST
# mode: julia
	AAA*3.0
# time: 2023-04-04 01:32:57 CEST
# mode: julia
	using Gridap.TensorValues
# time: 2023-04-04 01:33:14 CEST
# mode: julia
	SymTensorValue(1.0,2.0,3.0)
# time: 2023-04-04 01:33:20 CEST
# mode: julia
	x = SymTensorValue(1.0,2.0,3.0)
# time: 2023-04-04 01:33:21 CEST
# mode: julia
	x
# time: 2023-04-04 01:33:28 CEST
# mode: julia
	Matrix(x)
# time: 2023-04-04 01:33:43 CEST
# mode: julia
	get_array(x)
# time: 2023-04-04 01:34:28 CEST
# mode: julia
	x = SymTensorValue(1.0,2.0,3.0,4.0,5.0,6.0)
# time: 2023-04-04 01:34:30 CEST
# mode: julia
	get_array(x)
# time: 2023-04-04 01:38:13 CEST
# mode: julia
	x = SymTensorValue(1.0)
# time: 2023-04-04 01:54:18 CEST
# mode: julia
	Γ = Boundary(model)
# time: 2023-04-04 01:54:37 CEST
# mode: julia
	num_cells(Γ)
# time: 2023-04-04 01:54:54 CEST
# mode: julia
	get_cell_maps(Γ)
# time: 2023-04-04 01:54:59 CEST
# mode: julia
	cell_maps(Γ)
# time: 2023-04-04 01:55:07 CEST
# mode: julia
	using Gridap.Geometry
# time: 2023-04-04 01:55:08 CEST
# mode: julia
	cell_maps(Γ)
# time: 2023-04-04 01:55:11 CEST
# mode: julia
	get_cell_maps(Γ)
# time: 2023-04-04 01:55:26 CEST
# mode: julia
	get_cell_maps(dΩ)
# time: 2023-04-04 01:55:42 CEST
# mode: julia
	get_cell_maps(Ω)
# time: 2023-04-04 01:55:48 CEST
# mode: julia
	get_cell_maps(model)
# time: 2023-04-04 01:56:33 CEST
# mode: julia
	get_cell_map(Γ)
# time: 2023-04-04 01:56:41 CEST
# mode: julia
	cmaps = get_cell_map(Γ)
# time: 2023-04-04 01:56:54 CEST
# mode: julia
	m = cmaps[1]
# time: 2023-04-04 01:57:22 CEST
# mode: julia
	m(VectorValue(1.0))
# time: 2023-04-04 01:57:30 CEST
# mode: julia
	m(VectorValue(0.0))
# time: 2023-04-04 01:57:42 CEST
# mode: julia
	m(VectorValue(0.0,1.0))
# time: 2023-04-04 01:57:46 CEST
# mode: julia
	m(VectorValue(0.0))
# time: 2023-04-04 01:58:46 CEST
# mode: julia
	model = CartesianDiscreteModel((0,1,0,1),(2,2))
# time: 2023-04-04 01:59:14 CEST
# mode: julia
	Γ = Boundary(model)
# time: 2023-04-04 01:59:24 CEST
# mode: julia
	Λ = Skeleton(model)
# time: 2023-04-04 01:59:36 CEST
# mode: julia
	Ω = Triangulation(model)
# time: 2023-04-04 02:00:46 CEST
# mode: julia
	cmaps = get_cell_map(Γ)
# time: 2023-04-04 02:01:04 CEST
# mode: julia
	m = cmaps[1]
# time: 2023-04-04 02:01:17 CEST
# mode: julia
	m(VectorValue(0.0))
# time: 2023-04-04 02:01:23 CEST
# mode: julia
	m(VectorValue(1.0))
# time: 2023-04-04 02:01:31 CEST
# mode: julia
	m(VectorValue(0.5))
# time: 2023-04-04 02:01:49 CEST
# mode: julia
	using Gridap.Arrays
# time: 2023-04-04 02:02:01 CEST
# mode: julia
	lazy_map(inverse_map,cmaps)
# time: 2023-04-04 02:04:47 CEST
# mode: julia
	dΓ = Measure(Γ,2)
# time: 2023-04-04 02:05:06 CEST
# mode: julia
	get_cell_point(dΓ.quad)
# time: 2023-04-04 02:05:09 CEST
# mode: julia
	get_cell_points(dΓ.quad)
# time: 2023-04-04 02:05:14 CEST
# mode: julia
	pts = get_cell_points(dΓ.quad)
# time: 2023-04-04 02:05:32 CEST
# mode: julia
	pts.cell_phys_point
# time: 2023-04-04 02:05:59 CEST
# mode: julia
	pts.cell_ref_point
# time: 2023-04-04 02:09:08 CEST
# mode: julia
	u = get_trial_basis(U)
# time: 2023-04-04 02:09:12 CEST
# mode: julia
	u = get_trial_fe_basis(U)
# time: 2023-04-04 02:09:49 CEST
# mode: julia
	CellData.get_data(u)
# time: 2023-04-04 02:09:58 CEST
# mode: julia
	using Gridap.CellData
# time: 2023-04-04 02:10:00 CEST
# mode: julia
	CellData.get_data(u)
# time: 2023-04-04 02:10:20 CEST
# mode: julia
	evaluate(CellData.get_data(u),pts.cell_ref_point)
# time: 2023-04-04 02:10:44 CEST
# mode: julia
	lazy_map(evaluate,CellData.get_data(u),pts.cell_ref_point)
# time: 2023-04-04 02:11:07 CEST
# mode: julia
	u.basis_style
# time: 2023-04-04 02:11:13 CEST
# mode: julia
	u.domain_style
# time: 2023-04-04 02:11:26 CEST
# mode: julia
	u.cell_basis
# time: 2023-04-04 02:11:29 CEST
# mode: julia
	u.cell_basis[1]
# time: 2023-04-04 02:11:33 CEST
# mode: julia
	u.cell_basis[1][1]
# time: 2023-04-04 02:11:40 CEST
# mode: julia
	u.cell_basis[1][1](VectorValue(1.0))
# time: 2023-04-04 02:11:45 CEST
# mode: julia
	u.cell_basis[1][1](VectorValue(1.0,0.0))
# time: 2023-04-04 02:21:12 CEST
# mode: julia
	get_tangent(Γ)
# time: 2023-04-04 02:21:15 CEST
# mode: julia
	get_tangents(Γ)
# time: 2023-04-04 02:21:20 CEST
# mode: julia
	get_tangent_vectors(Γ)
# time: 2023-04-04 02:21:21 CEST
# mode: julia
	get_tangent_vector(Γ)
# time: 2023-04-04 02:23:23 CEST
# mode: julia
	get_normal_vector(Γ)
# time: 2023-04-04 02:23:48 CEST
# mode: julia
	e2_cf = get_normal_vector(Γ)
# time: 2023-04-04 02:24:03 CEST
# mode: julia
	e1_cf = get_tangent_vector(Γ)
# time: 2023-04-04 02:29:25 CEST
# mode: julia
	cmaps
# time: 2023-04-04 02:30:04 CEST
# mode: julia
	m = cmaps(1)
# time: 2023-04-04 02:30:09 CEST
# mode: julia
	m = cmaps[1]
# time: 2023-04-04 02:30:36 CEST
# mode: julia
	m(VectorValue(1.0)) - m(VectorValue(0.0))
# time: 2023-04-04 02:31:07 CEST
# mode: julia
	lazy_map(m->m(VectorValue(1.0)) - m(VectorValue(0.0)),cmaps)
# time: 2023-04-04 02:32:08 CEST
# mode: julia
	lazy_map(t->t/norm(t),lazy_map(m->m(VectorValue(1.0)) - m(VectorValue(0.0)),cmaps))
# time: 2023-04-04 02:32:29 CEST
# mode: julia
	CellField(lazy_map(t->t/norm(t),lazy_map(m->m(VectorValue(1.0)) - m(VectorValue(0.0)),cmaps)),Γ)
# time: 2023-04-04 02:32:43 CEST
# mode: julia
	tangent_cf = CellField(lazy_map(t->t/norm(t),lazy_map(m->m(VectorValue(1.0)) - m(VectorValue(0.0)),cmaps)),Γ)
# time: 2023-04-04 02:33:50 CEST
# mode: julia
	x = [1,2,3]
# time: 2023-04-04 02:34:02 CEST
# mode: julia
	x .+ 1
# time: 2023-04-04 02:34:39 CEST
# mode: julia
	map(xi->xi+1,x)
# time: 2023-04-04 02:35:43 CEST
# mode: julia
	my_tangent(m)
# time: 2023-04-04 02:44:48 CEST
# mode: julia
	m = cmaps[1]
# time: 2023-04-04 03:10:20 CEST
# mode: julia
	ones(2)
# time: 2023-04-04 03:10:38 CEST
# mode: julia
	ones(ones(2)
	)
# time: 2023-04-04 03:10:51 CEST
# mode: julia
	tr(ones(2))
# time: 2023-04-04 03:10:54 CEST
# mode: julia
	trc(ones(2))
# time: 2023-04-04 03:11:22 CEST
# mode: julia
	tr(ones(2))
# time: 2023-04-04 03:11:36 CEST
# mode: julia
	using Gridap.\Field
# time: 2023-04-04 03:11:39 CEST
# mode: julia
	using Gridap.Field
# time: 2023-04-04 03:11:42 CEST
# mode: julia
	using Gridap.Fieldd
# time: 2023-04-04 03:11:45 CEST
# mode: julia
	using Gridap.Fields
# time: 2023-04-04 03:11:48 CEST
# mode: julia
	tr(ones(2))
# time: 2023-04-04 03:32:20 CEST
# mode: julia
	num_components(SymFourthOrderTensorValue(1111,1121,1122, 2111,2121,2122, 2211,2221,2222)) == 16
# time: 2023-04-04 03:32:24 CEST
# mode: julia
	num_components(SymFourthOrderTensorValue(1111,1121,1122, 2111,2121,2122, 2211,2221,2222))
# time: 2023-04-04 03:32:30 CEST
# mode: julia
	num_components(SymFourthOrderTensorValue(1111))
# time: 2023-04-04 03:33:28 CEST
# mode: julia
	rand(SymFourthOrderTensorValue{3,Int})
# time: 2023-04-04 03:33:40 CEST
# mode: julia
	rand(SymFourthOrderTensorValue{3,Int})
# time: 2023-04-04 03:33:58 CEST
# mode: julia
	num_components( rand(SymFourthOrderTensorValue{3,Int}) )
# time: 2023-04-04 03:34:03 CEST
# mode: julia
	num_components( rand(SymFourthOrderTensorValue{1,Int}) )
# time: 2023-04-04 03:34:06 CEST
# mode: julia
	num_components( rand(SymFourthOrderTensorValue{2,Int}) )
# time: 2023-04-04 03:34:45 CEST
# mode: julia
	num_components( rand(SymFourthOrderTensorValue{3,Int}) )
# time: 2023-04-04 04:01:06 CEST
# mode: julia
	CTf
# time: 2023-04-04 04:01:12 CEST
# mode: julia
	CTs
# time: 2023-04-04 04:50:39 CEST
# mode: julia
	SymFourthOrderTensorValue( (3.0) )
# time: 2023-04-04 04:50:45 CEST
# mode: julia
	A = SymFourthOrderTensorValue( (3.0) )
# time: 2023-04-04 04:50:53 CEST
# mode: julia
	toarray(A)
# time: 2023-04-04 04:51:19 CEST
# mode: julia
	get_array(A)
# time: 2023-04-04 05:10:24 CEST
# mode: julia
	matmat = [6.41555556e6 0.0; 0.0 2.24822222e6]
# time: 2023-04-04 05:10:37 CEST
# mode: julia
	tuptup = (0.0, 0.0)
# time: 2023-04-04 05:10:43 CEST
# mode: julia
	matmat\cdot
# time: 2023-04-04 05:11:05 CEST
# mode: julia
	arrDd⋅ε
# time: 2023-04-04 05:11:15 CEST
# mode: julia
	matmat⋅tuptup
# time: 2023-04-04 05:11:30 CEST
# mode: julia
	lennght(tuptup)
# time: 2023-04-04 05:11:34 CEST
# mode: julia
	lenght(tuptup)
# time: 2023-04-04 05:11:50 CEST
# mode: julia
	shape(tuptup)
# time: 2023-04-04 05:11:53 CEST
# mode: julia
	shapes(tuptup)
# time: 2023-04-04 05:11:59 CEST
# mode: julia
	sizes(tuptup)
# time: 2023-04-04 05:12:03 CEST
# mode: julia
	tuptup
# time: 2023-04-04 07:18:22 CEST
# mode: julia
	exit()
# time: 2023-04-04 07:50:29 CEST
# mode: julia
	CTf
# time: 2023-04-04 09:37:58 CEST
# mode: julia
	CTf[4]
# time: 2023-04-04 09:38:04 CEST
# mode: julia
	CTf[4][1]
# time: 2023-04-04 09:38:16 CEST
# mode: julia
	ct1
# time: 2023-04-04 09:38:42 CEST
# mode: julia
	ct2
# time: 2023-04-05 04:37:23 CEST
# mode: julia
	tangent_cf.cell_field.maps
# time: 2023-04-05 04:37:30 CEST
# mode: julia
	tangent_cf.cell_field.maps[1]
# time: 2023-04-05 04:37:36 CEST
# mode: julia
	tangent_cf[1]
# time: 2023-04-05 04:42:40 CEST
# mode: julia
	[1.0 2.0; 3.0 4.0]
# time: 2023-04-05 04:43:02 CEST
# mode: julia
	[[1.0, 2.0],[3.0, 4.0]]
# time: 2023-04-05 08:44:36 CEST
# mode: julia
	A = first(contr.dict)
# time: 2023-04-05 08:44:42 CEST
# mode: julia
	A[1]
# time: 2023-04-05 08:44:45 CEST
# mode: julia
	A[1,1]
# time: 2023-04-05 08:44:50 CEST
# mode: julia
	A
# time: 2023-04-05 08:45:39 CEST
# mode: julia
	contr
# time: 2023-04-05 08:45:43 CEST
# mode: julia
	contr.dict
# time: 2023-04-05 08:46:01 CEST
# mode: julia
	first(contr.dict)
# time: 2023-04-05 08:47:15 CEST
# mode: julia
	second(contr.dict)
# time: 2023-04-05 08:47:28 CEST
# mode: julia
	A.second
# time: 2023-04-05 08:47:32 CEST
# mode: julia
	A.second[1]
# time: 2023-04-05 08:47:44 CEST
# mode: julia
	A.second[1][1,1]
# time: 2023-04-05 08:47:48 CEST
# mode: julia
	A.second[1][3,1]
# time: 2023-04-05 08:47:58 CEST
# mode: julia
	A.second[1][2,2]
# time: 2023-04-05 08:48:11 CEST
# mode: julia
	A.second[1][:,:]
# time: 2023-04-05 08:48:16 CEST
# mode: julia
	A.second[1][:]
# time: 2023-04-05 08:48:20 CEST
# mode: julia
	A.second[1][2,2]
# time: 2023-04-05 08:48:24 CEST
# mode: julia
	A.second[1]
# time: 2023-04-05 08:48:40 CEST
# mode: julia
	A.second[1][3,1]
# time: 2023-04-05 08:48:45 CEST
# mode: julia
	A.second[1][2,2]
# time: 2023-04-05 08:48:49 CEST
# mode: julia
	A.second[1][3,2]
# time: 2023-04-05 08:49:16 CEST
# mode: julia
	A.second[1][1,3]
# time: 2023-04-05 08:49:22 CEST
# mode: julia
	A.second[1][2,3]
# time: 2023-04-05 09:04:20 CEST
# mode: julia
	A.second[1][3,1]
# time: 2023-04-05 09:33:24 CEST
# mode: julia
	elements[1]
# time: 2023-04-05 09:33:37 CEST
# mode: julia
	elements.second
# time: 2023-04-05 09:33:40 CEST
# mode: julia
	elements.second[1]
# time: 2023-04-05 09:33:52 CEST
# mode: julia
	elements.second[1][1,1]
# time: 2023-04-05 09:33:56 CEST
# mode: julia
	elements.second[1][2,2]
# time: 2023-04-05 09:34:07 CEST
# mode: julia
	elements.second[1][2,3]
# time: 2023-04-05 09:34:11 CEST
# mode: julia
	elements.second[1][3,2]
# time: 2023-04-05 09:34:14 CEST
# mode: julia
	elements.second[1][3,3]
# time: 2023-04-05 09:35:54 CEST
# mode: julia
	ct1
# time: 2023-04-05 09:37:44 CEST
# mode: julia
	elements.second[1][1,1]
# time: 2023-04-05 09:37:47 CEST
# mode: julia
	elements.second[1][2,2]
# time: 2023-04-05 09:38:27 CEST
# mode: julia
	elements.second[1][2,3]
# time: 2023-04-05 09:38:32 CEST
# mode: julia
	elements.second[1][3,2]
# time: 2023-04-05 09:38:34 CEST
# mode: julia
	elements.second[1][3,3]
# time: 2023-04-05 09:47:57 CEST
# mode: julia
	elements.second[1][1,1]
# time: 2023-04-05 09:48:08 CEST
# mode: julia
	elements.second[1][2,2]
# time: 2023-04-05 09:48:19 CEST
# mode: julia
	elements.second[1][2,3]
# time: 2023-04-05 09:48:35 CEST
# mode: julia
	elements.second[1][3,2]
# time: 2023-04-05 09:48:44 CEST
# mode: julia
	elements.second[1][3,3]
# time: 2023-04-05 09:49:40 CEST
# mode: julia
	ct1
# time: 2023-04-05 09:56:55 CEST
# mode: julia
	elements.second[1][3,3]
# time: 2023-04-05 09:57:09 CEST
# mode: julia
	elements.second[1][3,2]
# time: 2023-04-05 09:58:07 CEST
# mode: julia
	elements.second[1][1,1]
# time: 2023-04-05 09:58:14 CEST
# mode: julia
	elements.second[1][3,3]
# time: 2023-04-05 09:59:20 CEST
# mode: julia
	elements.second[1][2,2]
# time: 2023-04-05 09:59:23 CEST
# mode: julia
	elements.second[1][2,3]
# time: 2023-04-05 09:59:27 CEST
# mode: julia
	elements.second[1][3,2]
# time: 2023-04-05 09:59:29 CEST
# mode: julia
	elements.second[1][3,3]
# time: 2023-04-05 10:19:17 CEST
# mode: julia
	element1
# time: 2023-04-05 10:19:20 CEST
# mode: julia
	element1[3,3]
# time: 2023-04-05 10:29:19 CEST
# mode: julia
	elementa[3,3]
# time: 2023-04-05 10:29:28 CEST
# mode: julia
	elementD[3,3]
# time: 2023-04-05 10:29:44 CEST
# mode: julia
	elementD
# time: 2023-04-05 10:30:04 CEST
# mode: julia
	elementS
# time: 2023-04-05 10:31:03 CEST
# mode: julia
	elementD
# time: 2023-04-05 10:31:09 CEST
# mode: julia
	element
# time: 2023-04-05 10:31:12 CEST
# mode: julia
	element[3,3]
# time: 2023-04-05 10:31:47 CEST
# mode: julia
	elementS
# time: 2023-04-05 10:31:53 CEST
# mode: julia
	elementD
# time: 2023-04-05 10:33:16 CEST
# mode: julia
	elementS[2,2]
# time: 2023-04-05 10:33:26 CEST
# mode: julia
	elementS[2,3]
# time: 2023-04-05 10:33:31 CEST
# mode: julia
	elementS[3,3]
# time: 2023-04-05 10:34:31 CEST
# mode: julia
	elementS[3,3]-element[3,3]
# time: 2023-04-05 10:34:39 CEST
# mode: julia
	elementS[3,3]-element
# time: 2023-04-05 10:34:47 CEST
# mode: julia
	elementS[3,3]+element
# time: 2023-04-05 10:35:03 CEST
# mode: julia
	elementD
# time: 2023-04-05 10:35:09 CEST
# mode: julia
	elementS[3,3]+elementD
# time: 2023-04-05 10:35:38 CEST
# mode: julia
	element
# time: 2023-04-05 13:25:02 CEST
# mode: julia
	ct1
# time: 2023-04-05 14:08:10 CEST
# mode: julia
	tenA = SymFourthOrderTensorValue{2, Float64, 9}(791208.7912087912, 0.0, 237362.63736263735, 0.0, 276923.07692307694, 0.0, 237362.63736263735, 0.0, 791208.7912087912)
# time: 2023-04-05 14:09:01 CEST
# mode: julia
	tenA = SymFourthOrderTensorValue((791208.7912087912, 0.0, 237362.63736263735, 0.0, 276923.07692307694, 0.0, 237362.63736263735, 0.0, 791208.7912087912))
# time: 2023-04-05 14:09:35 CEST
# mode: julia
	tenB = SymFourthOrderTensorValue((0.01, 0.0, 0.0))
# time: 2023-04-05 14:10:08 CEST
# mode: julia
	tenB = SymFourthOrderTensorValue(0.01, 0.0, 0.0)
# time: 2023-04-05 14:10:17 CEST
# mode: julia
	tenB = SymFourthOrderTensorValue((0.01, 0.0, 0.0))
# time: 2023-04-05 14:10:30 CEST
# mode: julia
	tenB = SymTensorValue((0.01, 0.0, 0.0))
# time: 2023-04-05 14:11:04 CEST
# mode: julia
	tenA ⊙ tenB
# time: 2023-04-05 14:11:31 CEST
# mode: julia
	tenB = SymTensorValue((0.01, 0.0, -0.003))
# time: 2023-04-05 14:11:33 CEST
# mode: julia
	tenA ⊙ tenB
# time: 2023-04-05 14:11:45 CEST
# mode: julia
	tenB = SymTensorValue((0.0, 0.01, 0.0))
# time: 2023-04-05 14:11:48 CEST
# mode: julia
	tenA ⊙ tenB
# time: 2023-04-06 01:15:49 CEST
# mode: julia
	tf
# time: 2023-04-06 01:15:58 CEST
# mode: julia
	tf.cell_field
# time: 2023-04-06 01:16:04 CEST
# mode: julia
	tf.cell_field.args
# time: 2023-04-06 01:16:13 CEST
# mode: julia
	tf.cell_field.args[1]
# time: 2023-04-06 01:19:08 CEST
# mode: julia
	tf.cell_field[1]
# time: 2023-04-06 01:19:13 CEST
# mode: julia
	tf.cell_field[1].value
# time: 2023-04-06 01:19:26 CEST
# mode: julia
	getL(tf.cell_field[1])
# time: 2023-04-06 01:19:35 CEST
# mode: julia
	getL(tf)
# time: 2023-04-06 01:19:39 CEST
# mode: julia
	B = getL(tf)
# time: 2023-04-06 01:19:45 CEST
# mode: julia
	B.args
# time: 2023-04-06 01:19:56 CEST
# mode: julia
	B.args[1]
# time: 2023-04-06 01:36:58 CEST
# mode: julia
	B = getL(tf)
# time: 2023-04-06 01:37:59 CEST
# mode: julia
	tf
# time: 2023-04-06 01:38:07 CEST
# mode: julia
	tf.trian
# time: 2023-04-06 01:38:18 CEST
# mode: julia
	tf.cell_field.args
# time: 2023-04-06 01:38:30 CEST
# mode: julia
	tf.cell_field.maps
# time: 2023-04-06 01:38:43 CEST
# mode: julia
	tf.cell_field.args
# time: 2023-04-06 01:39:36 CEST
# mode: julia
	B = getL(tf)
# time: 2023-04-06 01:40:45 CEST
# mode: julia
	B = ∘(getL(tf))
# time: 2023-04-06 01:41:01 CEST
# mode: julia
	B[1]
# time: 2023-04-06 01:41:34 CEST
# mode: julia
	B = ∘( getL(tf) )
# time: 2023-04-06 01:41:39 CEST
# mode: julia
	BB = ∘( getL(tf) )
# time: 2023-04-06 02:04:23 CEST
# mode: julia
	SymFourthOrderTensorValue( (1,2,3,4,5,6,7,8,9) )
# time: 2023-04-06 02:04:50 CEST
# mode: julia
	A = SymFourthOrderTensorValue( (1,2,3,4,5,6,7,8,9) )
# time: 2023-04-06 02:04:57 CEST
# mode: julia
	ten33 = SymFourthOrderTensorValue( (1,2,3,4,5,6,7,8,9) )
# time: 2023-04-06 02:05:11 CEST
# mode: julia
	toarray(ten33)
# time: 2023-04-06 02:05:18 CEST
# mode: julia
	toArray(ten33)
# time: 2023-04-06 02:05:43 CEST
# mode: julia
	get_array(ten33)
# time: 2023-04-06 02:05:53 CEST
# mode: julia
	ten33 = SymFourthOrderTensorValue( (1,2,3,4,5,6,7,8,9) )
# time: 2023-04-06 02:06:17 CEST
# mode: julia
	ten33.data
# time: 2023-04-06 02:06:33 CEST
# mode: julia
	get_array(ten33)
# time: 2023-04-06 02:06:47 CEST
# mode: julia
	TensorValue(ten33)
# time: 2023-04-06 02:06:55 CEST
# mode: julia
	ten22 = TensorValue(ten33)
# time: 2023-04-06 02:06:58 CEST
# mode: julia
	get_array(ten22)
# time: 2023-04-06 02:07:22 CEST
# mode: julia
	Matrix(ten33)
# time: 2023-04-06 02:10:34 CEST
# mode: julia
	A = SymFourthOrderTensorValue( (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) )
# time: 2023-04-06 02:11:33 CEST
# mode: julia
	A = SymFourthOrderTensorValue( (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21) )
# time: 2023-04-06 02:16:20 CEST
# mode: julia
	A = SymFourthOrderTensorValue( (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36) )
# time: 2023-04-06 02:16:29 CEST
# mode: julia
	Aaaaa = SymFourthOrderTensorValue( (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36) )
# time: 2023-04-06 02:17:03 CEST
# mode: julia
	Matrix(Aaaaa)
# time: 2023-04-06 02:17:10 CEST
# mode: julia
	get_array(Aaaaa)
# time: 2023-04-06 02:17:54 CEST
# mode: julia
	SymTensorValue(1,2,3)
# time: 2023-04-06 02:18:06 CEST
# mode: julia
	Matrix(SymTensorValue(1,2,3))
# time: 2023-04-06 02:18:21 CEST
# mode: julia
	get_array(SymTensorValue(1,2,3))
# time: 2023-04-06 02:18:38 CEST
# mode: julia
	get_array(SymTensorValue(1,2,3,4,5,6,7.8.9))
# time: 2023-04-06 02:18:42 CEST
# mode: julia
	get_array(SymTensorValue(1,2,3,4,5,6,7,8,9))
# time: 2023-04-06 02:18:53 CEST
# mode: julia
	get_array(SymTensorValue(1,2,3,4,5,6))
# time: 2023-04-06 02:48:18 CEST
# mode: julia
	uf
# time: 2023-04-06 02:48:32 CEST
# mode: julia
	ε(uf)
# time: 2023-04-06 02:48:43 CEST
# mode: julia
	epsf = ε(uf)
# time: 2023-04-06 02:48:52 CEST
# mode: julia
	epsf.args
# time: 2023-04-06 02:50:29 CEST
# mode: julia
	u_data = Gridap.CellData.get_data(u)
# time: 2023-04-06 02:50:40 CEST
# mode: julia
	u_data = Gridap.CellData.get_data(uf)
# time: 2023-04-06 02:50:48 CEST
# mode: julia
	u_data.args
# time: 2023-04-06 02:50:54 CEST
# mode: julia
	u_data.args[1]
# time: 2023-04-06 02:51:01 CEST
# mode: julia
	u_data.args[1].values
# time: 2023-04-06 02:51:58 CEST
# mode: julia
	tf_data = Gridap.CellData.get_data(tf)
# time: 2023-04-06 02:52:55 CEST
# mode: julia
	tf.cell_field.args
# time: 2023-04-06 02:53:22 CEST
# mode: julia
	g2l(tf,uf)
# time: 2023-04-06 02:53:33 CEST
# mode: julia
	GGG = g2l(tf,uf)
# time: 2023-04-06 04:17:51 CEST
# mode: julia
	u = FEFunction(0.0,U)
# time: 2023-04-06 04:18:02 CEST
# mode: julia
	u = FEFunction(x->0.0,U)
# time: 2023-04-06 04:18:13 CEST
# mode: julia
	u = zero(U)
# time: 2023-04-06 04:18:21 CEST
# mode: julia
	u = zero(U)[1]
# time: 2023-04-06 04:18:31 CEST
# mode: julia
	gl2(tf,u)
# time: 2023-04-06 04:18:41 CEST
# mode: julia
	g2l(tf,u)
# time: 2023-04-06 04:18:51 CEST
# mode: julia
	tfu = g2l(tf,u)
# time: 2023-04-06 04:19:02 CEST
# mode: julia
	∂(tfu)
# time: 2023-04-06 04:19:14 CEST
# mode: julia
	dtfu = ∂(tfu)
# time: 2023-04-06 04:20:47 CEST
# mode: julia
	CTf
# time: 2023-04-06 04:21:19 CEST
# mode: julia
	dtfu
# time: 2023-04-06 04:21:38 CEST
# mode: julia
	σ(CTF[1],dtfu)
# time: 2023-04-06 04:21:43 CEST
# mode: julia
	σ(CTf[1],dtfu)
# time: 2023-04-06 04:22:07 CEST
# mode: julia
	τ(CTf[1],dtfu)
# time: 2023-04-06 04:22:09 CEST
# mode: julia
	σ(CTf[1],dtfu)
# time: 2023-04-06 04:22:54 CEST
# mode: julia
	typeof(dtfu)
# time: 2023-04-06 04:23:08 CEST
# mode: julia
	typeof(CTf[1])
# time: 2023-04-06 04:23:19 CEST
# mode: julia
	σ(CTf[1],dtfu)
# time: 2023-04-06 04:25:36 CEST
# mode: julia
	dtfu
# time: 2023-04-06 04:25:44 CEST
# mode: julia
	tfu
# time: 2023-04-06 04:26:01 CEST
# mode: julia
	pts = get_cell_points(dΩ)
# time: 2023-04-06 04:26:05 CEST
# mode: julia
	tfu(pts)
# time: 2023-04-06 04:26:20 CEST
# mode: julia
	aux = tfu(pts)
# time: 2023-04-06 04:26:22 CEST
# mode: julia
	aux[1]
# time: 2023-04-06 04:26:44 CEST
# mode: julia
	aux[2]
# time: 2023-04-06 04:26:47 CEST
# mode: julia
	aux[3]
# time: 2023-04-06 04:27:00 CEST
# mode: julia
	aux = tfu(pts)
# time: 2023-04-06 04:27:16 CEST
# mode: julia
	tfu
# time: 2023-04-06 04:27:24 CEST
# mode: julia
	aux = tfu(pts)
# time: 2023-04-06 04:27:29 CEST
# mode: julia
	aux[1]
# time: 2023-04-06 04:28:44 CEST
# mode: julia
	tfu = g2l(tf,u)
# time: 2023-04-06 04:28:52 CEST
# mode: julia
	tfu(aux)
# time: 2023-04-06 04:28:57 CEST
# mode: julia
	tfu(pts)
# time: 2023-04-06 04:30:55 CEST
# mode: julia
	tfu = g2l(tf,u)
# time: 2023-04-06 04:30:58 CEST
# mode: julia
	tfu(pts)
# time: 2023-04-06 04:31:56 CEST
# mode: julia
	σ(CTf[1],dtfu)
# time: 2023-04-06 04:32:03 CEST
# mode: julia
	dtfu = ∂(tfu)
# time: 2023-04-06 04:32:25 CEST
# mode: julia
	tfu = g2l(tf,u)
# time: 2023-04-06 04:32:26 CEST
# mode: julia
	dtfu = ∂(tfu)
# time: 2023-04-06 04:33:41 CEST
# mode: julia
	tfu = g2l(tf,u)
# time: 2023-04-06 04:33:43 CEST
# mode: julia
	dtfu = ∂(tfu)
# time: 2023-04-06 04:39:55 CEST
# mode: julia
	du_p = tangential_derivative(u,tf)
# time: 2023-04-06 04:40:02 CEST
# mode: julia
	du_p(pts)
# time: 2023-04-06 04:40:15 CEST
# mode: julia
	aux = du_p(pts)
# time: 2023-04-06 04:40:18 CEST
# mode: julia
	aux[1]
# time: 2023-04-06 04:46:19 CEST
# mode: julia
	du_p = g2l(u,tf)
# time: 2023-04-06 04:46:25 CEST
# mode: julia
	aux = du_p(pts)
# time: 2023-04-06 04:50:42 CEST
# mode: julia
	du_p = g2l(u,tf)
# time: 2023-04-06 04:50:52 CEST
# mode: julia
	CTF[1]
# time: 2023-04-06 04:50:55 CEST
# mode: julia
	CTf[1]
# time: 2023-04-06 04:51:12 CEST
# mode: julia
	τ(CTf[1],dtu)
# time: 2023-04-06 04:51:21 CEST
# mode: julia
	τ(CTf[1],du_p)
# time: 2023-04-06 04:51:35 CEST
# mode: julia
	σ(CTf[1],du_p)
# time: 2023-04-06 04:51:51 CEST
# mode: julia
	a = CTf[1](pts)[1]
# time: 2023-04-06 04:51:56 CEST
# mode: julia
	aux1 = CTf[1](pts)[1]
# time: 2023-04-06 04:52:20 CEST
# mode: julia
	aux2 = du_p(pts)[1]
# time: 2023-04-06 04:52:24 CEST
# mode: julia
	aux2 = du_p(pts)[1][1]
# time: 2023-04-06 04:52:27 CEST
# mode: julia
	aux1 = CTf[1](pts)[1][1]
# time: 2023-04-06 04:52:57 CEST
# mode: julia
	aux1⊙aux2
# time: 2023-04-06 04:53:04 CEST
# mode: julia
	aux
# time: 2023-04-06 04:53:08 CEST
# mode: julia
	aux1
# time: 2023-04-06 04:53:14 CEST
# mode: julia
	aux2
# time: 2023-04-06 04:53:41 CEST
# mode: julia
	SymTensorValue(aux1)⊙aux2
# time: 2023-04-06 04:53:53 CEST
# mode: julia
	SymTensorValue(aux1.data)⊙aux2
# time: 2023-04-06 04:54:04 CEST
# mode: julia
	SymTensorValue(aux1.data)⋅aux2
# time: 2023-04-06 04:54:15 CEST
# mode: julia
	SymTensorValue(aux1)⋅aux2
# time: 2023-04-06 04:56:37 CEST
# mode: julia
	aux1⊙aux2
# time: 2023-04-06 04:56:41 CEST
# mode: julia
	aux1⋅aux2
# time: 2023-04-06 04:56:56 CEST
# mode: julia
	SymTensorValue(aux1.data)⋅aux2
# time: 2023-04-06 04:59:57 CEST
# mode: julia
	TensorValue(aux1.data)⊙aux2
# time: 2023-04-06 05:00:09 CEST
# mode: julia
	VectorValue(aux1.data)⊙aux2
# time: 2023-04-06 05:13:51 CEST
# mode: julia
	6.5⋅SymTensorValue((2.0))
# time: 2023-04-06 05:15:52 CEST
# mode: julia
	TensorValue((2.0))⊙5.0
# time: 2023-04-06 05:16:58 CEST
# mode: julia
	TensorValue((2.0))⊙SymTensorValue((5.0))
# time: 2023-04-06 05:17:17 CEST
# mode: julia
	4.0⊙SymTensorValue((5.0))
# time: 2023-04-06 05:17:31 CEST
# mode: julia
	SymTensorValue((5.0))⊙5.0
# time: 2023-04-06 05:17:56 CEST
# mode: julia
	4.0⊙SymTensorValue((5.0))
# time: 2023-04-06 05:18:09 CEST
# mode: julia
	4.0⊙TensorValue((5.0))
# time: 2023-04-06 05:18:18 CEST
# mode: julia
	4.0⊙TensorValue((5.0)).data
# time: 2023-04-06 06:53:31 CEST
# mode: julia
	u = zero(U)
# time: 2023-04-06 06:53:35 CEST
# mode: julia
	u = zero(U)[1]
# time: 2023-04-06 06:56:50 CEST
# mode: julia
	pts = get_cell_point(dΩ)
# time: 2023-04-06 06:56:55 CEST
# mode: julia
	pts = get_cell_point(dΩ.quad)
# time: 2023-04-06 06:56:58 CEST
# mode: julia
	pts = get_cell_points(dΩ.quad)
# time: 2023-04-06 06:57:05 CEST
# mode: julia
	u(pts)
# time: 2023-04-06 06:57:23 CEST
# mode: julia
	∂(u)(pts)
# time: 2023-04-06 06:57:32 CEST
# mode: julia
	∂(u)(pts)[1]
# time: 2023-04-06 06:58:10 CEST
# mode: julia
	CTf
# time: 2023-04-06 06:58:18 CEST
# mode: julia
	CTf[1]
# time: 2023-04-06 06:58:39 CEST
# mode: julia
	σ(CTf[1],∂(u))
# time: 2023-04-06 06:58:51 CEST
# mode: julia
	u = zero(U)[1]
# time: 2023-04-06 06:58:55 CEST
# mode: julia
	pts = get_cell_points(dΩ.quad)
# time: 2023-04-06 06:58:57 CEST
# mode: julia
	σ(CTf[1],∂(u))
# time: 2023-04-06 06:59:05 CEST
# mode: julia
	σ(CTf[1],∂(u))(pts)
# time: 2023-04-06 06:59:13 CEST
# mode: julia
	σ(CTf[1],∂(u))(pts)[1]
# time: 2023-04-06 07:03:06 CEST
# mode: julia
	inner(CTf[1],∂(u))
# time: 2023-04-06 07:03:21 CEST
# mode: julia
	inner(∂(u),inner(CTf[1],∂(u)))
# time: 2023-04-06 07:03:27 CEST
# mode: julia
	inner(∂(u),inner(CTf[1],∂(u)))(pts)[1]
# time: 2023-04-06 07:11:50 CEST
# mode: julia
	u = zero(U)
# time: 2023-04-06 07:11:57 CEST
# mode: julia
	pts = get_cell_points(dΩ.quad)
# time: 2023-04-06 07:12:05 CEST
# mode: julia
	σ(CTf[1],∂(u))
# time: 2023-04-06 07:12:11 CEST
# mode: julia
	σ(CTf[1],∂(u))(pts)
# time: 2023-04-06 07:12:26 CEST
# mode: julia
	σ(CTf[1],∂(u))(pts)[1]
# time: 2023-04-06 07:12:49 CEST
# mode: julia
	pts = get_cell_points(dΩ.quad)
# time: 2023-04-06 07:12:53 CEST
# mode: julia
	σ(CTf[1],∂(u))(pts)
# time: 2023-04-06 07:13:00 CEST
# mode: julia
	u = zero(U)[1]
# time: 2023-04-06 07:13:04 CEST
# mode: julia
	σ(CTf[1],∂(u))(pts)
# time: 2023-04-06 07:13:19 CEST
# mode: julia
	inner(∂(u),inner(CTf[1],∂(u)))(pts)[1]
# time: 2023-04-06 08:14:31 CEST
# mode: julia
	ssss = tangential_derivative_1(u,tf)
# time: 2023-04-06 08:14:38 CEST
# mode: julia
	u = zero(U)[1]
# time: 2023-04-06 08:14:41 CEST
# mode: julia
	ssss = tangential_derivative_1(u,tf)
# time: 2023-04-06 08:15:03 CEST
# mode: julia
	ssss[1]
# time: 2023-04-06 08:26:01 CEST
# mode: julia
	VectorValue(1.0,2.0)
# time: 2023-04-06 08:26:11 CEST
# mode: julia
	tt = VectorValue(1.0,2.0)
# time: 2023-04-06 08:26:22 CEST
# mode: julia
	gradient(tt)
# time: 2023-04-06 08:26:33 CEST
# mode: julia
	yyy = gradient(tt)
# time: 2023-04-06 08:26:35 CEST
# mode: julia
	yyy
# time: 2023-04-06 08:26:41 CEST
# mode: julia
	yyy.f
# time: 2023-04-06 08:26:59 CEST
# mode: julia
	gradient(tt)⋅gradient(tt)
# time: 2023-04-06 08:27:21 CEST
# mode: julia
	gradient(tt)⊙gradient(tt)
# time: 2023-04-06 08:27:50 CEST
# mode: julia
	ε(tt)⊙ε(tt)
# time: 2023-04-06 08:27:55 CEST
# mode: julia
	ε(tt)\cdotε(tt)
# time: 2023-04-06 08:27:58 CEST
# mode: julia
	ε(tt)⋅ε(tt)
# time: 2023-04-06 08:39:55 CEST
# mode: julia
	ε(tt)∈ε(tt)
# time: 2023-04-06 08:40:02 CEST
# mode: julia
	gradient(tt)⊙gradient(tt)
# time: 2023-04-06 08:40:22 CEST
# mode: julia
	gradient(tt)⋅VectorValue(1.0)
# time: 2023-04-06 08:40:35 CEST
# mode: julia
	tt = VectorValue(1.0,2.0)
# time: 2023-04-06 08:40:47 CEST
# mode: julia
	tt = VectorValue(1.0)
# time: 2023-04-06 08:41:02 CEST
# mode: julia
	gradient(tt)
# time: 2023-04-06 08:41:36 CEST
# mode: julia
	gradient(tt)∘gradient(tt)
# time: 2023-04-06 08:41:47 CEST
# mode: julia
	ii = gradient(tt)∘gradient(tt)
# time: 2023-04-06 08:43:15 CEST
# mode: julia
	VV = get_fe_basis(Vt)
# time: 2023-04-06 08:43:27 CEST
# mode: julia
	gradient(VV)
# time: 2023-04-06 08:43:38 CEST
# mode: julia
	ww = gradient(VV)
# time: 2023-04-06 08:43:49 CEST
# mode: julia
	ww.cell_basis.args
# time: 2023-04-06 08:43:55 CEST
# mode: julia
	ww.cell_basis.args[1]
# time: 2023-04-06 08:43:59 CEST
# mode: julia
	ww.cell_basis.args[1][1]
# time: 2023-04-12 03:39:31 CEST
# mode: julia
	VectorValue{Int8, 4}()
# time: 2023-04-12 03:39:39 CEST
# mode: julia
	AAA = VectorValue{Int8, 4}()
# time: 2023-04-12 03:39:46 CEST
# mode: julia
	AAA = VectorValue{Int8, 4}(3.0, 2.0)
# time: 2023-04-12 03:39:52 CEST
# mode: julia
	AAA = VectorValue{Int8, 2}(3.0, 2.0)
# time: 2023-04-12 03:40:01 CEST
# mode: julia
	AAA = VectorValue{Int8, 2}(3, 2)
# time: 2023-04-12 03:40:41 CEST
# mode: julia
	AAA = VectorValue{2, Int8}(3, 2)
# time: 2023-04-12 03:40:46 CEST
# mode: julia
	AAA
# time: 2023-04-12 09:01:35 CEST
# mode: julia
	import Pkg; Pkg.add("DataStructures")
# time: 2023-04-17 03:59:33 CEST
# mode: julia
	a = MMatrix{2,2}(1,2,3,4)
# time: 2023-04-17 03:59:37 CEST
# mode: julia
	t = TensorValue(a)
# time: 2023-04-17 03:59:42 CEST
# mode: julia
	using Gridap
# time: 2023-04-17 03:59:52 CEST
# mode: julia
	a = MMatrix{2,2}(1,2,3,4)
# time: 2023-04-17 04:09:35 CEST
# mode: julia
	t = TensorValue(1.0,2.0,3.0,4.0)
# time: 2023-04-17 04:09:46 CEST
# mode: julia
	t't
# time: 2023-04-17 04:09:53 CEST
# mode: julia
	t'⋅t
# time: 2023-04-17 04:10:05 CEST
# mode: julia
	transpose(t)⋅t
# time: 2023-04-17 04:10:13 CEST
# mode: julia
	tt = transpose(t)⋅t
# time: 2023-04-17 04:10:27 CEST
# mode: julia
	SymTensorValue(tt)
# time: 2023-04-17 04:11:23 CEST
# mode: julia
	SymTensorValue(get_array(tt))
# time: 2023-04-17 04:11:34 CEST
# mode: julia
	tt = SymTensorValue(get_array(tt))
# time: 2023-04-17 04:11:40 CEST
# mode: julia
	tt = transpose(t)⋅t
# time: 2023-04-17 04:11:49 CEST
# mode: julia
	ttt = SymTensorValue(get_array(tt))
# time: 2023-04-17 04:11:52 CEST
# mode: julia
	display(ttt)
# time: 2023-04-17 04:11:59 CEST
# mode: julia
	Matrix(ttt)
# time: 2023-04-17 04:12:09 CEST
# mode: julia
	get_matrix(ttt)
# time: 2023-04-17 04:12:12 CEST
# mode: julia
	ttt.data
# time: 2023-04-17 04:14:01 CEST
# mode: julia
	t
# time: 2023-04-17 04:14:06 CEST
# mode: julia
	symmetrize(t)
# time: 2023-04-17 04:15:24 CEST
# mode: julia
	symmetrize(TensorValue(0.3,0.0,0.5,0.0))
# time: 2023-04-17 04:57:20 CEST
# mode: julia
	TensorValue([ 1.0 2.0; 3.0 4.0; 5.0 6.0 ])
# time: 2023-04-17 04:57:31 CEST
# mode: julia
	myt = TensorValue([ 1.0 2.0; 3.0 4.0; 5.0 6.0 ])
# time: 2023-04-17 04:57:39 CEST
# mode: julia
	myt[1,1]
# time: 2023-04-17 04:57:42 CEST
# mode: julia
	myt[1,2]
# time: 2023-04-17 04:57:44 CEST
# mode: julia
	myt[1,3]
# time: 2023-04-17 04:57:48 CEST
# mode: julia
	myt[3,1]
# time: 2023-04-17 04:57:50 CEST
# mode: julia
	myt[3,2]
# time: 2023-04-17 04:58:23 CEST
# mode: julia
	myt2 = TensorValue([ 1.0 2.0 9.0; 3.0 4.0 7.0; 5.0 6.0 8.0 ])
# time: 2023-04-17 04:58:34 CEST
# mode: julia
	myt2 ⋅ myt
# time: 2023-04-17 06:20:04 CEST
# mode: julia
	zero(SymTensorValue{3,Int})
# time: 2023-04-17 06:20:54 CEST
# mode: julia
	rand(SymTensorValue{3,Int})
# time: 2023-04-17 06:21:04 CEST
# mode: julia
	AAA = rand(SymTensorValue{2,Int})
# time: 2023-04-17 06:21:07 CEST
# mode: julia
	AAA
# time: 2023-04-17 06:22:53 CEST
# mode: julia
	AAA = rand(SymTensorValue{3,Int})
# time: 2023-04-17 06:23:14 CEST
# mode: julia
	AAA = rand(SymFourOrderTensorValue{3,Int})
# time: 2023-04-17 06:23:26 CEST
# mode: julia
	AAA = rand(SymFourthOrderTensorValue{3,Int})
# time: 2023-04-17 06:23:35 CEST
# mode: julia
	AAA = rand(SymFourthOrderTensorValue{2,Int})
# time: 2023-04-17 06:23:43 CEST
# mode: julia
	AAA = rand(SymSecondOrderTensorValue{2,Int})
# time: 2023-04-17 06:23:53 CEST
# mode: julia
	AAA = rand(SymThirdOrderTensorValue{2,Int})
# time: 2023-04-17 06:25:16 CEST
# mode: julia
	AAA = rand(SymFourthOrderTensorValue{2,Int})
# time: 2023-04-17 06:25:24 CEST
# mode: julia
	AAA = rand(SymFourthOrderTensorValue{1,Int})
# time: 2023-04-17 06:25:40 CEST
# mode: julia
	AAA = rand(SymFourthOrderTensorValue{2,Int})
# time: 2023-04-17 06:27:04 CEST
# mode: julia
	tensor = SymFourthOrderTensorValue( (5.0, 9.0, 1.0,
	                                     8.0, 8.0, 2.0, 
	                                     6.0, 7.0, 3.0,) )
# time: 2023-04-17 06:36:13 CEST
# mode: julia
	elementA[1,1]
# time: 2023-04-17 06:36:18 CEST
# mode: julia
	elementA[2,2]
# time: 2023-04-17 06:36:24 CEST
# mode: julia
	elementA[1,2]
# time: 2023-04-17 08:38:25 CEST
# mode: julia
	TensorValue([1.0 2.0; 3.0 4.0])
# time: 2023-04-17 08:38:30 CEST
# mode: julia
	TT = TensorValue([1.0 2.0; 3.0 4.0])
# time: 2023-04-17 08:38:34 CEST
# mode: julia
	TT = TensorValue([1.0 2.0; 3.0 4.0])[:,:]
# time: 2023-04-17 08:38:41 CEST
# mode: julia
	TT = TensorValue([1.0 2.0; 3.0 4.0])[1,1]
# time: 2023-04-17 08:38:47 CEST
# mode: julia
	TT = TensorValue([1.0 2.0; 3.0 4.0])
# time: 2023-04-17 08:38:50 CEST
# mode: julia
	TT[1,1]
# time: 2023-04-17 08:38:55 CEST
# mode: julia
	TT[:,:]
# time: 2023-04-17 08:39:00 CEST
# mode: julia
	TT[end,end]
# time: 2023-04-17 08:39:14 CEST
# mode: julia
	TT[1:2,1:2]
# time: 2023-04-17 09:32:05 CEST
# mode: julia
	u = get_fe_basis(U)
# time: 2023-04-17 09:32:11 CEST
# mode: julia
	tf
# time: 2023-04-17 09:32:40 CEST
# mode: julia
	du = ∇(u)
# time: 2023-04-17 09:33:00 CEST
# mode: julia
	u,θ = get_fe_basis(U)
# time: 2023-04-17 09:33:11 CEST
# mode: julia
	du = ∇(u)
# time: 2023-04-17 09:33:35 CEST
# mode: julia
	dude = du⋅tf
# time: 2023-04-17 09:34:10 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-04-17 09:34:14 CEST
# mode: julia
	dude(pts)
# time: 2023-04-17 09:34:18 CEST
# mode: julia
	dude(pts)[1]
# time: 2023-04-17 09:34:42 CEST
# mode: julia
	u
# time: 2023-04-17 09:34:49 CEST
# mode: julia
	θ
# time: 2023-04-17 09:35:09 CEST
# mode: julia
	u.cell_basis
# time: 2023-04-17 09:35:25 CEST
# mode: julia
	u.single_field
# time: 2023-04-17 09:35:33 CEST
# mode: julia
	u = u.single_field
# time: 2023-04-17 09:35:38 CEST
# mode: julia
	du = ∇(u)
# time: 2023-04-17 09:35:45 CEST
# mode: julia
	dude = du⋅tf
# time: 2023-04-17 09:35:58 CEST
# mode: julia
	dude(pts)[1]
# time: 2023-04-17 09:36:16 CEST
# mode: julia
	du(pts)[1]
# time: 2023-04-17 09:36:43 CEST
# mode: julia
	θ = θ.single_field
# time: 2023-04-17 09:36:54 CEST
# mode: julia
	dθ = ∇(θ)
# time: 2023-04-17 09:37:09 CEST
# mode: julia
	dθde = dθ⋅tf
# time: 2023-04-17 09:37:21 CEST
# mode: julia
	tf
# time: 2023-04-17 09:37:23 CEST
# mode: julia
	df
# time: 2023-04-17 09:37:27 CEST
# mode: julia
	dθ
# time: 2023-04-17 09:37:36 CEST
# mode: julia
	dθde = dθ⋅tf
# time: 2023-04-17 09:37:55 CEST
# mode: julia
	dθ
# time: 2023-04-17 09:37:58 CEST
# mode: julia
	dθ(pts)[1]
# time: 2023-04-17 09:38:49 CEST
# mode: julia
	dθde = tf⋅dθ
# time: 2023-04-17 09:39:18 CEST
# mode: julia
	dθde = lazy_map(Broadcast(⋅),tf,dθ)
# time: 2023-04-17 09:39:24 CEST
# mode: julia
	dθde = lazy_map(Broadcasting(⋅),tf,dθ)
# time: 2023-04-17 09:39:34 CEST
# mode: julia
	dθde = lazy_array(Broadcasting(⋅),tf,dθ)
# time: 2023-04-17 09:39:45 CEST
# mode: julia
	using Gridap.Arrays
# time: 2023-04-17 09:39:48 CEST
# mode: julia
	dθde = lazy_array(Broadcasting(⋅),tf,dθ)
# time: 2023-04-17 09:39:56 CEST
# mode: julia
	dθde = lazy_map(Broadcasting(⋅),tf,dθ)
# time: 2023-04-17 09:40:30 CEST
# mode: julia
	dθde = tf⋅dθ
# time: 2023-04-17 09:41:12 CEST
# mode: julia
	dθde = tf⋅dθ'
# time: 2023-04-17 09:41:24 CEST
# mode: julia
	dθde(pts)[1]
# time: 2023-04-17 09:41:43 CEST
# mode: julia
	dθde = tf⋅dθ'
# time: 2023-04-17 09:42:04 CEST
# mode: julia
	dθde = dθ'⋅tf
# time: 2023-04-17 09:42:19 CEST
# mode: julia
	dθde(pts)[1]
# time: 2023-04-17 09:58:31 CEST
# mode: julia
	dθde = dθ'⋅tf
# time: 2023-04-17 09:58:55 CEST
# mode: julia
	dθde1 = dθd1⋅tfθ
# time: 2023-04-17 09:59:00 CEST
# mode: julia
	dθde1 = dθde⋅tfθ
# time: 2023-04-17 09:59:31 CEST
# mode: julia
	dθde1(pts)[1]]
# time: 2023-04-17 09:59:32 CEST
# mode: julia
	dθde1(pts)[1]
# time: 2023-04-17 10:03:09 CEST
# mode: julia
	elementA[1,1]
# time: 2023-04-17 10:03:13 CEST
# mode: julia
	elementA[1,2]
# time: 2023-04-17 10:03:15 CEST
# mode: julia
	elementA[2,2]
# time: 2023-04-19 04:10:16 CEST
# mode: julia
	element[1,1]
# time: 2023-04-19 04:10:27 CEST
# mode: julia
	element[1,2\]
# time: 2023-04-19 04:10:30 CEST
# mode: julia
	element[1,2]
# time: 2023-04-19 04:10:33 CEST
# mode: julia
	element[1,1]
# time: 2023-04-19 04:10:35 CEST
# mode: julia
	element[1,2]
# time: 2023-04-19 04:10:37 CEST
# mode: julia
	element[2,2]
# time: 2023-04-19 04:16:40 CEST
# mode: julia
	elementS
# time: 2023-04-19 04:18:59 CEST
# mode: julia
	elementS[1,1]]
# time: 2023-04-19 04:19:03 CEST
# mode: julia
	elementS[1,1]
# time: 2023-04-19 04:19:10 CEST
# mode: julia
	elementS[1,2]
# time: 2023-04-19 04:19:21 CEST
# mode: julia
	elementS[2,2]
# time: 2023-05-01 03:08:14 CEST
# mode: julia
	id_face
# time: 2023-05-01 03:08:35 CEST
# mode: julia
	get_face_labeling(model)
# time: 2023-05-01 03:09:29 CEST
# mode: julia
	get_tag_to_name(get_face_labeling(model))
# time: 2023-05-01 03:51:04 CEST
# mode: julia
	get_cell_coordinates(model)
# time: 2023-05-01 03:51:29 CEST
# mode: julia
	c2n_map
# time: 2023-05-01 03:56:00 CEST
# mode: julia
	a < b
# time: 2023-05-01 04:06:13 CEST
# mode: julia
	map(x->x+1,a)
# time: 2023-05-01 04:07:26 CEST
# mode: julia
	ceil(0.1)
# time: 2023-05-01 04:07:36 CEST
# mode: julia
	floor(0.1)
# time: 2023-05-01 04:08:15 CEST
# mode: julia
	com_coords
# time: 2023-05-01 04:10:49 CEST
# mode: julia
	int_coords
# time: 2023-05-01 04:21:36 CEST
# mode: julia
	interface_coords[p]
# time: 2023-05-01 04:26:08 CEST
# mode: julia
	using DataStructures
# time: 2023-05-01 04:26:14 CEST
# mode: julia
	counter
# time: 2023-05-01 04:26:24 CEST
# mode: julia
	counter(interface_coords)
# time: 2023-05-01 04:27:26 CEST
# mode: julia
	counter(map(c->c[3],interface_coords))
# time: 2023-05-01 04:29:08 CEST
# mode: julia
	interface_coords[p]
# time: 2023-05-01 04:36:14 CEST
# mode: julia
	interface_faces
# time: 2023-05-01 04:49:34 CEST
# mode: julia
	interface_coords[p]
# time: 2023-05-01 04:49:53 CEST
# mode: julia
	Table(data,y_ptrs)
# time: 2023-05-01 05:18:45 CEST
# mode: julia
	length(tabl)
# time: 2023-05-01 05:57:03 CEST
# mode: julia
	Vector{Int}(ones, 4)
# time: 2023-05-01 05:57:16 CEST
# mode: julia
	Vector{Int}(zeros, 4)
# time: 2023-05-01 05:57:30 CEST
# mode: julia
	Vector{Int}(0, 4)
# time: 2023-05-01 05:57:46 CEST
# mode: julia
	Vector{Int}(4)
# time: 2023-05-01 05:58:33 CEST
# mode: julia
	Vector{Int}(1, 4)
# time: 2023-05-01 05:59:40 CEST
# mode: julia
	Vector{Int}(4)
# time: 2023-05-01 05:59:44 CEST
# mode: julia
	Vector{Integer}(4)
# time: 2023-05-01 06:00:36 CEST
# mode: julia
	zeros(Int, 4)
# time: 2023-05-01 06:08:02 CEST
# mode: julia
	o2n_cells[1]
# time: 2023-05-01 06:08:04 CEST
# mode: julia
	o2n_cells[2]
# time: 2023-05-01 06:08:07 CEST
# mode: julia
	o2n_cells[2] = 2
# time: 2023-05-01 06:08:10 CEST
# mode: julia
	o2n_cells
# time: 2023-05-01 06:11:11 CEST
# mode: julia
	child_ids
# time: 2023-05-04 03:46:57 CEST
# mode: julia
	exit
# time: 2023-05-04 03:50:47 CEST
# mode: julia
	Γ
# time: 2023-05-04 03:51:16 CEST
# mode: julia
	Γ.glue.face_to_bgface
# time: 2023-05-04 03:52:15 CEST
# mode: julia
	o2n_cells
# time: 2023-05-04 03:52:45 CEST
# mode: julia
	child_ids
# time: 2023-05-04 03:54:13 CEST
# mode: julia
	UnstructuredDiscreteModel(Γ)
# time: 2023-05-04 03:54:24 CEST
# mode: julia
	mm = UnstructuredDiscreteModel(Γ)
# time: 2023-05-04 03:54:37 CEST
# mode: julia
	mm.grid_topology.polytopes
# time: 2023-05-04 03:54:47 CEST
# mode: julia
	mm.grid_topology.vertex_coordinates
# time: 2023-05-04 03:55:24 CEST
# mode: julia
	mm.grid
# time: 2023-05-04 03:55:32 CEST
# mode: julia
	mm.grid.node_coordinates
# time: 2023-05-04 03:55:39 CEST
# mode: julia
	mm.grid.reffes
# time: 2023-05-04 03:56:05 CEST
# mode: julia
	mm.grid_topology.n_m_to_nface_to_mfaces
# time: 2023-05-04 03:56:17 CEST
# mode: julia
	get_faces(model)
# time: 2023-05-04 03:56:24 CEST
# mode: julia
	num_faces(model)
# time: 2023-05-04 03:56:28 CEST
# mode: julia
	num_faces(model,Val(1))
# time: 2023-05-04 03:56:34 CEST
# mode: julia
	num_faces(model,1)
# time: 2023-05-04 03:56:36 CEST
# mode: julia
	num_faces(model,0)
# time: 2023-05-04 03:56:38 CEST
# mode: julia
	num_faces(model,2)
# time: 2023-05-04 03:56:54 CEST
# mode: julia
	num_faces(model,3)
# time: 2023-05-04 03:57:11 CEST
# mode: julia
	mm
# time: 2023-05-04 03:57:23 CEST
# mode: julia
	mm.grid_topology.polytopes
# time: 2023-05-04 03:57:33 CEST
# mode: julia
	num_faces(mm,1)
# time: 2023-05-04 03:57:36 CEST
# mode: julia
	num_faces(mm,0)
# time: 2023-05-04 03:57:41 CEST
# mode: julia
	num_faces(mm,2)
# time: 2023-05-04 03:57:44 CEST
# mode: julia
	num_faces(mm,3)
# time: 2023-05-04 03:57:48 CEST
# mode: julia
	num_faces(mm,1)
# time: 2023-05-04 04:02:13 CEST
# mode: julia
	rr.ref_grid
# time: 2023-05-04 04:02:58 CEST
# mode: julia
	get_faces(rr.ref_grid,2,0)
# time: 2023-05-04 04:03:07 CEST
# mode: julia
	Gridap.Geometry.get_faces(rr.ref_grid,2,0)
# time: 2023-05-04 04:03:25 CEST
# mode: julia
	Gridap.Geometry.get_faces(rr.ref_grid.grid_topology,2,0)
# time: 2023-05-04 04:04:24 CEST
# mode: julia
	get_cell_map(rr.ref_grid)
# time: 2023-05-04 04:05:26 CEST
# mode: julia
	AdaptivityGlue
# time: 2023-05-04 04:08:23 CEST
# mode: pkg
	add FillArrays
# time: 2023-05-04 04:12:23 CEST
# mode: pkg
	status
# time: 2023-05-04 04:13:48 CEST
# mode: pkg
	compat
# time: 2023-05-04 04:14:24 CEST
# mode: pkg
	up
# time: 2023-05-04 04:15:52 CEST
# mode: pkg
	compat
# time: 2023-05-04 04:16:19 CEST
# mode: pkg
	status
# time: 2023-05-04 04:18:54 CEST
# mode: julia
	glue.o2n_faces_map
# time: 2023-05-04 04:19:09 CEST
# mode: julia
	tabl
# time: 2023-05-04 04:28:46 CEST
# mode: julia
	cell_dims(coarse_face_model)
# time: 2023-05-04 04:28:49 CEST
# mode: julia
	cell_dim(coarse_face_model)
# time: 2023-05-04 04:29:50 CEST
# mode: julia
	coarse_face_model.grid_topology
# time: 2023-05-04 04:29:55 CEST
# mode: julia
	coarse_face_model.grid_topology.cell_type
# time: 2023-05-04 04:30:00 CEST
# mode: julia
	coarse_face_model.grid_topology.polytopes
# time: 2023-05-04 04:30:08 CEST
# mode: julia
	coarse_face_model.grid_topology.vertex_coordinates
# time: 2023-05-04 04:34:13 CEST
# mode: julia
	mm.face_labeling.tag_to_entities
# time: 2023-05-04 04:34:24 CEST
# mode: julia
	mm.face_labeling.tag_to_name
# time: 2023-05-04 04:53:33 CEST
# mode: julia
	λf
# time: 2023-05-04 04:53:59 CEST
# mode: julia
	CellField()
# time: 2023-05-04 04:54:12 CEST
# mode: julia
	CellField(λf,Γ)
# time: 2023-05-04 04:54:16 CEST
# mode: julia
	CellField(λ,Γ)
# time: 2023-05-04 04:56:44 CEST
# mode: julia
	CellField(λ,Γc)
# time: 2023-05-04 04:57:35 CEST
# mode: julia
	GenericCellField(λf,Γ)
# time: 2023-05-04 04:57:42 CEST
# mode: julia
	Gridap.CellData.GenericCellField(λf,Γ)
# time: 2023-05-04 04:58:02 CEST
# mode: julia
	Gridap.CellData.GenericCellField(λf,Γ,Gridap.CellData.ReferenceDomain())
# time: 2023-05-04 05:01:21 CEST
# mode: julia
	contr.dict
# time: 2023-05-04 05:02:17 CEST
# mode: julia
	last(contr.dict)
# time: 2023-05-04 05:02:20 CEST
# mode: julia
	lastindex(contr.dict)
# time: 2023-05-04 05:02:27 CEST
# mode: julia
	first(contr.dict)
# time: 2023-05-04 05:02:38 CEST
# mode: julia
	contr.dict
# time: 2023-05-04 05:02:41 CEST
# mode: julia
	contr.dict[2]
# time: 2023-05-04 05:03:08 CEST
# mode: julia
	contr.dict[keys(contr.dict)[2]]
# time: 2023-05-04 05:03:17 CEST
# mode: julia
	keys(contr.dict)
# time: 2023-05-04 05:03:21 CEST
# mode: julia
	keys(contr.dict)[2]
# time: 2023-05-04 05:03:37 CEST
# mode: julia
	keys(contr.dict[BoundaryTriangulation()])
# time: 2023-05-11 07:47:32 CEST
# mode: julia
	aaa = ones(Int,4)
# time: 2023-05-11 07:47:57 CEST
# mode: julia
	append!(aaa, ones(Int, 2))
# time: 2023-05-11 07:48:30 CEST
# mode: julia
	ones(Int, 2)*3
# time: 2023-05-15 04:29:13 CEST
# mode: pkg
	status
# time: 2023-05-15 04:29:45 CEST
# mode: pkg
	add Gridap#glued-triangulations
# time: 2023-05-15 04:31:55 CEST
# mode: julia
	using Gridap
# time: 2023-05-15 04:32:09 CEST
# mode: julia
	using Gridap.Adaptivity
# time: 2023-05-15 04:32:19 CEST
# mode: julia
	Adaptivity.GluedTriangulation
# time: 2023-05-16 01:42:27 CEST
# mode: julia
	aa = tr_Γf(μ)
# time: 2023-05-16 02:29:08 CEST
# mode: pkg
	up
# time: 2023-05-16 02:31:14 CEST
# mode: pkg
	status
# time: 2023-05-16 02:31:35 CEST
# mode: julia
	using GridapGmesh
# time: 2023-05-16 02:31:48 CEST
# mode: julia
	using GridapGmsh
# time: 2023-05-16 02:33:19 CEST
# mode: julia
	upgrade_manifest
# time: 2023-05-16 02:33:26 CEST
# mode: julia
	Pkg().upgrade_manifest
# time: 2023-05-16 02:33:34 CEST
# mode: pkg
	upgrade_manifest
# time: 2023-05-16 02:33:43 CEST
# mode: pkg
	upgrade manifest
# time: 2023-05-16 02:34:19 CEST
# mode: julia
	Pkg.upgrade_manifest()
# time: 2023-05-16 02:34:23 CEST
# mode: julia
	using Pkg
# time: 2023-05-16 02:34:24 CEST
# mode: julia
	Pkg.upgrade_manifest()
# time: 2023-05-16 02:35:18 CEST
# mode: julia
	exit()
# time: 2023-05-16 02:36:18 CEST
# mode: pkg
	instantiate
# time: 2023-05-16 02:36:42 CEST
# mode: pkg
	add Gridap#glued-triangulations
# time: 2023-05-16 02:36:58 CEST
# mode: pkg
	precompile
# time: 2023-05-16 02:37:07 CEST
# mode: pkg
	up
# time: 2023-05-16 02:37:14 CEST
# mode: pkg
	build
# time: 2023-05-16 02:38:45 CEST
# mode: julia
	using GridapGmsh
# time: 2023-05-16 02:40:47 CEST
# mode: pkg
	status
# time: 2023-05-16 02:40:58 CEST
# mode: pkg
	up
# time: 2023-05-16 02:41:08 CEST
# mode: pkg
	upgrade
# time: 2023-05-16 02:41:10 CEST
# mode: pkg
	update
# time: 2023-05-16 02:41:41 CEST
# mode: pkg
	rm GridapGmsh
# time: 2023-05-16 02:41:53 CEST
# mode: pkg
	up
# time: 2023-05-16 02:42:05 CEST
# mode: pkg
	add GridapGmsh
# time: 2023-05-16 02:44:34 CEST
# mode: pkg
	add GridapGmsh#v0.6.1
# time: 2023-05-16 02:45:28 CEST
# mode: pkg
	compat
# time: 2023-05-16 02:45:50 CEST
# mode: pkg
	up
# time: 2023-05-16 02:45:59 CEST
# mode: pkg
	status
# time: 2023-05-16 02:46:26 CEST
# mode: pkg
	add GridapDistributed
# time: 2023-05-16 02:47:24 CEST
# mode: pkg
	status
# time: 2023-05-16 02:56:19 CEST
# mode: pkg
	exit()
# time: 2023-05-16 02:56:27 CEST
# mode: julia
	ecit()
# time: 2023-05-16 02:56:29 CEST
# mode: julia
	exit()
# time: 2023-05-16 03:13:06 CEST
# mode: pkg
	status
# time: 2023-05-16 03:27:49 CEST
# mode: pkg
	add GridapGmsh#fix-gridap-typo
# time: 2023-05-16 04:41:12 CEST
# mode: julia
	56⁵
# time: 2023-05-16 04:41:17 CEST
# mode: julia
	56^5
# time: 2023-05-16 04:41:50 CEST
# mode: julia
	a = Int(5)
# time: 2023-05-16 04:41:54 CEST
# mode: julia
	aaa = Int(5)
# time: 2023-05-16 04:42:00 CEST
# mode: julia
	45^aaa
# time: 2023-05-17 04:27:13 CEST
# mode: julia
	get_cell_points(Γ)
# time: 2023-05-17 04:27:23 CEST
# mode: julia
	pts = get_cell_points(Γ)
# time: 2023-05-17 04:27:26 CEST
# mode: julia
	pts
# time: 2023-05-17 07:47:51 CEST
# mode: julia
	Λe
# time: 2023-05-17 07:48:03 CEST
# mode: julia
	num_dims(Λe)
# time: 2023-05-17 08:02:56 CEST
# mode: julia
	get_zero_free_dof_values(Ue)
# time: 2023-05-17 08:02:59 CEST
# mode: julia
	get_zero_free_values(Ue)
# time: 2023-05-17 08:03:05 CEST
# mode: julia
	using Gridap.FESpaces
# time: 2023-05-17 08:03:06 CEST
# mode: julia
	get_zero_free_values(Ue)
# time: 2023-05-17 08:03:11 CEST
# mode: julia
	zero_free_values(Ue)
# time: 2023-05-17 08:03:55 CEST
# mode: julia
	num_free_dofs(Ue)
# time: 2023-05-17 08:06:07 CEST
# mode: julia
	sol
# time: 2023-05-17 08:07:01 CEST
# mode: julia
	soluh
# time: 2023-05-17 08:07:02 CEST
# mode: julia
	uh
# time: 2023-05-18 02:38:44 CEST
# mode: pkg
	status
# time: 2023-05-19 04:37:31 CEST
# mode: julia
	SingleFieldFE
# time: 2023-05-19 04:37:38 CEST
# mode: julia
	SingleFieldFESpace
# time: 2023-05-19 06:02:05 CEST
# mode: julia
	contrA
# time: 2023-05-19 06:05:05 CEST
# mode: julia
	contrA[1,1]]
# time: 2023-05-19 06:05:07 CEST
# mode: julia
	contrA[1,1]
# time: 2023-05-19 06:05:12 CEST
# mode: julia
	contrA
# time: 2023-05-19 06:05:19 CEST
# mode: julia
	elementA[1,1]
# time: 2023-05-19 06:05:22 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:07:09 CEST
# mode: julia
	contrA
# time: 2023-05-19 06:07:10 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:07:51 CEST
# mode: julia
	contrA
# time: 2023-05-19 06:07:52 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:09:43 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:10:51 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:10:53 CEST
# mode: julia
	elementA[1,1]
# time: 2023-05-19 06:16:38 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:17:41 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:18:25 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:19:23 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:19:52 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:20:28 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:21:25 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:21:57 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:22:15 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:22:36 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:22:38 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:23:40 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:24:20 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:26:38 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:27:16 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:28:12 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:37:19 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:37:24 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 06:39:36 CEST
# mode: julia
	A
# time: 2023-05-19 06:39:51 CEST
# mode: julia
	contrA.dict
# time: 2023-05-19 06:41:03 CEST
# mode: julia
	contrl[1]
# time: 2023-05-19 06:41:13 CEST
# mode: julia
	elementl[1]
# time: 2023-05-19 06:41:39 CEST
# mode: julia
	elementl[2]
# time: 2023-05-19 06:43:21 CEST
# mode: julia
	elementl[1]
# time: 2023-05-19 06:43:22 CEST
# mode: julia
	elementl[2]
# time: 2023-05-19 06:46:32 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 06:46:40 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 07:19:39 CEST
# mode: pkg
	status
# time: 2023-05-19 07:20:16 CEST
# mode: julia
	add Gridap#master
# time: 2023-05-19 07:20:25 CEST
# mode: pkg
	status
# time: 2023-05-19 07:20:36 CEST
# mode: pkg
	add Gridap#master
# time: 2023-05-19 07:30:05 CEST
# mode: pkg
	status
# time: 2023-05-19 07:37:46 CEST
# mode: julia
	elementA[1,2]
# time: 2023-05-19 07:37:49 CEST
# mode: pkg
	status
# time: 2023-05-19 07:37:59 CEST
# mode: julia
	elementA[2,1]
# time: 2023-05-19 08:12:07 CEST
# mode: pkg
	status
# time: 2023-05-19 08:12:33 CEST
# mode: pkg
	add Gridap#glued-triangulations
# time: 2023-05-19 08:14:27 CEST
# mode: pkg
	up
# time: 2023-05-19 08:15:48 CEST
# mode: pkg
	status
# time: 2023-05-19 08:31:57 CEST
# mode: julia
	xe_X0
# time: 2023-05-22 05:57:08 CEST
# mode: julia
	TensorValue(1[get_i∘(1,u),get_i∘(1,u),get_i∘(2,u)])
# time: 2023-05-22 05:57:16 CEST
# mode: julia
	TensorValue(1)
# time: 2023-05-22 05:57:20 CEST
# mode: julia
	TensorValue(1,2)
# time: 2023-05-22 05:57:22 CEST
# mode: julia
	TensorValue(1 2)
# time: 2023-05-22 05:57:26 CEST
# mode: julia
	TensorValue(1,2)
# time: 2023-05-22 05:57:31 CEST
# mode: julia
	TensorValue([1 2])
# time: 2023-05-22 05:59:18 CEST
# mode: julia
	TensorValue([3 2])
# time: 2023-05-22 05:59:24 CEST
# mode: julia
	A = TensorValue([1 2])
# time: 2023-05-22 05:59:28 CEST
# mode: julia
	B = TensorValue([3 2])
# time: 2023-05-22 05:59:31 CEST
# mode: julia
	A⋅B
# time: 2023-05-22 05:59:37 CEST
# mode: julia
	A*B
# time: 2023-05-22 05:59:52 CEST
# mode: julia
	A⊙ B
# time: 2023-05-22 05:59:53 CEST
# mode: julia
	A⊙B
# time: 2023-05-23 01:40:35 CEST
# mode: julia
	aΓa(Uα,Vα)
# time: 2023-05-23 01:44:52 CEST
# mode: julia
	Uα
# time: 2023-05-23 01:45:02 CEST
# mode: julia
	Tuple(Uα...)
# time: 2023-05-23 01:45:13 CEST
# mode: julia
	[Uα...]
# time: 2023-05-23 02:01:40 CEST
# mode: julia
	u
# time: 2023-05-23 02:09:01 CEST
# mode: julia
	f1(intrfA)
# time: 2023-05-23 02:10:17 CEST
# mode: julia
	E1,E2,E3 = f1(intrfA)
# time: 2023-05-23 02:10:32 CEST
# mode: julia
	CellData.get_data(E1)
# time: 2023-05-23 02:10:37 CEST
# mode: julia
	Gridap.CellData.get_data(E1)
# time: 2023-05-23 02:11:16 CEST
# mode: julia
	get_cell_points(Ω)
# time: 2023-05-23 02:11:20 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-05-23 02:11:23 CEST
# mode: julia
	E1(pts)
# time: 2023-05-23 02:11:27 CEST
# mode: julia
	E1(pts)[1]
# time: 2023-05-23 02:11:55 CEST
# mode: julia
	pts = get_cell_points(Γ)
# time: 2023-05-23 02:11:59 CEST
# mode: julia
	pts = get_cell_points(Γa)
# time: 2023-05-23 02:12:02 CEST
# mode: julia
	E1(pts)[1]
# time: 2023-05-23 02:12:32 CEST
# mode: julia
	Gridap.TensorValue()
# time: 2023-05-23 02:13:07 CEST
# mode: julia
	Gridap.TensorValue(zeros(2,3))
# time: 2023-05-23 02:13:20 CEST
# mode: julia
	Gridap.TensorValue(zeros(3,2))
# time: 2023-05-23 02:13:56 CEST
# mode: julia
	Gridap.TensorValue(3,2,[1,2,3,4,5,6])
# time: 2023-05-23 02:14:01 CEST
# mode: julia
	Gridap.TensorValue((3,2),[1,2,3,4,5,6])
# time: 2023-05-23 02:14:07 CEST
# mode: julia
	Gridap.TensorValue((3,2),Tuple([1,2,3,4,5,6]))
# time: 2023-05-23 02:15:32 CEST
# mode: julia
	Gridap.TensorValue{3,2}(Tuple([1,2,3,4,5,6]))
# time: 2023-05-23 02:15:42 CEST
# mode: julia
	Gridap.TensorValue{3,2}([1,2,3,4,5,6])
# time: 2023-05-23 02:15:50 CEST
# mode: julia
	Gridap.TensorValue{3,2}((1,2,3,4,5,6))
# time: 2023-05-23 02:16:04 CEST
# mode: julia
	x = Gridap.TensorValue{3,2}((1,2,3,4,5,6))
# time: 2023-05-23 02:16:19 CEST
# mode: julia
	to_array(x)
# time: 2023-05-23 02:16:30 CEST
# mode: julia
	Gridap.TensorValues.to_array(x)
# time: 2023-05-23 02:16:52 CEST
# mode: julia
	Gridap.TensorValues.get_array(x)
# time: 2023-05-23 02:18:20 CEST
# mode: julia
	E = f1(intrfA)
# time: 2023-05-23 02:18:36 CEST
# mode: julia
	E(pts)[1]
# time: 2023-05-23 02:19:09 CEST
# mode: julia
	E = f1(intrfA)
# time: 2023-05-23 02:19:10 CEST
# mode: julia
	E(pts)[1]
# time: 2023-05-23 02:20:35 CEST
# mode: julia
	λ⋅E⋅u
# time: 2023-05-23 02:20:44 CEST
# mode: julia
	res = λ⋅E⋅u
# time: 2023-05-23 02:20:49 CEST
# mode: julia
	res(pts)[1]
# time: 2023-05-23 02:22:43 CEST
# mode: julia
	res = μ⋅E⋅u
# time: 2023-05-23 02:22:48 CEST
# mode: julia
	res(pts)[1]
# time: 2023-05-23 03:09:11 CEST
# mode: julia
	intrfA
# time: 2023-05-23 03:09:19 CEST
# mode: julia
	intrfA.Ef
# time: 2023-05-23 03:09:23 CEST
# mode: julia
	intrfA.Ef[1]
# time: 2023-05-23 03:09:31 CEST
# mode: julia
	intrfA.Ef(pts)[1]
# time: 2023-05-23 03:09:48 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-05-23 03:09:51 CEST
# mode: julia
	intrfA.Ef(pts)[1]
# time: 2023-05-23 03:10:03 CEST
# mode: julia
	intrfA.zf(pts)[1]
# time: 2023-05-23 03:10:18 CEST
# mode: julia
	pts = get_cell_points(Γ)
# time: 2023-05-23 03:10:30 CEST
# mode: julia
	ptsΓa = get_cell_points(Γa)
# time: 2023-05-23 03:10:51 CEST
# mode: julia
	E1(ptsΓa)[1]
# time: 2023-05-23 03:10:58 CEST
# mode: julia
	intrfA.E1(ptsΓa)[1]
# time: 2023-05-23 03:11:04 CEST
# mode: julia
	intrfA.Ef(ptsΓa)[1]
# time: 2023-05-23 03:21:02 CEST
# mode: julia
	print_info(intrfA)
# time: 2023-05-23 03:21:22 CEST
# mode: julia
	intrfA.Ef
# time: 2023-05-23 03:21:33 CEST
# mode: julia
	intrfA.Ef(pts)
# time: 2023-05-23 03:21:42 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-05-23 03:21:46 CEST
# mode: julia
	intrfA.Ef(pts)[1]
# time: 2023-05-23 03:24:20 CEST
# mode: julia
	element[2,1]
# time: 2023-05-23 03:44:16 CEST
# mode: julia
	E1,E2,E3 = f1(intrfA)
# time: 2023-05-23 03:45:50 CEST
# mode: julia
	E = f1(intrfA)
# time: 2023-05-23 03:45:59 CEST
# mode: julia
	E(pts)[1]
# time: 2023-05-23 03:46:05 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-05-23 03:46:07 CEST
# mode: julia
	E(pts)[1]
# time: 2023-05-23 03:47:12 CEST
# mode: julia
	E = f1(intrfA)
# time: 2023-05-23 03:47:43 CEST
# mode: julia
	E(pts)[1]
# time: 2023-05-23 03:48:30 CEST
# mode: julia
	print_info(intrfA)
# time: 2023-05-23 03:49:36 CEST
# mode: julia
	ptsΓa = get_cell_points(Γa)
# time: 2023-05-23 03:49:54 CEST
# mode: julia
	Ef(pts)[:]
# time: 2023-05-23 03:50:13 CEST
# mode: julia
	Ef(ptsΓa)[:]
# time: 2023-05-23 06:00:53 CEST
# mode: julia
	Set(Integer[72, 89])
# time: 2023-05-23 06:00:58 CEST
# mode: julia
	a = Set(Integer[72, 89])
# time: 2023-05-23 06:01:07 CEST
# mode: julia
	b = Set(Integer[89, 72])
# time: 2023-05-23 06:01:10 CEST
# mode: julia
	a == b
# time: 2023-05-23 06:27:38 CEST
# mode: julia
	print_info(intrfA)
# time: 2023-05-23 06:28:34 CEST
# mode: julia
	CTf[1]
# time: 2023-05-23 06:28:44 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-05-23 06:28:54 CEST
# mode: julia
	CTf[1](pts)[1]
# time: 2023-05-23 06:30:40 CEST
# mode: julia
	CTf[1](pts)[5]
# time: 2023-05-23 06:30:43 CEST
# mode: julia
	CTf[1](pts)[50]
# time: 2023-05-23 06:30:47 CEST
# mode: julia
	CTf[1](pts)[10]
# time: 2023-05-23 06:30:50 CEST
# mode: julia
	CTf[1](pts)[80]
# time: 2023-05-23 06:35:38 CEST
# mode: julia
	CTf[1](pts)
# time: 2023-05-23 06:36:30 CEST
# mode: julia
	for ct in CT(pts)
	  println(ct)
	end
# time: 2023-05-23 06:36:41 CEST
# mode: julia
	for ct in CTf(pts)
	  println(ct)
	end
# time: 2023-05-23 06:36:56 CEST
# mode: julia
	for ct in [CTf(pts)]
	  println(ct)
	end
# time: 2023-05-23 06:37:02 CEST
# mode: julia
	for ct in [CTf(pts)...]
	  println(ct)
	end
# time: 2023-05-23 06:37:19 CEST
# mode: julia
	CTf[1](pts)
# time: 2023-05-23 06:57:31 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-05-23 06:57:37 CEST
# mode: julia
	Ef(pts)
# time: 2023-05-23 06:58:14 CEST
# mode: julia
	Ef(pts)[80]
# time: 2023-05-23 06:58:17 CEST
# mode: julia
	Ef(pts)[81]
# time: 2023-05-23 06:58:32 CEST
# mode: julia
	Ef(pts)[160]
# time: 2023-05-23 06:58:34 CEST
# mode: julia
	Ef(pts)[161]
# time: 2023-05-23 07:04:18 CEST
# mode: julia
	tags
# time: 2023-05-23 07:05:05 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-05-23 07:05:09 CEST
# mode: julia
	Ef(pts)[80]
# time: 2023-05-23 07:05:15 CEST
# mode: julia
	Ef(pts)[160]
# time: 2023-05-23 07:05:20 CEST
# mode: julia
	Ef(pts)[240]
# time: 2023-05-23 07:11:25 CEST
# mode: julia
	Ef(pts)[80]
# time: 2023-05-23 07:11:42 CEST
# mode: julia
	Ef(pts)[160]
# time: 2023-05-23 07:11:49 CEST
# mode: julia
	Ef(pts)[240]
# time: 2023-05-23 07:14:33 CEST
# mode: julia
	Gridap.CellData.change_domain(Ef,Γa)
# time: 2023-05-23 07:14:54 CEST
# mode: julia
	_Ef = Gridap.CellData.change_domain(Ef,Γa,DomainStyle(Ef))
# time: 2023-05-23 07:15:04 CEST
# mode: julia
	_Ef.cell_field[1]
# time: 2023-05-23 07:15:26 CEST
# mode: julia
	pts = get_cell_points(Γa)
# time: 2023-05-23 07:15:29 CEST
# mode: julia
	pts = get_cell_points(dΓa)
# time: 2023-05-23 07:15:40 CEST
# mode: julia
	pts.cell_phys_point
# time: 2023-05-23 07:16:22 CEST
# mode: julia
	_Ef(pts)
# time: 2023-05-23 07:17:02 CEST
# mode: julia
	_pts = get_cell_points(Ω)
# time: 2023-05-23 07:17:11 CEST
# mode: julia
	_pts.cell_phys_point
# time: 2023-05-23 07:18:56 CEST
# mode: julia
	Ef
# time: 2023-05-23 07:19:10 CEST
# mode: julia
	Ef(VectorValue(0.0,1.0))
# time: 2023-05-23 07:19:29 CEST
# mode: julia
	Ef(VectorValue(0.0,-5))
# time: 2023-05-23 07:19:32 CEST
# mode: julia
	Ef(VectorValue(0.0,-4))
# time: 2023-05-23 07:19:34 CEST
# mode: julia
	Ef(VectorValue(0.0,-3))
# time: 2023-05-23 07:19:37 CEST
# mode: julia
	Ef(VectorValue(0.0,-2))
# time: 2023-05-23 07:19:41 CEST
# mode: julia
	Ef(VectorValue(0.0,-1))
# time: 2023-05-23 07:22:22 CEST
# mode: julia
	Ef
# time: 2023-05-23 07:23:37 CEST
# mode: julia
	_Ef
# time: 2023-05-23 07:24:03 CEST
# mode: julia
	Γa.glue
# time: 2023-05-23 07:24:12 CEST
# mode: julia
	Γa.glue.face_to_bgface
# time: 2023-05-23 07:25:24 CEST
# mode: julia
	Γa.glue.face_to_cell
# time: 2023-05-23 07:25:34 CEST
# mode: julia
	labels
# time: 2023-05-23 07:25:52 CEST
# mode: julia
	labels.d_to_dface_to_entity[3]
# time: 2023-05-23 07:26:14 CEST
# mode: julia
	labels.d_to_dface_to_entity[3][Γa.glue.face_to_cell]
# time: 2023-05-23 07:27:24 CEST
# mode: julia
	labels.tag_to_name[labels.d_to_dface_to_entity[3][Γa.glue.face_to_cell]]
# time: 2023-05-23 07:27:48 CEST
# mode: julia
	labels.d_to_dface_to_entity[3][Γa.glue.face_to_cell]
# time: 2023-05-23 07:27:55 CEST
# mode: julia
	labels.tag_to_entities
# time: 2023-05-23 07:28:05 CEST
# mode: julia
	labels.tag_to_name
# time: 2023-05-23 07:28:59 CEST
# mode: julia
	sum(∫(Ef)*dΓa)
# time: 2023-05-23 07:29:39 CEST
# mode: julia
	sum(∫(Ef)*dΓb)
# time: 2023-05-23 07:29:57 CEST
# mode: julia
	sum(∫(Ef)*Measure(Γb,3))
# time: 2023-05-23 07:30:07 CEST
# mode: julia
	sum(∫(Ef)*Measure(Ω,3))
# time: 2023-05-23 07:31:26 CEST
# mode: julia
	sum(∫(Ef)*Measure(Γa,3))
# time: 2023-05-23 07:32:10 CEST
# mode: julia
	CT1
# time: 2023-05-23 07:32:14 CEST
# mode: julia
	CT2
# time: 2023-05-23 07:32:16 CEST
# mode: julia
	CT3
# time: 2023-05-23 07:33:58 CEST
# mode: julia
	(1000*3 + 100*4 + 300*3)*20
# time: 2023-05-23 07:34:38 CEST
# mode: julia
	Γa.glue.face_to_cell
# time: 2023-05-23 07:35:04 CEST
# mode: julia
	Ef.cell_field[Γa.glue.face_to_cell]
# time: 2023-05-23 07:35:43 CEST
# mode: julia
	lazy_map(evaluate,Ef.cell_field[Γa.glue.face_to_cell],Fill(VectorValue(0.0,0.0),12))
# time: 2023-05-23 07:35:50 CEST
# mode: julia
	lazy_map(evaluate,Ef.cell_field[Γa.glue.face_to_cell],fill(VectorValue(0.0,0.0),12))
# time: 2023-05-23 07:37:08 CEST
# mode: julia
	_Ef = CellField(Ef.cell_field,Ω,PhysicalDomain())
# time: 2023-05-23 07:37:42 CEST
# mode: julia
	_Ef = CellField(Ef,Ω,PhysicalDomain())
# time: 2023-05-23 07:38:01 CEST
# mode: julia
	_Ef = GenericCellField(Ef.cell_field,Ω,PhysicalDomain())
# time: 2023-05-23 07:38:12 CEST
# mode: julia
	_Ef = Gridap.CellData.GenericCellField(Ef.cell_field,Ω,PhysicalDomain())
# time: 2023-05-23 07:38:33 CEST
# mode: julia
	_Ef(pts)
# time: 2023-05-23 07:42:19 CEST
# mode: julia
	Ef
# time: 2023-05-23 07:42:24 CEST
# mode: julia
	Ef.cell_field
# time: 2023-05-24 04:52:18 CEST
# mode: julia
	typeof()
# time: 2023-05-24 04:52:24 CEST
# mode: julia
	typeof(Γ_X0)
# time: 2023-05-24 06:51:33 CEST
# mode: julia
	typeof(ga)
# time: 2023-05-26 04:53:51 CEST
# mode: julia
	intrfA.dΓ
# time: 2023-05-26 04:53:56 CEST
# mode: julia
	intrfA.Γ
# time: 2023-05-26 07:21:51 CEST
# mode: julia
	true
# time: 2023-05-26 07:21:56 CEST
# mode: julia
	typeof(true)
# time: 2023-05-26 07:24:05 CEST
# mode: julia
	c2f_faces
# time: 2023-05-26 07:24:18 CEST
# mode: julia
	reverse(c2f_faces)
# time: 2023-05-26 07:25:42 CEST
# mode: julia
	reverse!(c2f_faces)
# time: 2023-05-26 07:26:01 CEST
# mode: julia
	c2f_faces = reverse(c2f_faces)
# time: 2023-05-29 01:50:13 CEST
# mode: julia
	uF(0)
# time: 2023-05-29 01:50:28 CEST
# mode: julia
	uF(0.0)
# time: 2023-05-29 01:50:36 CEST
# mode: julia
	uF(Point(0))
# time: 2023-05-29 01:50:39 CEST
# mode: julia
	uF(Point(2))
# time: 2023-05-29 01:50:42 CEST
# mode: julia
	uF(Point(1))
# time: 2023-05-29 01:50:46 CEST
# mode: julia
	uF(Point(0.99))
# time: 2023-05-29 01:50:56 CEST
# mode: julia
	uF(Point(0.1))
# time: 2023-05-29 01:50:58 CEST
# mode: julia
	uF(Point(0.2))
# time: 2023-05-29 01:51:00 CEST
# mode: julia
	uF(Point(0.3))
# time: 2023-05-29 01:52:03 CEST
# mode: julia
	uF(Point(0.9))
# time: 2023-05-29 01:52:08 CEST
# mode: julia
	uF(Point(0.5))
# time: 2023-05-29 01:52:11 CEST
# mode: julia
	uF(Point(0.6))
# time: 2023-05-29 01:52:14 CEST
# mode: julia
	uF(Point(0.7))
# time: 2023-05-29 02:09:00 CEST
# mode: pkg
	status
# time: 2023-05-29 02:21:19 CEST
# mode: pkg
	add GridapGmsh>multidimensional
# time: 2023-05-29 02:36:17 CEST
# mode: pkg
	remove GridapGmsh
# time: 2023-05-29 02:36:26 CEST
# mode: pkg
	add /home/deck/Documents/GridapGmsh.jl
# time: 2023-05-29 02:37:36 CEST
# mode: pkg
	up
# time: 2023-05-30 06:28:29 CEST
# mode: julia
	labeling
# time: 2023-05-30 06:28:42 CEST
# mode: julia
	get_face_mask(labeling,tags,0)
# time: 2023-05-30 06:28:57 CEST
# mode: julia
	mask = get_face_mask(labeling,tags,0)
# time: 2023-05-30 06:29:16 CEST
# mode: julia
	node_ids = findall(mask)
# time: 2023-05-30 06:29:31 CEST
# mode: julia
	mask = get_face_mask(labeling,tags,1)
# time: 2023-05-30 06:29:35 CEST
# mode: julia
	node_ids = findall(mask)
# time: 2023-05-30 06:30:20 CEST
# mode: julia
	get_node_coordinates(model)
# time: 2023-05-30 06:30:25 CEST
# mode: julia
	get_vertex_coordinates(model)
# time: 2023-05-30 06:30:34 CEST
# mode: julia
	get_nodes_coordinates(model)
# time: 2023-05-30 06:30:56 CEST
# mode: julia
	get_node_coordinates(model)
# time: 2023-05-30 06:31:03 CEST
# mode: julia
	using Gridap.Geometry
# time: 2023-05-30 06:31:04 CEST
# mode: julia
	get_node_coordinates(model)
# time: 2023-05-30 06:31:10 CEST
# mode: julia
	Geometry.get_node_coordinates(model)
# time: 2023-05-30 06:31:27 CEST
# mode: julia
	coordinates = Geometry.get_node_coordinates(model)
# time: 2023-05-30 06:31:35 CEST
# mode: julia
	mask = get_face_mask(labeling,tags,0)
# time: 2023-05-30 06:31:50 CEST
# mode: julia
	coordinates[nodes]
# time: 2023-05-30 06:31:56 CEST
# mode: julia
	coordinates[mask]
# time: 2023-05-30 06:32:10 CEST
# mode: julia
	nodes = findall(mask)
# time: 2023-05-30 06:32:33 CEST
# mode: julia
	lazy_map(Reindex(coordinates),nodes)
# time: 2023-05-30 06:33:19 CEST
# mode: julia
	coordinates[nodes]
# time: 2023-05-30 06:34:05 CEST
# mode: julia
	mask = get_face_mask(labeling,tags,0)
# time: 2023-05-30 06:36:41 CEST
# mode: julia
	get_faces()
# time: 2023-05-30 06:37:03 CEST
# mode: julia
	Geometry.get_faces(topo,1,0)
# time: 2023-05-31 02:25:30 CEST
# mode: julia
	face_B1_pos
# time: 2023-05-31 04:34:12 CEST
# mode: julia
	ct1
# time: 2023-05-31 04:34:35 CEST
# mode: julia
	ct1.toarray()
# time: 2023-05-31 04:34:40 CEST
# mode: julia
	toarray(ct1)
# time: 2023-05-31 04:34:46 CEST
# mode: julia
	to_array(ct1)
# time: 2023-05-31 04:34:50 CEST
# mode: julia
	array(ct1)
# time: 2023-05-31 04:34:52 CEST
# mode: julia
	Array(ct1)
# time: 2023-05-31 04:35:00 CEST
# mode: julia
	arr = Array(ct1)
# time: 2023-05-31 04:35:01 CEST
# mode: julia
	arr
# time: 2023-05-31 04:35:31 CEST
# mode: julia
	arr = get_array(ct1)
# time: 2023-05-31 04:35:35 CEST
# mode: julia
	arr
# time: 2023-05-31 04:35:44 CEST
# mode: julia
	arr[1,1]
# time: 2023-05-31 04:37:48 CEST
# mode: julia
	ct1
# time: 2023-05-31 04:37:58 CEST
# mode: julia
	BB = get_array(ct1)
# time: 2023-05-31 04:40:52 CEST
# mode: julia
	AAAA = ThirdOrderTensorValue{2,2,2,Int}(1,2.0,3,4,5,6,7,8)
# time: 2023-05-31 04:41:14 CEST
# mode: julia
	using Gridap.TensorValues
# time: 2023-05-31 04:41:15 CEST
# mode: julia
	AAAA = ThirdOrderTensorValue{2,2,2,Int}(1,2.0,3,4,5,6,7,8)
# time: 2023-05-31 04:41:22 CEST
# mode: julia
	get_array(AAAA)
# time: 2023-05-31 04:41:38 CEST
# mode: julia
	get_matrix(AAAA)
# time: 2023-05-31 04:41:43 CEST
# mode: julia
	AAAA
# time: 2023-05-31 04:41:49 CEST
# mode: julia
	AAAA[1,1,1]
# time: 2023-05-31 04:42:22 CEST
# mode: julia
	arrct2
# time: 2023-05-31 04:42:23 CEST
# mode: julia
	ct2
# time: 2023-05-31 04:42:30 CEST
# mode: julia
	get_array(ct2)ct2
# time: 2023-05-31 04:42:34 CEST
# mode: julia
	get_array(ct2)
# time: 2023-05-31 04:42:56 CEST
# mode: julia
	get_array(ct2)[1,1]
# time: 2023-05-31 04:43:00 CEST
# mode: julia
	get_array(ct2)[1,1,1]
# time: 2023-05-31 04:43:03 CEST
# mode: julia
	get_array(ct2)[1,1,1,1]
# time: 2023-05-31 04:43:05 CEST
# mode: julia
	get_array(ct2)[1,1,1,1,1]
# time: 2023-05-31 04:43:13 CEST
# mode: julia
	get_array(ct2)
# time: 2023-05-31 04:44:29 CEST
# mode: julia
	AAAA = TensorValue{2,2}(1,2,3,4)
# time: 2023-05-31 04:44:40 CEST
# mode: julia
	get_array(AAAA)
# time: 2023-05-31 04:44:54 CEST
# mode: julia
	ct1[1,2,3]
# time: 2023-05-31 04:45:04 CEST
# mode: julia
	ct1[1,2,]
# time: 2023-05-31 04:45:05 CEST
# mode: julia
	ct1[1,2]
# time: 2023-05-31 04:45:13 CEST
# mode: julia
	ct1
# time: 2023-05-31 04:45:18 CEST
# mode: julia
	ct1[1,1,1,]
# time: 2023-05-31 04:45:20 CEST
# mode: julia
	ct1[1,1,1,1]
# time: 2023-05-31 04:45:23 CEST
# mode: julia
	ct1[1]
# time: 2023-05-31 04:45:27 CEST
# mode: julia
	ct1[1,2]
# time: 2023-05-31 04:45:30 CEST
# mode: julia
	ct1[1,1]
# time: 2023-05-31 04:45:36 CEST
# mode: julia
	ct1[1,1,1]
# time: 2023-05-31 04:45:40 CEST
# mode: julia
	ct1[1,1,1,1,1,1]
# time: 2023-05-31 04:46:01 CEST
# mode: julia
	ct1⋅VectorValue(0,1)
# time: 2023-05-31 04:46:07 CEST
# mode: julia
	ct1⋅VectorValue(0,1,0)
# time: 2023-05-31 04:46:10 CEST
# mode: julia
	ct1⋅VectorValue(0,1,0,0)
# time: 2023-05-31 04:46:13 CEST
# mode: julia
	ct1⋅VectorValue(0,1,0,0,0)
# time: 2023-05-31 04:46:23 CEST
# mode: julia
	VectorValue(0,1)
# time: 2023-05-31 04:46:29 CEST
# mode: julia
	ct2⋅VectorValue(0,1,0,0,0)
# time: 2023-05-31 04:47:01 CEST
# mode: julia
	ct2⋅TensorValue{2,2}(0,1,0,0)
# time: 2023-05-31 05:01:10 CEST
# mode: julia
	ct1
# time: 2023-05-31 05:01:16 CEST
# mode: julia
	to_array(ct1)
# time: 2023-05-31 05:01:21 CEST
# mode: julia
	get_array(ct1)
# time: 2023-05-31 05:01:41 CEST
# mode: julia
	ct1[1,1]
# time: 2023-05-31 05:01:44 CEST
# mode: julia
	ct1[1,2]
# time: 2023-05-31 05:01:48 CEST
# mode: julia
	ct1[1,1]
# time: 2023-05-31 05:01:50 CEST
# mode: julia
	ct1[2,1]
# time: 2023-05-31 05:01:55 CEST
# mode: julia
	ct1[2,2]
# time: 2023-05-31 05:01:59 CEST
# mode: julia
	get_array(ct1)[1,1]
# time: 2023-05-31 05:02:27 CEST
# mode: julia
	TensorValue{2,2}(1,2,3,4)
# time: 2023-05-31 05:02:31 CEST
# mode: julia
	AAAA = TensorValue{2,2}(1,2,3,4)
# time: 2023-05-31 05:02:55 CEST
# mode: julia
	get_array(AAAA)
# time: 2023-05-31 05:06:52 CEST
# mode: julia
	get_array(ct1)
# time: 2023-05-31 05:06:55 CEST
# mode: julia
	ct1
# time: 2023-05-31 05:12:40 CEST
# mode: julia
	[0]
# time: 2023-05-31 05:12:47 CEST
# mode: julia
	ct1[1]
# time: 2023-05-31 05:12:55 CEST
# mode: julia
	get_array(ct1[1])
# time: 2023-05-31 05:13:06 CEST
# mode: julia
	get_array(ct2[1])
# time: 2023-05-31 07:46:36 CEST
# mode: julia
	TensorValue{2,2}(1,-5,6,3)
# time: 2023-05-31 07:46:55 CEST
# mode: julia
	VectorValue(1,4)
# time: 2023-05-31 07:47:01 CEST
# mode: julia
	vvv = TensorValue{2,2}(1,-5,6,3)
# time: 2023-05-31 07:47:06 CEST
# mode: julia
	vv = VectorValue(1,4)
# time: 2023-05-31 07:47:13 CEST
# mode: julia
	vvv⋅vv
# time: 2023-06-01 01:38:04 CEST
# mode: julia
	TensorValue{2,2}(1,2,3,4)
# time: 2023-06-01 01:38:39 CEST
# mode: julia
	A = TensorValue{2,2}(1,2,3,4)
# time: 2023-06-01 01:39:04 CEST
# mode: julia
	B = TensorValue{3,3}(0,0,0,0,0,0,0,0,0)
# time: 2023-06-01 01:39:18 CEST
# mode: julia
	B[1:2,1:2] = A
# time: 2023-06-01 01:39:40 CEST
# mode: julia
	B[1,1] = 1
# time: 2023-06-01 01:42:32 CEST
# mode: julia
	Barr = zeros(3,3)
# time: 2023-06-01 01:42:42 CEST
# mode: julia
	Barr[1:2,1:2] = A
# time: 2023-06-01 01:42:47 CEST
# mode: julia
	Barr[1:2,1:2] = getarray(A)
# time: 2023-06-01 01:42:51 CEST
# mode: julia
	Barr[1:2,1:2] = get_array(A)
# time: 2023-06-01 01:42:56 CEST
# mode: julia
	Barr
# time: 2023-06-01 02:13:38 CEST
# mode: julia
	Array{2,2}(1.0, 2.0, 3.0, 1.0)
# time: 2023-06-01 02:14:02 CEST
# mode: julia
	Array{Float64 , (2,2)}(1.0, 2.0, 3.0, 1.0)
# time: 2023-06-01 02:14:15 CEST
# mode: julia
	Array{Float64 , 2)}(1.0, 2.0, 3.0, 1.0)
# time: 2023-06-01 02:14:19 CEST
# mode: julia
	Array{Float64 , 2}(1.0, 2.0, 3.0, 1.0)
# time: 2023-06-01 02:14:24 CEST
# mode: julia
	Array{Float64 , 2}(1.0 1.0, 2.0, 3.0, 1.0)
# time: 2023-06-01 02:14:28 CEST
# mode: julia
	Array{Float64 , 2}(1.0, 1.0, 2.0, 3.0, 1.0)
# time: 2023-06-01 02:14:37 CEST
# mode: julia
	A
# time: 2023-06-01 02:14:43 CEST
# mode: julia
	get_array(A)
# time: 2023-06-01 02:17:32 CEST
# mode: julia
	Aarr = [1.0, 2.0; 3.0, 4.0]
# time: 2023-06-01 02:17:41 CEST
# mode: julia
	Aarr = [1.0 2.0, 3.0 4.0]
# time: 2023-06-01 02:17:56 CEST
# mode: julia
	Aarr = [[1.0, 2.0], [3.0, 4.0]]
# time: 2023-06-01 02:59:29 CEST
# mode: julia
	cf1
# time: 2023-06-01 02:59:36 CEST
# mode: julia
	ct1
# time: 2023-06-01 02:59:40 CEST
# mode: julia
	ct1[1]
# time: 2023-06-01 03:00:35 CEST
# mode: julia
	invtesD
# time: 2023-06-01 03:01:59 CEST
# mode: julia
	CTf_2D
# time: 2023-06-01 03:02:01 CEST
# mode: julia
	CTf_2D[1]
# time: 2023-06-01 03:02:26 CEST
# mode: julia
	CTf_2D[1](Point(1))
# time: 2023-06-01 03:02:32 CEST
# mode: julia
	CTf_2D[1](Point(1,1,1))
# time: 2023-06-01 03:02:48 CEST
# mode: julia
	ct1_2D
# time: 2023-06-01 03:02:52 CEST
# mode: julia
	ct1_2D[1]
# time: 2023-06-01 03:08:31 CEST
# mode: julia
	Matrix{2,2}
# time: 2023-06-01 03:08:41 CEST
# mode: julia
	Matrix{2,2}(1,2,3,4)
# time: 2023-06-01 03:08:49 CEST
# mode: julia
	Matrix{2}(1,2,3,4)
# time: 2023-06-01 03:15:43 CEST
# mode: julia
	AAA = zeros(2,2)
# time: 2023-06-01 03:15:58 CEST
# mode: julia
	AAA[1,1]=5
# time: 2023-06-01 03:16:03 CEST
# mode: julia
	AAA[2,2]=7
# time: 2023-06-01 03:16:15 CEST
# mode: julia
	bbb = VectorValue(1,2)
# time: 2023-06-01 03:16:22 CEST
# mode: julia
	AAA⋅bbb
# time: 2023-06-01 03:16:35 CEST
# mode: julia
	AAA*bbb
# time: 2023-06-01 03:16:50 CEST
# mode: julia
	typeof(AAA)
# time: 2023-06-01 03:17:05 CEST
# mode: julia
	typeof(bbb)
# time: 2023-06-01 03:17:17 CEST
# mode: julia
	get_array(AAA)
# time: 2023-06-01 03:57:27 CEST
# mode: julia
	CellField(3.0,Ω)
# time: 2023-06-01 04:23:26 CEST
# mode: julia
	_my_tensor∘(intrf₀.zf, intrf₀.CTf_2D, intrf₀.invD)
# time: 2023-06-01 04:23:44 CEST
# mode: julia
	Interface._my_tensor∘(intrf₀.zf, intrf₀.CTf_2D, intrf₀.invD)
# time: 2023-06-01 04:23:52 CEST
# mode: julia
	modInterface._my_tensor∘(intrf₀.zf, intrf₀.CTf_2D, intrf₀.invD)
# time: 2023-06-01 04:24:14 CEST
# mode: julia
	pts = get_cell_points(Ω)
# time: 2023-06-01 04:24:25 CEST
# mode: julia
	aux = modInterface._my_tensor∘(intrf₀.zf, intrf₀.CTf_2D, intrf₀.invD)(pts)
# time: 2023-06-01 04:25:20 CEST
# mode: julia
	invD = CellField(intrf₀.invD,intrf₀.Ω)
# time: 2023-06-01 04:25:27 CEST
# mode: julia
	aux = modInterface._my_tensor∘(intrf₀.zf, intrf₀.CTf_2D, invD)(pts)
# time: 2023-06-01 04:25:37 CEST
# mode: julia
	aux = modInterface._my_tensor∘(intrf₀.zf, intrf₀.CTf_2D, invD)
# time: 2023-06-01 04:25:44 CEST
# mode: julia
	aux(pts)
# time: 2023-06-01 04:25:48 CEST
# mode: julia
	aux(pts)[1]
# time: 2023-06-01 04:25:51 CEST
# mode: julia
	aux(pts)[1][1]
# time: 2023-06-01 04:25:59 CEST
# mode: julia
	get_array(aux(pts)[1][1])
# time: 2023-06-01 04:28:58 CEST
# mode: julia
	aux = contribute_matrix(intrf₀, (u,λ₀,λ₁,λ₂,λ₃), (v,μ₀,μ₁,μ₂,μ₃), 1, 2)
# time: 2023-06-01 04:29:29 CEST
# mode: julia
	x = get_trial_fe_basis(X)
# time: 2023-06-01 04:29:36 CEST
# mode: julia
	x = get_trial_fe_basis(U)
# time: 2023-06-01 04:29:45 CEST
# mode: julia
	y = get_fe_basis(V)
# time: 2023-06-01 04:29:56 CEST
# mode: julia
	aux = contribute_matrix(intrf₀, x, y, 1, 2)
# time: 2023-06-01 04:30:14 CEST
# mode: julia
	pts = get_cell_points(Γf)
# time: 2023-06-01 04:30:33 CEST
# mode: julia
	pts = get_cell_points(intrf₀.Γf)
# time: 2023-06-01 04:30:37 CEST
# mode: julia
	aux(pts)
# time: 2023-06-01 04:31:01 CEST
# mode: julia
	get_array(aux)
# time: 2023-06-01 04:32:16 CEST
# mode: julia
	x
# time: 2023-06-01 04:32:19 CEST
# mode: julia
	x[1]
# time: 2023-06-01 04:32:26 CEST
# mode: julia
	x[1](pts)
# time: 2023-06-01 04:32:29 CEST
# mode: julia
	x[1](pts)[1]
# time: 2023-06-01 04:32:59 CEST
# mode: julia
	x[1](pts)[1][1]
# time: 2023-06-01 04:33:04 CEST
# mode: julia
	x[1](pts)[1]
# time: 2023-06-01 04:33:16 CEST
# mode: julia
	x[1](pts)[1][1]
# time: 2023-06-01 04:33:18 CEST
# mode: julia
	x[1](pts)[1][1][1]
# time: 2023-06-02 02:45:42 CEST
# mode: julia
	aaa = zeros(2,2)
# time: 2023-06-02 02:45:58 CEST
# mode: julia
	aaa = [[1.0,2.0],[3.0,5.0]]
# time: 2023-06-02 02:46:07 CEST
# mode: julia
	aaa = zeros(2,2)
# time: 2023-06-02 02:46:31 CEST
# mode: julia
	aaa[:,:] = [1 2;3 4]
# time: 2023-06-02 02:48:51 CEST
# mode: julia
	arrayA
# time: 2023-06-02 06:28:55 CEST
# mode: julia
	CTf_2D
# time: 2023-06-02 06:28:59 CEST
# mode: julia
	CTf_2D[1]
# time: 2023-06-02 06:29:28 CEST
# mode: julia
	pts
# time: 2023-06-02 06:29:32 CEST
# mode: julia
	pts[1]
# time: 2023-06-02 06:29:34 CEST
# mode: julia
	pts
# time: 2023-06-02 06:36:18 CEST
# mode: julia
	pts[1]
# time: 2023-06-02 06:36:21 CEST
# mode: julia
	pts
# time: 2023-06-02 06:39:28 CEST
# mode: julia
	CTf_2D
# time: 2023-06-02 06:39:31 CEST
# mode: julia
	CTf_2D[1]
# time: 2023-06-02 06:47:41 CEST
# mode: julia
	_proj(CTf_2D[1])
# time: 2023-06-02 06:48:19 CEST
# mode: julia
	Ψ
# time: 2023-06-02 06:48:24 CEST
# mode: julia
	Ψ.model
# time: 2023-06-02 06:48:30 CEST
# mode: julia
	Ψ.tface_to_mface
# time: 2023-06-02 06:51:16 CEST
# mode: julia
	Ψ
# time: 2023-06-02 06:51:50 CEST
# mode: julia
	Ψ.tface_to_mface
# time: 2023-06-02 06:53:08 CEST
# mode: julia
	_proj(CTf_2D[1])
# time: 2023-06-02 06:55:01 CEST
# mode: julia
	CTf_2D[1]
# time: 2023-06-02 06:59:21 CEST
# mode: julia
	_proj(CTf_2D[1])
# time: 2023-06-02 07:02:08 CEST
# mode: julia
	Ψ
# time: 2023-06-02 07:02:11 CEST
# mode: julia
	Ψ.model
# time: 2023-06-02 07:02:16 CEST
# mode: julia
	model == Ψ.model
# time: 2023-06-02 07:02:51 CEST
# mode: julia
	step_field(zf_l,z_val,Ω)*CTf_2D[1]
# time: 2023-06-02 07:03:02 CEST
# mode: julia
	step_field(zf_l,0.0,Ω)*CTf_2D[1]
# time: 2023-06-02 07:03:13 CEST
# mode: julia
	aux = step_field(zf_l,0.0,Ω)*CTf_2D[1]
# time: 2023-06-02 07:04:14 CEST
# mode: julia
	pts = get_cell_points(Ψ)
# time: 2023-06-02 07:04:17 CEST
# mode: julia
	aux(pts)
# time: 2023-06-02 07:04:48 CEST
# mode: julia
	change_domain(pts,Ω,PhysicalDomain())
# time: 2023-06-02 07:07:43 CEST
# mode: julia
	pts
# time: 2023-06-02 07:07:49 CEST
# mode: julia
	pts.cell_phys_point
# time: 2023-06-02 07:08:26 CEST
# mode: julia
	tglue = get_glue(Ψ)
# time: 2023-06-02 07:08:33 CEST
# mode: julia
	tglue = Geometry.get_glue(Ψ)
# time: 2023-06-02 07:09:33 CEST
# mode: julia
	tglue = get_glue(Ψ,Val(2))
# time: 2023-06-02 07:09:46 CEST
# mode: julia
	tglue
# time: 2023-06-02 07:09:51 CEST
# mode: julia
	tglue = get_glue(Ψ,Val(1))
# time: 2023-06-02 07:12:17 CEST
# mode: julia
	tglue.mface_to_tface
# time: 2023-06-02 07:12:26 CEST
# mode: julia
	tglue.tface_to_mface
# time: 2023-06-02 07:13:14 CEST
# mode: julia
	tglue.tface_to_mface_map
# time: 2023-06-02 07:22:50 CEST
# mode: julia
	ψ
# time: 2023-06-02 07:22:57 CEST
# mode: julia
	Ψ
# time: 2023-06-02 07:23:00 CEST
# mode: julia
	Ψ.glue
# time: 2023-06-02 07:23:24 CEST
# mode: julia
	aux = step_field(zf_l,0.0,Ω)*CTf_2D[1]
# time: 2023-06-02 07:23:33 CEST
# mode: julia
	pts = get_cell_points(Ψ)
# time: 2023-06-02 07:23:40 CEST
# mode: julia
	aux(pts)
# time: 2023-06-02 07:23:54 CEST
# mode: julia
	aux(pts)[1]
# time: 2023-06-02 07:24:25 CEST
# mode: julia
	pts
# time: 2023-06-02 07:24:30 CEST
# mode: julia
	pts.cell_phys_point
# time: 2023-06-02 07:24:51 CEST
# mode: julia
	aux(pts)
# time: 2023-06-05 05:05:01 CEST
# mode: julia
	CTf_2D[1]
# time: 2023-06-05 05:05:08 CEST
# mode: julia
	typeof(CTf_2D[1])
# time: 2023-06-05 05:16:24 CEST
# mode: julia
	intrf₀.Aa
# time: 2023-06-05 06:01:49 CEST
# mode: julia
	arrayD
# time: 2023-06-05 06:01:57 CEST
# mode: julia
	intrf₀.invD
# time: 2023-06-05 06:02:10 CEST
# mode: julia
	invtesD
# time: 2023-06-05 06:03:33 CEST
# mode: julia
	intrf₀.Da
# time: 2023-06-05 06:04:02 CEST
# mode: julia
	arrayDa
# time: 2023-06-05 06:04:12 CEST
# mode: julia
	arrayDb
# time: 2023-06-05 06:04:14 CEST
# mode: julia
	arrayDd
# time: 2023-06-05 06:04:17 CEST
# mode: julia
	intrf₀.Dd
# time: 2023-06-05 06:04:46 CEST
# mode: julia
	intrf₀.Aa
# time: 2023-06-05 06:05:00 CEST
# mode: julia
	intrf₀.Ab
# time: 2023-06-05 06:05:58 CEST
# mode: julia
	intrf₀.I
# time: 2023-06-05 06:06:01 CEST
# mode: julia
	intrf₀.L
# time: 2023-06-05 06:06:14 CEST
# mode: julia
	intrf₀.h
# time: 2023-06-05 06:10:44 CEST
# mode: julia
	intrf₀.l
# time: 2023-06-05 06:12:00 CEST
# mode: julia
	intrf₀.D1
# time: 2023-06-05 06:12:02 CEST
# mode: julia
	intrf₀.Da
# time: 2023-06-05 06:12:54 CEST
# mode: julia
	intrf₀.Db
# time: 2023-06-05 06:12:56 CEST
# mode: julia
	intrf₀.Dd
# time: 2023-06-05 06:13:12 CEST
# mode: julia
	intrf₀.invD
# time: 2023-06-05 06:13:42 CEST
# mode: julia
	inv(get_array(intrf₀.invD))
# time: 2023-06-05 06:13:58 CEST
# mode: julia
	inv(get_array(intrf₀.invA))
# time: 2023-06-05 06:14:49 CEST
# mode: julia
	get_array(intrf₀.Aa)
# time: 2023-06-05 06:14:58 CEST
# mode: julia
	get_array(intrf₀.Ab)
# time: 2023-06-05 06:15:48 CEST
# mode: julia
	inv(get_array(intrf₀.invA))
# time: 2023-06-05 07:04:41 CEST
# mode: julia
	(i == 1)
# time: 2023-06-05 07:04:44 CEST
# mode: julia
	(2 == 1)
# time: 2023-06-05 07:04:59 CEST
# mode: julia
	(2 == 1) | (1 == 1)
# time: 2023-06-05 08:11:41 CEST
# mode: julia
	iext_iint_to_dof
# time: 2023-06-05 08:11:49 CEST
# mode: julia
	iext_iint_to_dof[1]
# time: 2023-06-05 08:12:01 CEST
# mode: julia
	nodes_loc
# time: 2023-06-06 01:53:08 CEST
# mode: julia
	A = [0]*4
# time: 2023-06-06 01:53:13 CEST
# mode: julia
	A
# time: 2023-06-06 01:54:47 CEST
# mode: julia
	5%3
# time: 2023-06-06 01:54:50 CEST
# mode: julia
	5%5
# time: 2023-06-06 01:54:55 CEST
# mode: julia
	5%1
# time: 2023-06-06 01:54:57 CEST
# mode: julia
	5%5
# time: 2023-06-06 01:55:00 CEST
# mode: julia
	5%6
# time: 2023-06-06 02:00:27 CEST
# mode: julia
	AAA = [1,2,3,4,5]
# time: 2023-06-06 02:00:31 CEST
# mode: julia
	AAA[end]
# time: 2023-06-06 02:00:34 CEST
# mode: julia
	AAA[end+1]
# time: 2023-06-06 02:00:43 CEST
# mode: julia
	AAA[end-1]
# time: 2023-06-06 02:12:29 CEST
# mode: julia
	5%4
# time: 2023-06-06 02:12:47 CEST
# mode: julia
	6%4
# time: 2023-06-06 02:13:19 CEST
# mode: julia
	(4+1)%4
# time: 2023-06-06 02:13:26 CEST
# mode: julia
	(4+4)%4
# time: 2023-06-06 02:18:45 CEST
# mode: julia
	(1+4)%4
# time: 2023-06-06 02:19:33 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 02:19:42 CEST
# mode: julia
	iext_dofs
# time: 2023-06-06 02:21:59 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 02:22:02 CEST
# mode: julia
	iext_dofs
# time: 2023-06-06 02:22:17 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 02:22:32 CEST
# mode: julia
	iext_dofs
# time: 2023-06-06 02:22:34 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 02:23:35 CEST
# mode: julia
	iext_dofs
# time: 2023-06-06 02:23:37 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 02:31:35 CEST
# mode: julia
	(1-1)%4
# time: 2023-06-06 02:31:44 CEST
# mode: julia
	(1-2)%4
# time: 2023-06-06 02:32:47 CEST
# mode: julia
	(4+(1-1))%4
# time: 2023-06-06 02:32:51 CEST
# mode: julia
	(4+(1-1))%4 +1
# time: 2023-06-06 02:32:56 CEST
# mode: julia
	(4+(2-1))%4 +1
# time: 2023-06-06 02:33:15 CEST
# mode: julia
	(4+(1-1))
# time: 2023-06-06 02:33:18 CEST
# mode: julia
	(4+(2-1))
# time: 2023-06-06 02:34:37 CEST
# mode: julia
	(1+4)%4+1
# time: 2023-06-06 02:34:54 CEST
# mode: julia
	(1+5)%4+1
# time: 2023-06-06 02:35:02 CEST
# mode: julia
	(2+5)%4+1
# time: 2023-06-06 02:35:07 CEST
# mode: julia
	(3+5)%4+1
# time: 2023-06-06 02:35:20 CEST
# mode: julia
	(1+5)%4+1
# time: 2023-06-06 02:35:23 CEST
# mode: julia
	(1+3)%4+1
# time: 2023-06-06 02:35:30 CEST
# mode: julia
	(1+5)%4+1
# time: 2023-06-06 02:35:33 CEST
# mode: julia
	(1+6)%4+1
# time: 2023-06-06 02:35:36 CEST
# mode: julia
	(2+6)%4+1
# time: 2023-06-06 02:35:39 CEST
# mode: julia
	(3+6)%4+1
# time: 2023-06-06 02:35:41 CEST
# mode: julia
	(4+6)%4+1
# time: 2023-06-06 02:35:45 CEST
# mode: julia
	(1+6)%4+1
# time: 2023-06-06 02:36:56 CEST
# mode: julia
	interfaces[end-(1-1)]
# time: 2023-06-06 02:37:03 CEST
# mode: julia
	interfaces[end-(2-1)]
# time: 2023-06-06 02:37:20 CEST
# mode: julia
	AAA = [1,2,3,4]
# time: 2023-06-06 02:37:25 CEST
# mode: julia
	AAA[end-(2-1)]
# time: 2023-06-06 02:37:28 CEST
# mode: julia
	AAA[end-(1-1)]
# time: 2023-06-06 02:37:55 CEST
# mode: julia
	interfaces[end-(1-1)]
# time: 2023-06-06 02:38:01 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 03:00:49 CEST
# mode: julia
	iext_dofs
# time: 2023-06-06 03:01:05 CEST
# mode: julia
	iext_dofs[1]
# time: 2023-06-06 03:02:14 CEST
# mode: julia
	vectorsB
# time: 2023-06-06 03:05:01 CEST
# mode: julia
	[1,2,3]
# time: 2023-06-06 03:05:03 CEST
# mode: julia
	[1,2,3]-1
# time: 2023-06-06 03:05:07 CEST
# mode: julia
	[1,2,3].-1
# time: 2023-06-06 03:05:36 CEST
# mode: julia
	vectorsB
# time: 2023-06-06 03:41:57 CEST
# mode: julia
	mapping = 1:ext_discrt[i]
# time: 2023-06-06 03:42:06 CEST
# mode: julia
	mapping = 1:ext_discrt[1] + 1
# time: 2023-06-06 03:42:20 CEST
# mode: julia
	mapping[1]
# time: 2023-06-06 03:42:23 CEST
# mode: julia
	mapping[2]
# time: 2023-06-06 03:42:25 CEST
# mode: julia
	mapping[3]
# time: 2023-06-06 03:52:06 CEST
# mode: julia
	intrf_map
# time: 2023-06-06 03:52:22 CEST
# mode: julia
	mapping1
# time: 2023-06-06 03:57:58 CEST
# mode: julia
	intrf_map
# time: 2023-06-06 04:00:28 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 04:00:49 CEST
# mode: julia
	intrf_map
# time: 2023-06-06 04:02:04 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 04:06:30 CEST
# mode: julia
	intrf_map
# time: 2023-06-06 04:24:52 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 04:46:50 CEST
# mode: julia
	vectorsB
# time: 2023-06-06 04:48:21 CEST
# mode: julia
	intrf_map[1]
# time: 2023-06-06 04:48:25 CEST
# mode: julia
	intrf_map[1,2]
# time: 2023-06-06 04:48:29 CEST
# mode: julia
	intrf_map[1][2]
# time: 2023-06-06 04:49:53 CEST
# mode: julia
	intrf_map[1,2]
# time: 2023-06-06 04:49:56 CEST
# mode: julia
	iext_intrf
# time: 2023-06-06 04:50:17 CEST
# mode: julia
	vectorsB
# time: 2023-06-06 05:13:50 CEST
# mode: julia
	vectorsB[1,1]
# time: 2023-06-06 05:14:27 CEST
# mode: julia
	vectorsB[1][1]
# time: 2023-06-06 05:14:41 CEST
# mode: julia
	vectorsB[1][4]
# time: 2023-06-06 05:25:36 CEST
# mode: julia
	λh
# time: 2023-06-06 05:25:40 CEST
# mode: julia
	λh[1]
# time: 2023-06-06 05:25:58 CEST
# mode: julia
	uh
# time: 2023-06-06 05:26:03 CEST
# mode: julia
	λh
# time: 2023-06-06 06:25:12 CEST
# mode: julia
	λ1
# time: 2023-06-06 06:25:22 CEST
# mode: julia
	λ1.fe_space
# time: 2023-06-06 06:25:28 CEST
# mode: julia
	λ1.minus
# time: 2023-06-06 06:25:34 CEST
# mode: julia
	λ1[1]
# time: 2023-06-06 06:25:39 CEST
# mode: julia
	λ1
# time: 2023-06-06 06:25:48 CEST
# mode: julia
	λ1[1]
# time: 2023-06-06 06:27:07 CEST
# mode: julia
	λ₀
# time: 2023-06-07 09:41:31 CEST
# mode: julia
	contraA
# time: 2023-06-07 09:41:37 CEST
# mode: julia
	contrA
# time: 2023-06-07 09:41:48 CEST
# mode: julia
	elementA = first(contrA.dict).second[1]
# time: 2023-06-07 09:41:57 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1]]
# time: 2023-06-07 09:42:00 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1]
# time: 2023-06-07 09:42:32 CEST
# mode: julia
	elementA = first(contrA.dict).second[1]
# time: 2023-06-07 09:42:44 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][3,3]
# time: 2023-06-07 09:43:00 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1]
# time: 2023-06-07 09:47:50 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1][1,9]
# time: 2023-06-07 09:47:58 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1][1,8]
# time: 2023-06-07 09:48:03 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1][1,5]
# time: 2023-06-07 09:48:14 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1][1,1]
# time: 2023-06-07 09:48:18 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1][1,2]
# time: 2023-06-07 09:48:23 CEST
# mode: julia
	elementA = first(contrA.dict).second[1][1,1][1,:]
# time: 2023-06-13 04:46:14 CEST
# mode: julia
	uh
# time: 2023-06-13 04:46:25 CEST
# mode: julia
	uh.cell_dof_values
# time: 2023-06-13 04:46:41 CEST
# mode: julia
	λ₀.cell_dof_values
# time: 2023-06-13 04:47:16 CEST
# mode: julia
	λ₀.free_values
# time: 2023-06-13 05:14:33 CEST
# mode: julia
	vectorsB
# time: 2023-06-13 05:17:09 CEST
# mode: julia
	λ₀.free_values
# time: 2023-06-13 05:18:48 CEST
# mode: julia
	vectorsB
# time: 2023-06-13 05:22:25 CEST
# mode: julia
	vectorsB[1:30,:]
# time: 2023-06-13 05:22:34 CEST
# mode: julia
	vectorsB[1:5,:]
# time: 2023-06-13 05:22:40 CEST
# mode: julia
	vectorsB[6:10,:]
# time: 2023-06-13 05:22:48 CEST
# mode: julia
	vectorsB[11:15,:]
# time: 2023-06-13 05:23:22 CEST
# mode: julia
	vectorsB[16:20,:]
# time: 2023-06-13 05:25:13 CEST
# mode: julia
	vectorsB[1:5,:]
# time: 2023-06-13 05:27:27 CEST
# mode: julia
	vectorsB[6:10,:]
# time: 2023-06-13 05:27:33 CEST
# mode: julia
	vectorsB[11:15,:]
# time: 2023-06-13 05:32:39 CEST
# mode: julia
	vectorsB[1:5,:]
# time: 2023-06-13 06:07:07 CEST
# mode: julia
	vectorsB[1:30,:]
# time: 2023-06-15 03:32:45 CEST
# mode: julia
	intrf₀
# time: 2023-06-15 03:32:59 CEST
# mode: julia
	intrf₀.Γc
# time: 2023-06-15 03:35:48 CEST
# mode: julia
	Γc
# time: 2023-06-15 03:35:54 CEST
# mode: julia
	num_cells(Γc)
# time: 2023-06-15 03:44:10 CEST
# mode: julia
	cface_model
# time: 2023-06-15 03:44:19 CEST
# mode: julia
	cface_model.grid_topology
# time: 2023-06-15 03:44:39 CEST
# mode: julia
	get_faces(cface_model.grid_topology,2,2)
# time: 2023-06-15 03:44:48 CEST
# mode: julia
	Geometry.get_faces(cface_model.grid_topology,2,2)
# time: 2023-06-15 03:45:07 CEST
# mode: julia
	Geometry.get_faces(cface_model.grid_topology,2,1)
# time: 2023-06-15 03:45:54 CEST
# mode: julia
	Geometry.get_faces(cface_model.grid_topology,2,0)
# time: 2023-06-15 04:16:47 CEST
# mode: julia
	to__
# time: 2023-06-16 04:17:31 CEST
# mode: julia
	vectorsB
# time: 2023-06-16 04:35:59 CEST
# mode: julia
	(true)
# time: 2023-06-16 04:36:08 CEST
# mode: julia
	typeof((true))
# time: 2023-06-16 04:44:54 CEST
# mode: julia
	uC
# time: 2023-06-16 04:45:01 CEST
# mode: julia
	uC.free_values
# time: 2023-06-16 05:24:39 CEST
# mode: julia
	glue₀
# time: 2023-06-16 05:24:44 CEST
# mode: julia
	glue₀.is_refined
# time: 2023-06-16 07:04:55 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:05:18 CEST
# mode: julia
	child_ids
# time: 2023-06-16 07:05:24 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:05:32 CEST
# mode: julia
	child_ids
# time: 2023-06-16 07:05:38 CEST
# mode: julia
	rrules
# time: 2023-06-16 07:05:44 CEST
# mode: julia
	rrules.value
# time: 2023-06-16 07:05:48 CEST
# mode: julia
	rrules.axes
# time: 2023-06-16 07:05:54 CEST
# mode: julia
	rrules
# time: 2023-06-16 07:06:18 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:07:24 CEST
# mode: julia
	AAA = zeros(4)
# time: 2023-06-16 07:07:29 CEST
# mode: julia
	AAA += 2
# time: 2023-06-16 07:07:39 CEST
# mode: julia
	AAA .+= 2
# time: 2023-06-16 07:08:13 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:08:26 CEST
# mode: julia
	child_ids
# time: 2023-06-16 07:12:31 CEST
# mode: julia
	c2f_faces
# time: 2023-06-16 07:15:45 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:16:58 CEST
# mode: julia
	child_ids
# time: 2023-06-16 07:17:28 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:17:30 CEST
# mode: julia
	child_ids
# time: 2023-06-16 07:21:48 CEST
# mode: julia
	c2f_faces
# time: 2023-06-16 07:23:23 CEST
# mode: julia
	n2o_faces_g
# time: 2023-06-16 07:23:30 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:24:14 CEST
# mode: julia
	n2o_faces_g
# time: 2023-06-16 07:24:17 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:33:04 CEST
# mode: julia
	n2o_faces_g
# time: 2023-06-16 07:41:46 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:41:51 CEST
# mode: julia
	child_ids
# time: 2023-06-16 07:42:45 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:42:47 CEST
# mode: julia
	child_ids
# time: 2023-06-16 07:43:02 CEST
# mode: julia
	n2o_faces
# time: 2023-06-16 07:43:09 CEST
# mode: julia
	n2o_faces_g
# time: 2023-06-16 07:43:53 CEST
# mode: julia
	child_ids_g
# time: 2023-06-16 07:44:24 CEST
# mode: julia
	child_ids_g[30:40]
# time: 2023-06-16 07:44:32 CEST
# mode: julia
	child_ids_g[1:8]
# time: 2023-06-16 07:44:57 CEST
# mode: julia
	child_ids_g
# time: 2023-06-16 07:45:39 CEST
# mode: julia
	c2f_faces
# time: 2023-06-16 07:46:13 CEST
# mode: julia
	append!(c2f_faces, c2f_faces)
# time: 2023-06-16 07:48:12 CEST
# mode: julia
	c2f_faces
# time: 2023-06-16 07:48:16 CEST
# mode: julia
	child_ids_g
# time: 2023-06-16 07:48:26 CEST
# mode: julia
	n2o_faces_g
# time: 2023-06-16 08:00:17 CEST
# mode: julia
	c2f_faces_list
# time: 2023-06-16 08:01:00 CEST
# mode: julia
	c2f_faces_list...
# time: 2023-06-16 08:03:37 CEST
# mode: julia
	c2f_faces_list
# time: 2023-06-16 08:03:40 CEST
# mode: julia
	c2f_faces_list[1]
# time: 2023-06-16 08:04:27 CEST
# mode: julia
	append_tables_globally(c2f_faces_list[1])
# time: 2023-06-16 08:04:34 CEST
# mode: julia
	append_tables_globally(c2f_faces_list[1], c2f_faces_list[2])
# time: 2023-06-16 08:04:48 CEST
# mode: julia
	append_tables_globally(c2f_faces_list[1], c2f_faces_list[2], c2f_faces_list[3])
# time: 2023-06-16 08:04:56 CEST
# mode: julia
	append_tables_globally(c2f_faces_list[1], c2f_faces_list[3])
# time: 2023-06-16 08:05:03 CEST
# mode: julia
	append_tables_globally(c2f_faces_list[3])
# time: 2023-06-16 08:05:09 CEST
# mode: julia
	c2f_faces_list[3]
# time: 2023-06-16 08:10:27 CEST
# mode: julia
	Table(c2f_faces_list[3])
# time: 2023-06-16 08:13:29 CEST
# mode: julia
	sortperm([2,3,4])
# time: 2023-06-16 08:13:41 CEST
# mode: julia
	sortperm([2,4,1,3])
# time: 2023-06-16 08:13:51 CEST
# mode: julia
	sortperm([2,4,1,3],rev)
# time: 2023-06-16 08:13:56 CEST
# mode: julia
	sortperm([2,4,1,3],:rev)
# time: 2023-06-16 08:14:04 CEST
# mode: julia
	methods(sortperm)
# time: 2023-06-16 08:14:19 CEST
# mode: julia
	sortperm([2,4,1,3],rev)
# time: 2023-06-16 08:14:53 CEST
# mode: julia
	sortperm([2,4,1,3],by=rev)
# time: 2023-06-16 08:14:56 CEST
# mode: julia
	sortperm([2,4,1,3],by=:rev)
# time: 2023-06-16 08:15:09 CEST
# mode: julia
	sortperm([2,4,1,3],rev=true)
# time: 2023-06-16 08:16:08 CEST
# mode: julia
	aux
# time: 2023-06-16 08:17:22 CEST
# mode: julia
	glue
# time: 2023-06-16 08:17:32 CEST
# mode: julia
	glue.n2o_cell_to_child_id
# time: 2023-06-16 08:21:34 CEST
# mode: julia
	glue.n2o_faces_map
# time: 2023-06-16 08:21:54 CEST
# mode: julia
	glue.n2o_faces_map[3]
# time: 2023-06-16 08:22:30 CEST
# mode: julia
	glue.refinement_rules
# time: 2023-06-16 08:29:42 CEST
# mode: julia
	Γf.parent_trian
# time: 2023-06-16 08:41:06 CEST
# mode: julia
	u = get_trial_fe_basis(U)
# time: 2023-06-16 08:41:10 CEST
# mode: julia
	u = get_trial_fe_basis(Uu)
# time: 2023-06-16 08:41:55 CEST
# mode: julia
	v = get_fe_basis(Uu)
# time: 2023-06-16 08:42:18 CEST
# mode: julia
	dΩ = Measure(Ω,3)
# time: 2023-06-16 08:42:34 CEST
# mode: julia
	c = ∫(u⋅v)dΩ
# time: 2023-06-16 08:42:58 CEST
# mode: julia
	c += ∫(u⋅v)dΩ
# time: 2023-06-16 09:08:32 CEST
# mode: julia
	totalColumns
# time: 2023-06-16 09:23:19 CEST
# mode: julia
	intrf₀.rot_arr
# time: 2023-06-16 09:38:05 CEST
# mode: julia
	totalColumns
# time: 2023-06-16 09:56:49 CEST
# mode: julia
	DomainStyle(λ)
# time: 2023-06-16 10:01:14 CEST
# mode: julia
	λf.cell_basis
# time: 2023-06-16 10:01:32 CEST
# mode: julia
	CellData.get_data(λf)
# time: 2023-06-16 10:01:43 CEST
# mode: julia
	typeof(CellData.get_data(λf))
# time: 2023-06-16 10:02:22 CEST
# mode: julia
	typeof(CellData.get_data(λf)) <: AbstractArray{<:Gridap.Field}
# time: 2023-06-16 10:02:27 CEST
# mode: julia
	typeof(CellData.get_data(λf)) <: AbstractArray{<:Gridap.Fields.Field}
# time: 2023-06-16 10:02:41 CEST
# mode: julia
	typeof(CellData.get_data(λf)) <: AbstractArray{<:AbstractArray{<:Gridap.Fields.Field}}
# time: 2023-06-16 10:02:52 CEST
# mode: julia
	CellData.get_data(λf)[1]
# time: 2023-06-16 10:04:01 CEST
# mode: julia
	u
# time: 2023-06-16 10:04:08 CEST
# mode: julia
	u.single_field
# time: 2023-06-16 10:04:14 CEST
# mode: julia
	u.cell_basis
# time: 2023-06-16 10:04:23 CEST
# mode: julia
	u.cell_basis[1]
# time: 2023-06-16 10:18:17 CEST
# mode: julia
	Γf
# time: 2023-06-16 10:18:22 CEST
# mode: julia
	Γf.glue
# time: 2023-06-16 10:18:39 CEST
# mode: julia
	Γf.glue.n2o_faces_map
# time: 2023-06-16 10:18:44 CEST
# mode: julia
	Γf.glue.n2o_faces_map[3]
# time: 2023-06-16 10:19:38 CEST
# mode: julia
	aux
# time: 2023-06-16 10:19:46 CEST
# mode: julia
	aux = λf
# time: 2023-06-16 10:19:51 CEST
# mode: julia
	aux = λf.cell_basis
# time: 2023-06-16 10:19:58 CEST
# mode: julia
	aux.args
# time: 2023-06-16 10:20:10 CEST
# mode: julia
	aux.args[1]
# time: 2023-06-16 10:20:14 CEST
# mode: julia
	aux.args[2]
# time: 2023-06-16 10:20:23 CEST
# mode: julia
	aux.maps
# time: 2023-06-16 10:20:58 CEST
# mode: julia
	aux.args[1]
# time: 2023-06-16 10:21:10 CEST
# mode: julia
	aux.args[1].args
# time: 2023-06-16 10:21:19 CEST
# mode: julia
	aux.args[1].args[1]
# time: 2023-06-16 10:21:50 CEST
# mode: julia
	aux.args[1].args[2]
# time: 2023-06-16 10:22:40 CEST
# mode: julia
	aux.args[1].maps
# time: 2023-06-16 10:23:00 CEST
# mode: julia
	aux.args[1][1]
# time: 2023-06-16 10:23:20 CEST
# mode: julia
	aux.args[1][1] <: AbstractArray{<:Field}
# time: 2023-06-16 10:23:28 CEST
# mode: julia
	aux.args[1][1] <: AbstractArray{<:Gridap.Fields.Field}
# time: 2023-06-16 10:23:45 CEST
# mode: julia
	typeof(aux.args[1][1]) <: AbstractArray{<:Gridap.Fields.Field}
# time: 2023-06-19 02:53:48 CEST
# mode: pkg
	up
# time: 2023-06-19 02:55:47 CEST
# mode: pkg
	status
# time: 2023-06-19 03:19:42 CEST
# mode: julia
	totalColumns
# time: 2023-06-19 03:43:37 CEST
# mode: julia
	aux[1]
# time: 2023-06-19 03:44:23 CEST
# mode: julia
	T1 = aux[1]
# time: 2023-06-19 03:44:37 CEST
# mode: julia
	v1 = VectorValue(1.0,2.0,3.0)
# time: 2023-06-19 03:44:47 CEST
# mode: julia
	v1⋅T1⋅v2
# time: 2023-06-19 03:44:49 CEST
# mode: julia
	v1⋅T1⋅v1
# time: 2023-06-19 03:45:09 CEST
# mode: julia
	T1 = aux[1][1]
# time: 2023-06-19 03:45:12 CEST
# mode: julia
	v1⋅T1⋅v1
# time: 2023-06-19 03:46:20 CEST
# mode: julia
	v1 = VectorValue(1.0,2.0,3.0)
# time: 2023-06-19 03:46:41 CEST
# mode: julia
	v2 = VectorValue(1.0,2.0,3.0,4.0,5.0)
# time: 2023-06-19 03:46:47 CEST
# mode: julia
	v2⋅T1⋅v1
# time: 2023-06-19 03:50:53 CEST
# mode: julia
	num_dofs(reffe)
# time: 2023-06-19 03:50:56 CEST
# mode: julia
	num_dofs(reffe_u)
# time: 2023-06-19 03:51:06 CEST
# mode: julia
	Gridap.ReferenceFEs.num_dofs(reffe_u)
# time: 2023-06-19 03:51:32 CEST
# mode: julia
	Vu.fe_basis
# time: 2023-06-19 03:51:41 CEST
# mode: julia
	Vu.fe_basis.cell_basis.value
# time: 2023-06-19 03:57:50 CEST
# mode: julia
	get_cell_dof_ids(Vλ)
# time: 2023-06-19 03:58:37 CEST
# mode: julia
	get_cell_dof_ids(Vλ,Γf)
# time: 2023-06-19 03:58:54 CEST
# mode: julia
	get_cell_dof_ids(Vλ,Γi)
# time: 2023-06-19 03:58:59 CEST
# mode: julia
	get_cell_dof_ids(Vλ,intrf.Γi)
# time: 2023-06-19 04:08:27 CEST
# mode: julia
	glue.n2o_cell_to_child_id[3]
# time: 2023-06-19 04:08:30 CEST
# mode: julia
	glue.n2o_cell_to_child_id
# time: 2023-06-19 04:08:38 CEST
# mode: julia
	glue.n2o_faces_map[3]
# time: 2023-06-19 04:38:44 CEST
# mode: julia
	b[1216:1220]
# time: 2023-06-19 04:39:40 CEST
# mode: julia
	maximum(b)
# time: 2023-06-19 04:39:48 CEST
# mode: julia
	minimum(b)
# time: 2023-06-19 04:44:12 CEST
# mode: julia
	xC
# time: 2023-06-19 04:44:53 CEST
# mode: julia
	Point(0.0)
# time: 2023-06-19 04:45:03 CEST
# mode: julia
	fun_uC()
# time: 2023-06-19 04:45:14 CEST
# mode: julia
	pointA = Point(0.0)
# time: 2023-06-19 04:45:19 CEST
# mode: julia
	fun_uC(pointA)
# time: 2023-06-19 04:45:36 CEST
# mode: julia
	pointB = Point(0.1)
# time: 2023-06-19 04:45:39 CEST
# mode: julia
	pointB = Point(0.01)
# time: 2023-06-19 04:45:42 CEST
# mode: julia
	fun_uC(pointB)
# time: 2023-06-19 04:47:27 CEST
# mode: julia
	sum(crs_dsc)
# time: 2023-06-19 04:48:12 CEST
# mode: julia
	xC
# time: 2023-06-19 04:48:42 CEST
# mode: julia
	fun_uC(pointB)
# time: 2023-06-19 04:48:48 CEST
# mode: julia
	fun_uC(pointA)
# time: 2023-06-19 04:48:56 CEST
# mode: julia
	fun_uC(Point(0.9))
# time: 2023-06-19 04:50:42 CEST
# mode: julia
	b[1216:1220]
# time: 2023-06-19 04:50:51 CEST
# mode: julia
	b[1211:1220]
# time: 2023-06-19 04:51:09 CEST
# mode: julia
	b[1206:1220]
# time: 2023-06-19 04:51:16 CEST
# mode: julia
	b[1216:1220]
# time: 2023-06-19 04:51:20 CEST
# mode: julia
	b[1216:1225]
# time: 2023-06-19 04:52:55 CEST
# mode: julia
	maximum(b)
# time: 2023-06-19 04:53:00 CEST
# mode: julia
	minimum(b)
# time: 2023-06-19 04:53:55 CEST
# mode: julia
	xC
# time: 2023-06-19 05:00:58 CEST
# mode: julia
	b[1216:end]
# time: 2023-06-19 05:05:51 CEST
# mode: julia
	modelΛc
# time: 2023-06-19 05:05:56 CEST
# mode: julia
	modelΛc.grid
# time: 2023-06-19 05:06:02 CEST
# mode: julia
	modelΛc.grid_topology
# time: 2023-06-19 05:06:08 CEST
# mode: julia
	modelΛc.grid_topology.cell_type
# time: 2023-06-19 05:06:16 CEST
# mode: julia
	modelΛc.grid_topology.vertex_coordinates
# time: 2023-06-19 05:07:51 CEST
# mode: julia
	VΛc
# time: 2023-06-19 05:08:05 CEST
# mode: julia
	VΛc.cell_dofs_ids
# time: 2023-06-19 05:09:03 CEST
# mode: julia
	VΛc.nfree
# time: 2023-06-19 05:10:30 CEST
# mode: julia
	Λf
# time: 2023-06-19 05:10:41 CEST
# mode: julia
	Λf.grid.node_coords
# time: 2023-06-19 05:11:13 CEST
# mode: julia
	totalColumns
# time: 2023-06-19 05:11:26 CEST
# mode: julia
	VΛf
# time: 2023-06-19 05:11:33 CEST
# mode: julia
	VΛf.cell_dofs_ids
# time: 2023-06-19 05:11:47 CEST
# mode: julia
	VΛf.nfree
# time: 2023-06-19 05:15:15 CEST
# mode: julia
	b
# time: 2023-06-19 05:16:49 CEST
# mode: julia
	bVu
# time: 2023-06-19 05:16:52 CEST
# mode: julia
	Vu
# time: 2023-06-19 05:16:59 CEST
# mode: julia
	Vu.nfree
# time: 2023-06-19 05:18:29 CEST
# mode: julia
	num_dofs(U)
# time: 2023-06-19 05:18:34 CEST
# mode: julia
	num_free_dofs(U)
# time: 2023-06-19 05:18:40 CEST
# mode: julia
	num_free_dofs(Uλ)
# time: 2023-06-19 05:20:59 CEST
# mode: julia
	num_free_dofs(Uu)
# time: 2023-06-19 05:23:08 CEST
# mode: julia
	num_free_dofs(Λc)
# time: 2023-06-19 05:23:13 CEST
# mode: julia
	num_free_dofs(VΛc)
# time: 2023-06-19 05:27:19 CEST
# mode: julia
	num_free_dofs(UΛc)
# time: 2023-06-19 05:27:26 CEST
# mode: julia
	modelΛc
# time: 2023-06-19 05:27:33 CEST
# mode: julia
	modelΛc.grid_topology
# time: 2023-06-19 05:27:46 CEST
# mode: julia
	Geometry.get_faces(modelΛc.grid_topology,0,0)
# time: 2023-06-19 05:27:58 CEST
# mode: julia
	Geometry.get_faces(modelΛc.grid_topology,0,1)
# time: 2023-06-19 05:28:41 CEST
# mode: julia
	modelΛc.grid_topology.vertex_coordinates
# time: 2023-06-19 05:30:29 CEST
# mode: julia
	get_cell_dof_ids(VΛc)
# time: 2023-06-19 06:16:00 CEST
# mode: julia
	det(A)
# time: 2023-06-19 06:22:10 CEST
# mode: julia
	contrA
# time: 2023-06-19 06:22:20 CEST
# mode: julia
	contrA.dict
# time: 2023-06-19 06:23:04 CEST
# mode: julia
	first(contrA.dict)
# time: 2023-06-19 06:23:26 CEST
# mode: julia
	first(contrA.dict).second
# time: 2023-06-19 06:23:43 CEST
# mode: julia
	elementA
# time: 2023-06-19 06:23:48 CEST
# mode: julia
	elementA[2,1]
# time: 2023-06-19 06:24:11 CEST
# mode: julia
	elementA[2,1][2,:]
# time: 2023-06-19 06:24:19 CEST
# mode: julia
	mazimum(elementA[2,1][2,:])
# time: 2023-06-19 06:24:23 CEST
# mode: julia
	maximum(elementA[2,1][2,:])
# time: 2023-06-19 06:24:37 CEST
# mode: julia
	maximum(elementA[2,1][1,:])
# time: 2023-06-19 06:24:39 CEST
# mode: julia
	maximum(elementA[2,1][2,:])
# time: 2023-06-19 06:24:41 CEST
# mode: julia
	maximum(elementA[2,1][3,:])
# time: 2023-06-19 06:24:43 CEST
# mode: julia
	maximum(elementA[2,1][4,:])
# time: 2023-06-19 06:24:45 CEST
# mode: julia
	maximum(elementA[2,1][5,:])
# time: 2023-06-19 06:25:19 CEST
# mode: julia
	first(contrA.dict).second[2]
# time: 2023-06-19 06:25:22 CEST
# mode: julia
	first(contrA.dict).second[1]
# time: 2023-06-19 06:25:36 CEST
# mode: julia
	first(contrA.dict).first[1]
# time: 2023-06-19 06:25:42 CEST
# mode: julia
	first(contrA.dict)
# time: 2023-06-19 06:26:25 CEST
# mode: julia
	A
# time: 2023-06-19 06:26:40 CEST
# mode: julia
	A[82:101]
# time: 2023-06-19 06:26:55 CEST
# mode: julia
	A[82:101,:]
# time: 2023-06-19 06:27:05 CEST
# mode: julia
	A[82:101,1]
# time: 2023-06-19 06:27:07 CEST
# mode: julia
	A[82:101,2]
# time: 2023-06-19 06:27:09 CEST
# mode: julia
	A[82:101,3]
# time: 2023-06-19 06:27:11 CEST
# mode: julia
	A[82:101,4]
# time: 2023-06-19 06:27:13 CEST
# mode: julia
	A[82:101,5]
# time: 2023-06-19 06:27:15 CEST
# mode: julia
	A[82:101,6]
# time: 2023-06-19 06:27:30 CEST
# mode: julia
	A[82,:]
# time: 2023-06-19 06:27:40 CEST
# mode: julia
	A[83,:]
# time: 2023-06-19 06:27:43 CEST
# mode: julia
	A[84,:]
# time: 2023-06-19 06:27:54 CEST
# mode: julia
	maximum( A[84,:] )
# time: 2023-06-19 06:27:57 CEST
# mode: julia
	maximum( A[85,:] )
# time: 2023-06-19 06:28:00 CEST
# mode: julia
	maximum( A[87,:] )
# time: 2023-06-19 06:28:04 CEST
# mode: julia
	maximum( A[88,:] )
# time: 2023-06-19 06:28:06 CEST
# mode: julia
	maximum( A[89,:] )
# time: 2023-06-19 06:28:08 CEST
# mode: julia
	maximum( A[80,:] )
# time: 2023-06-19 06:28:17 CEST
# mode: julia
	maximum( A[90,:] )
# time: 2023-06-19 06:28:19 CEST
# mode: julia
	maximum( A[91,:] )
# time: 2023-06-19 06:28:21 CEST
# mode: julia
	maximum( A[92,:] )
# time: 2023-06-19 06:28:23 CEST
# mode: julia
	maximum( A[93,:] )
# time: 2023-06-19 06:28:26 CEST
# mode: julia
	maximum( A[94,:] )
# time: 2023-06-19 06:28:28 CEST
# mode: julia
	maximum( A[95,:] )
# time: 2023-06-19 07:10:29 CEST
# mode: julia
	face_glue.n2o_faces_map
# time: 2023-06-19 07:10:37 CEST
# mode: julia
	face_glue.n2o_faces_map[3]
# time: 2023-06-19 07:14:03 CEST
# mode: julia
	det(A)
# time: 2023-06-19 07:27:32 CEST
# mode: julia
	Γᵢ_glues
# time: 2023-06-19 07:49:57 CEST
# mode: julia
	n2o_faces_local
# time: 2023-06-19 08:13:13 CEST
# mode: julia
	n2o_faces_g
# time: 2023-06-19 08:26:55 CEST
# mode: julia
	global_n2o_faces
# time: 2023-06-19 08:33:12 CEST
# mode: julia
	M = ones(3,3)
# time: 2023-06-19 08:33:35 CEST
# mode: julia
	M_vec = Fill(M,10)
# time: 2023-06-19 08:33:47 CEST
# mode: julia
	M_vec.value
# time: 2023-06-19 08:34:18 CEST
# mode: julia
	M
# time: 2023-06-19 08:34:26 CEST
# mode: julia
	M = ones(3,3)
# time: 2023-06-19 08:34:35 CEST
# mode: julia
	N = 2*ones(4,4)
# time: 2023-06-19 08:35:10 CEST
# mode: julia
	CompressedArray([M,N],[1,2,1,1,1,2,2,1,2])
# time: 2023-06-19 08:41:18 CEST
# mode: julia
	glue.n2o_faces_map[3]
# time: 2023-06-19 08:41:42 CEST
# mode: julia
	partial_glues[1].n2o_faces_map[3]
# time: 2023-06-19 08:41:47 CEST
# mode: julia
	partial_glues[2].n2o_faces_map[3]
# time: 2023-06-19 08:41:55 CEST
# mode: julia
	partial_glues[3].n2o_faces_map[3]
# time: 2023-06-19 08:42:05 CEST
# mode: julia
	partial_glues[4].n2o_faces_map[3]
# time: 2023-06-19 08:42:21 CEST
# mode: julia
	gg = partial_glues[4]
# time: 2023-06-19 08:42:27 CEST
# mode: julia
	gg.o2n_faces_map
# time: 2023-06-19 08:43:37 CEST
# mode: julia
	gg.n2o_cell_to_child_id
# time: 2023-06-19 08:45:34 CEST
# mode: julia
	λf
# time: 2023-06-19 08:45:42 CEST
# mode: julia
	λf.cell_basis
# time: 2023-06-19 08:46:44 CEST
# mode: julia
	Γf
# time: 2023-06-19 08:46:56 CEST
# mode: julia
	Γf.glue.n2o_cell_to_child_id
# time: 2023-06-19 08:47:23 CEST
# mode: julia
	Γf.glue.n2o_faces_map
# time: 2023-06-19 08:47:27 CEST
# mode: julia
	Γf.glue.n2o_faces_map[3]
# time: 2023-06-19 08:48:04 CEST
# mode: julia
	glue.n2o_faces_map[3]
# time: 2023-06-19 08:48:28 CEST
# mode: julia
	Γf.glue.n2o_faces_map[3]
# time: 2023-06-19 08:50:36 CEST
# mode: julia
	b[1216:]
# time: 2023-06-19 08:50:43 CEST
# mode: julia
	b[1216:end]
# time: 2023-06-20 04:29:46 CEST
# mode: julia
	exit()
# time: 2023-06-20 04:29:57 CEST
# mode: julia
	ls
# time: 2023-06-20 04:30:16 CEST
# mode: julia
	include("MDC_SL_Reissner_Inrtf_SuperElement.jl")
# time: 2023-06-20 06:35:37 CEST
# mode: julia
	paraview
# time: 2023-06-20 06:35:44 CEST
# mode: julia
	exit()
# time: 2023-06-20 10:18:49 CEST
# mode: julia
	Uu.nfree
# time: 2023-06-21 03:28:28 CEST
# mode: julia
	intrfA.Da
# time: 2023-06-21 03:37:40 CEST
# mode: julia
	zf
# time: 2023-06-21 04:21:49 CEST
# mode: julia
	intrfA
# time: 2023-06-22 02:30:43 CEST
# mode: pkg
	status
# time: 2023-06-22 02:38:47 CEST
# mode: pkg
	list
# time: 2023-06-22 02:38:52 CEST
# mode: pkg
	?
# time: 2023-06-22 02:39:27 CEST
# mode: pkg
	registry status
# time: 2023-06-22 02:39:37 CEST
# mode: pkg
	pin
# time: 2023-06-22 02:41:19 CEST
# mode: pkg
	status
# time: 2023-06-22 02:41:36 CEST
# mode: pkg
	status --outdated
# time: 2023-06-22 02:46:04 CEST
# mode: julia
	λh
# time: 2023-06-22 02:46:12 CEST
# mode: julia
	λh.free_values
# time: 2023-06-22 02:47:56 CEST
# mode: pkg
	version Gridap
# time: 2023-06-22 02:48:02 CEST
# mode: pkg
	status
# time: 2023-06-22 02:48:16 CEST
# mode: pkg
	help
# time: 2023-06-22 02:50:39 CEST
# mode: pkg
	?\
# time: 2023-06-22 02:50:41 CEST
# mode: pkg
	?
# time: 2023-06-22 02:51:17 CEST
# mode: pkg
	st gridap
# time: 2023-06-22 02:52:22 CEST
# mode: pkg
	exit()
# time: 2023-06-22 02:52:28 CEST
# mode: julia
	exit()
# time: 2023-06-23 05:00:20 CEST
# mode: julia
	a
# time: 2023-06-23 05:01:33 CEST
# mode: julia
	A
# time: 2023-06-23 05:01:40 CEST
# mode: julia
	A[1,1]
# time: 2023-06-23 05:01:47 CEST
# mode: julia
	A[1,27]
# time: 2023-06-23 05:01:51 CEST
# mode: julia
	A[1,28]
# time: 2023-06-23 05:01:55 CEST
# mode: julia
	A[1,27]
# time: 2023-06-23 05:04:26 CEST
# mode: julia
	A[:,:]
# time: 2023-06-26 02:36:57 CEST
# mode: julia
	x = collect(-1.0:0.001:1.0)
# time: 2023-06-26 02:37:23 CEST
# mode: julia
	y = exp(-x.^2)
# time: 2023-06-26 02:37:31 CEST
# mode: julia
	y = exp.(-x.^2)
# time: 2023-06-26 02:37:58 CEST
# mode: julia
	y = exp.(-x.^2/0.1)
# time: 2023-06-26 02:38:03 CEST
# mode: julia
	y = exp.(-x.^2/0.5)
# time: 2023-06-26 02:38:07 CEST
# mode: julia
	y = exp.(-x.^2/0.1)
# time: 2023-06-26 02:57:29 CEST
# mode: julia
	get_cell_dof_ids(Uu,Γ₀_glue)
# time: 2023-06-26 02:57:50 CEST
# mode: julia
	map(length,get_cell_dof_ids(Uu,Γ₀_glue))
# time: 2023-06-26 03:08:49 CEST
# mode: julia
	104-81
# time: 2023-06-26 03:10:14 CEST
# mode: julia
	dof_ids
# time: 2023-06-26 03:16:14 CEST
# mode: julia
	A[1090,1]
# time: 2023-06-26 03:16:17 CEST
# mode: julia
	A[1093,1]
# time: 2023-06-26 03:16:28 CEST
# mode: julia
	A[1093,92]
# time: 2023-06-26 03:16:30 CEST
# mode: julia
	A[1093,82]
# time: 2023-06-26 03:16:41 CEST
# mode: julia
	A[1093,4]
# time: 2023-06-26 04:04:54 CEST
# mode: julia
	reffe_u
# time: 2023-06-26 04:05:19 CEST
# mode: julia
	Uu.fe_dof_basis.cell_dof[1]
# time: 2023-06-26 04:05:46 CEST
# mode: julia
	get_cell_reffe(Uu)
# time: 2023-06-26 04:06:01 CEST
# mode: julia
	reffe_u
# time: 2023-06-26 04:06:44 CEST
# mode: julia
	reffe_u = LagrangianRefFE(HEX,VectorValue{3,Float64},2)
# time: 2023-06-26 04:06:54 CEST
# mode: julia
	using Gridap.ReferenceFEs
# time: 2023-06-26 04:06:55 CEST
# mode: julia
	reffe_u = LagrangianRefFE(HEX,VectorValue{3,Float64},2)
# time: 2023-06-26 04:07:04 CEST
# mode: julia
	reffe_u = LagrangianRefFE(VectorValue{3,Float64},HEX,2)
# time: 2023-06-26 04:07:17 CEST
# mode: julia
	get_face_dofs(reffe_u)
# time: 2023-06-26 04:09:54 CEST
# mode: julia
	get_face_dofs(reffe_u,0)
# time: 2023-06-26 04:09:58 CEST
# mode: julia
	get_face_dofs(reffe_u,1)
# time: 2023-06-26 04:10:01 CEST
# mode: julia
	get_face_dofs(reffe_u,0)
# time: 2023-06-26 04:10:11 CEST
# mode: julia
	face_dofs = get_face_dofs(reffe_u,0)
# time: 2023-06-26 04:10:29 CEST
# mode: julia
	cell_dof_ids = get_cell_dof_ids(Uu)
# time: 2023-06-26 04:11:03 CEST
# mode: julia
	cell_dof_ids = get_cell_dof_ids(Uu,Γ₀)
# time: 2023-06-26 04:12:04 CEST
# mode: julia
	map(ids->lazy_map(Reindex(ids),face_dofs[1]),cell_dof_ids)
# time: 2023-06-26 04:12:40 CEST
# mode: julia
	map(ids->lazy_map(Reindex(ids),face_dofs[2]),cell_dof_ids)
# time: 2023-06-26 04:13:01 CEST
# mode: julia
	face_dofs = get_face_dofs(reffe_u,1)
# time: 2023-06-26 04:13:04 CEST
# mode: julia
	map(ids->lazy_map(Reindex(ids),face_dofs[1]),cell_dof_ids)
# time: 2023-06-26 04:16:38 CEST
# mode: julia
	face_dofs = get_face_dofs(reffe_u,1)
# time: 2023-06-26 04:38:20 CEST
# mode: julia
	arr[1]
# time: 2023-06-26 04:38:31 CEST
# mode: julia
	arr[1][2,1]
# time: 2023-06-26 04:38:46 CEST
# mode: julia
	nnz(arr[1][2,1])
# time: 2023-06-26 04:38:57 CEST
# mode: julia
	count(arr[1][2,1])
# time: 2023-06-26 04:39:14 CEST
# mode: julia
	count(arr[1][2,1] .> tol)
# time: 2023-06-26 04:39:19 CEST
# mode: julia
	count(arr[1][1,2] .> tol)
# time: 2023-06-26 04:39:23 CEST
# mode: julia
	count(arr[2][1,2] .> tol)
# time: 2023-06-26 04:39:27 CEST
# mode: julia
	count(arr[3][1,2] .> tol)
# time: 2023-06-26 04:39:30 CEST
# mode: julia
	count(arr[4][1,2] .> tol)
# time: 2023-06-26 04:39:35 CEST
# mode: julia
	count(arr[5][1,2] .> tol)
# time: 2023-06-26 04:41:23 CEST
# mode: julia
	count(arr[5][1,2][1:3,1:81] .> tol)
# time: 2023-06-26 04:41:30 CEST
# mode: julia
	count(arr[5][2,1][1:3,1:81] .> tol)
# time: 2023-06-26 04:41:36 CEST
# mode: julia
	count(arr[5][2,1][4:6,1:81] .> tol)
# time: 2023-06-26 04:41:42 CEST
# mode: julia
	count(arr[1][2,1][4:6,1:81] .> tol)
# time: 2023-06-26 04:41:47 CEST
# mode: julia
	count(arr[1][2,1][1:3,1:81] .> tol)
# time: 2023-06-26 04:42:38 CEST
# mode: julia
	count(arr[1][2,1] .> tol)
# time: 2023-06-26 04:42:41 CEST
# mode: julia
	count(arr[2][2,1] .> tol)
# time: 2023-06-26 04:42:44 CEST
# mode: julia
	count(arr[3][2,1] .> tol)
# time: 2023-06-26 04:49:51 CEST
# mode: julia
	get_array(aΩ((u,λ),(v,μ)))
# time: 2023-06-26 04:49:54 CEST
# mode: julia
	get_array(aΩ((u,λ),(v,μ)))[1]
# time: 2023-06-26 04:51:46 CEST
# mode: julia
	Uu.dirichlet_
# time: 2023-06-26 04:51:49 CEST
# mode: julia
	Uu.dirichlet_cells
# time: 2023-06-26 04:51:54 CEST
# mode: julia
	Uλ.dirichlet_cells
# time: 2023-06-26 04:53:52 CEST
# mode: julia
	aΓ((u,λ),(v,μ))
# time: 2023-06-26 05:13:04 CEST
# mode: julia
	1000*1000
# time: 2023-06-26 05:38:21 CEST
# mode: julia
	u
# time: 2023-06-26 05:38:26 CEST
# mode: julia
	u.plus
# time: 2023-06-26 05:38:31 CEST
# mode: julia
	u.minus
# time: 2023-06-26 06:19:23 CEST
# mode: julia
	Γ₃.glue
# time: 2023-06-26 06:19:32 CEST
# mode: julia
	Γ₃.glue.face_to_bgface
# time: 2023-06-26 06:20:30 CEST
# mode: julia
	partial_glues[1]
# time: 2023-06-26 06:20:39 CEST
# mode: julia
	glue = partial_glues[1]
# time: 2023-06-26 06:20:47 CEST
# mode: julia
	glue.n2o_faces_map[3]
# time: 2023-06-26 06:22:32 CEST
# mode: julia
	Γf.glue
# time: 2023-06-26 06:22:50 CEST
# mode: julia
	Γf.trian.glue.face_to_bgface
# time: 2023-06-26 06:23:04 CEST
# mode: julia
	Γ₀.glue.face_to_bgface
# time: 2023-06-26 06:23:32 CEST
# mode: julia
	Γ1.glue.face_to_bgface
# time: 2023-06-26 06:23:40 CEST
# mode: julia
	Γ₁.glue.face_to_bgface
# time: 2023-06-26 06:25:59 CEST
# mode: julia
	Γ₉ .glue.
# time: 2023-06-26 06:26:10 CEST
# mode: julia
	Γ₉.glue.face_to_bgface
# time: 2023-06-26 06:27:17 CEST
# mode: julia
	Γ₄.glue.face_to_bgface
# time: 2023-06-26 06:27:21 CEST
# mode: julia
	Γ3.glue.face_to_bgface
# time: 2023-06-26 06:27:29 CEST
# mode: julia
	Γ₃.glue.face_to_bgface
# time: 2023-06-26 06:32:25 CEST
# mode: julia
	Γc
# time: 2023-06-26 06:32:48 CEST
# mode: julia
	get_cell_coordinates(Γc)
# time: 2023-06-26 06:39:12 CEST
# mode: julia
	Γ₃.glue.face_to_bgface
# time: 2023-06-26 06:43:13 CEST
# mode: julia
	global_n2o_faces
# time: 2023-06-26 06:43:25 CEST
# mode: julia
	glue.o2n_faces_map
# time: 2023-06-26 06:43:38 CEST
# mode: julia
	glue.n2o_faces_map
# time: 2023-06-26 06:44:21 CEST
# mode: julia
	glue.n2o_cell_to_child_id
# time: 2023-06-26 06:44:56 CEST
# mode: julia
	global_child_ids
# time: 2023-06-26 07:27:24 CEST
# mode: julia
	uh
# time: 2023-06-26 07:27:29 CEST
# mode: julia
	uh.free_values
# time: 2023-06-26 07:31:13 CEST
# mode: julia
	A
# time: 2023-06-26 07:34:10 CEST
# mode: julia
	op
# time: 2023-06-26 07:34:16 CEST
# mode: julia
	op.test
# time: 2023-06-26 07:34:22 CEST
# mode: julia
	op.test.spaces
# time: 2023-06-26 07:34:58 CEST
# mode: julia
	A[end,:]
# time: 2023-06-26 08:11:49 CEST
# mode: julia
	global_rrules
# time: 2023-06-26 08:11:59 CEST
# mode: julia
	global_rrules[1]
# time: 2023-06-26 08:12:07 CEST
# mode: julia
	rr = global_rrules[1]
# time: 2023-06-26 08:12:14 CEST
# mode: julia
	get_cell_map(rr)
# time: 2023-06-26 08:12:38 CEST
# mode: julia
	get_cell_map(rr.grid)
# time: 2023-06-26 08:12:42 CEST
# mode: julia
	get_cell_map(rr.ref_grid)
# time: 2023-06-26 08:12:49 CEST
# mode: julia
	get_cell_map(rr.ref_grid)[1]
# time: 2023-06-26 08:13:13 CEST
# mode: julia
	get_cell_map(rr.ref_grid.grid)
# time: 2023-06-26 08:13:40 CEST
# mode: julia
	rr.ref_grid.grid.cell_map[1].values
# time: 2023-06-26 08:19:26 CEST
# mode: julia
	ux
# time: 2023-06-26 08:19:38 CEST
# mode: julia
	ux[1]
# time: 2023-06-26 08:19:44 CEST
# mode: julia
	ux[1][1,1]
# time: 2023-06-26 08:20:29 CEST
# mode: julia
	ux[1][1:4]
# time: 2023-06-26 08:20:38 CEST
# mode: julia
	ux[1][:,1,1:4]
# time: 2023-06-26 08:20:51 CEST
# mode: julia
	ux[1][1,1][:,1,1:4]
# time: 2023-06-26 08:21:11 CEST
# mode: julia
	ux[1][1,1][:,1,1:8]
# time: 2023-06-26 08:22:35 CEST
# mode: julia
	λx[1][1,1]
# time: 2023-06-26 08:22:41 CEST
# mode: julia
	λx[1][2,2]
# time: 2023-06-26 08:22:47 CEST
# mode: julia
	λx[1]
# time: 2023-06-26 08:23:00 CEST
# mode: julia
	λx[1,2]
# time: 2023-06-26 08:23:06 CEST
# mode: julia
	λx[1][1,2]
# time: 2023-06-26 08:24:59 CEST
# mode: julia
	ux[1][1,1][:,1,1:8]
# time: 2023-06-26 08:26:31 CEST
# mode: julia
	λx[1][1,2]
# time: 2023-06-26 08:50:17 CEST
# mode: julia
	partial_glues
# time: 2023-06-26 08:50:21 CEST
# mode: julia
	partial_glues[1]
# time: 2023-06-26 08:50:28 CEST
# mode: julia
	g = partial_glues[1]
# time: 2023-06-26 08:50:34 CEST
# mode: julia
	g.refinement_rules
# time: 2023-06-26 08:56:48 CEST
# mode: julia
	global_rrules_ptrs
# time: 2023-06-26 09:08:35 CEST
# mode: julia
	num_free_dofs(Uu)
# time: 2023-06-26 09:35:18 CEST
# mode: julia
	v = [2]
# time: 2023-06-26 09:35:23 CEST
# mode: julia
	push!(v,1)
# time: 2023-06-26 09:35:27 CEST
# mode: julia
	push!(v,3)
# time: 2023-06-26 09:35:33 CEST
# mode: julia
	exit()
# time: 2023-06-26 09:51:24 CEST
# mode: julia
	num_cell_dims(Γc)
# time: 2023-06-26 09:51:32 CEST
# mode: julia
	num_point_dims(Γc)
# time: 2023-06-26 09:51:35 CEST
# mode: julia
	num_point_dims(Γf)
# time: 2023-06-26 09:51:39 CEST
# mode: julia
	num_dim_dims(Γf)
# time: 2023-06-26 09:51:45 CEST
# mode: julia
	num_cell_dims(Γf)
# time: 2023-06-29 02:35:04 CEST
# mode: julia
	pts.cell_ref_point[1]
# time: 2023-06-29 02:35:20 CEST
# mode: julia
	pts.cell_ref_point[2]
# time: 2023-06-29 02:35:22 CEST
# mode: julia
	pts.cell_ref_point[3]
# time: 2023-06-29 02:35:25 CEST
# mode: julia
	pts.cell_ref_point[4]
# time: 2023-06-29 02:46:09 CEST
# mode: julia
	get_cell_dof_ids(Γf)
# time: 2023-06-29 02:46:14 CEST
# mode: julia
	get_cell_dof_ids(Uu)
# time: 2023-06-29 02:46:33 CEST
# mode: julia
	get_cell_dof_ids(Uλ)
# time: 2023-06-29 02:56:40 CEST
# mode: julia
	glue.n2o_faces_map
# time: 2023-06-29 02:56:51 CEST
# mode: julia
	glue.o2n_faces_map
# time: 2023-06-29 02:57:05 CEST
# mode: julia
	vcat(glue.o2n_faces_map)
# time: 2023-06-29 02:57:08 CEST
# mode: julia
	vcat(glue.o2n_faces_map...)
# time: 2023-06-29 02:57:26 CEST
# mode: julia
	f2c_cell_map = vcat(glue.o2n_faces_map...)
# time: 2023-06-29 02:58:00 CEST
# mode: julia
	lazy_map(Reindex(get_cell_dof_ids(Uu)),f2c_cell_map)
# time: 2023-06-29 02:58:11 CEST
# mode: julia
	get_cell_dof_ids(Uλ)
# time: 2023-06-29 03:22:03 CEST
# mode: julia
	glue.refinement_rules
# time: 2023-06-29 03:29:02 CEST
# mode: julia
	partial_glues
# time: 2023-06-29 03:29:04 CEST
# mode: julia
	partial_glues[1]
# time: 2023-06-29 03:29:33 CEST
# mode: julia
	Γ[1]
# time: 2023-06-29 03:29:39 CEST
# mode: julia
	num_cells(Γ[1])
# time: 2023-06-29 03:32:04 CEST
# mode: julia
	glue
# time: 2023-06-29 03:32:27 CEST
# mode: julia
	get_n2o_reference_coordinate_map(glue)
# time: 2023-06-29 03:37:37 CEST
# mode: julia
	glue
# time: 2023-06-29 03:37:50 CEST
# mode: julia
	glue.n2o_faces_map[3]
# time: 2023-06-29 03:39:25 CEST
# mode: julia
	get_node_coordinates(Γ₉)
# time: 2023-06-29 03:39:33 CEST
# mode: julia
	get_vertex_coordinates(Γ₉)
# time: 2023-06-29 03:39:42 CEST
# mode: julia
	Geometry.get_node_coordinates(Γ₉)
# time: 2023-06-29 03:40:22 CEST
# mode: julia
	coords = Geometry.get_node_coordinates(Γ₉)
# time: 2023-06-29 03:40:29 CEST
# mode: julia
	coords = Geometry.get_cell_node_coordinates(Γ₉)
# time: 2023-06-29 03:40:39 CEST
# mode: julia
	Geometry.get_node_coordinates(Γ₉)
# time: 2023-06-29 03:40:46 CEST
# mode: julia
	get_cell_dof_ids
# time: 2023-06-29 03:41:55 CEST
# mode: julia
	Geometry.get_cell_coordinates(Γ₉)
# time: 2023-06-29 03:57:37 CEST
# mode: julia
	glue.o2n_faces_map
# time: 2023-06-29 04:07:50 CEST
# mode: julia
	get_o2n_faces_map(glue)
# time: 2023-06-29 04:07:57 CEST
# mode: julia
	Adaptivity.get_o2n_faces_map(glue)
# time: 2023-06-29 04:08:17 CEST
# mode: julia
	glue.o2n_faces_map
# time: 2023-06-29 04:12:11 CEST
# mode: julia
	partial_glues[1]
# time: 2023-06-29 04:12:22 CEST
# mode: julia
	partial_glues[1].n2o_faces_map[3]
# time: 2023-06-29 04:12:28 CEST
# mode: julia
	partial_glues[2].n2o_faces_map[3]
# time: 2023-06-29 04:12:32 CEST
# mode: julia
	partial_glues[3].n2o_faces_map[3]
# time: 2023-06-29 04:12:40 CEST
# mode: julia
	partial_glues[4].n2o_faces_map[3]
# time: 2023-06-29 04:14:26 CEST
# mode: julia
	Γ[4].glue.face_to_bgface
# time: 2023-06-29 04:14:39 CEST
# mode: julia
	boundary_faces
# time: 2023-06-29 05:09:36 CEST
# mode: julia
	CellData.DomainStyle(fun_ue)
# time: 2023-06-29 05:13:57 CEST
# mode: julia
	get_cell_dof_ids(V_sol)
# time: 2023-06-29 05:15:26 CEST
# mode: julia
	zero_free_dofs(U_sol)
# time: 2023-06-29 05:15:33 CEST
# mode: julia
	zero_free_dof_values(U_sol)
# time: 2023-06-29 05:15:48 CEST
# mode: julia
	zero_free_values(U_sol)
# time: 2023-06-29 06:42:53 CEST
# mode: julia
	xC
# time: 2023-06-29 06:46:56 CEST
# mode: julia
	intrf₀
# time: 2023-06-29 07:08:28 CEST
# mode: julia
	ue_c₀(PointValue(0.0))
# time: 2023-06-29 07:08:33 CEST
# mode: julia
	ue_c₀(Point(0.0))
# time: 2023-06-29 07:09:09 CEST
# mode: julia
	fun_ue₀(1.0)
# time: 2023-06-29 07:09:16 CEST
# mode: julia
	fun_ue₀(Point(1.0)
	)
# time: 2023-06-29 07:09:22 CEST
# mode: julia
	fun_ue₀(Point(.0)
	)
# time: 2023-06-29 07:09:29 CEST
# mode: julia
	fun_ue₀(Point(0.0))
# time: 2023-06-29 07:09:33 CEST
# mode: julia
	fun_ue₀(Point(0.1))
# time: 2023-06-29 07:09:35 CEST
# mode: julia
	fun_ue₀(Point(0.2))
# time: 2023-06-29 07:09:37 CEST
# mode: julia
	fun_ue₀(Point(0.3))
# time: 2023-06-29 07:09:40 CEST
# mode: julia
	fun_ue₀(Point(0.5))
# time: 2023-07-03 03:32:28 CEST
# mode: julia
	first(contr.dict)
# time: 2023-07-03 03:32:58 CEST
# mode: julia
	tt
# time: 2023-07-03 03:36:07 CEST
# mode: julia
	get_cell_dof_ids(Vλ,Ψ₉)
# time: 2023-07-03 03:36:19 CEST
# mode: julia
	dof_ids = get_cell_dof_ids(Vλ,Ψ₉)
# time: 2023-07-03 03:44:27 CEST
# mode: julia
	glue
# time: 2023-07-03 03:53:48 CEST
# mode: julia
	tt
# time: 2023-07-03 03:53:54 CEST
# mode: julia
	tt.parent
# time: 2023-07-03 03:54:06 CEST
# mode: julia
	tt.parent.dtrian
# time: 2023-07-03 03:54:13 CEST
# mode: julia
	tt.parent.dtrian.minus
# time: 2023-07-03 03:54:18 CEST
# mode: julia
	tt.parent.rtrian
# time: 2023-07-03 04:11:31 CEST
# mode: julia
	arr
# time: 2023-07-03 04:11:44 CEST
# mode: julia
	typeof(arr)
# time: 2023-07-03 04:11:58 CEST
# mode: julia
	arr[1]
# time: 2023-07-03 04:12:52 CEST
# mode: julia
	tt
# time: 2023-07-03 04:12:57 CEST
# mode: julia
	tt.parent
# time: 2023-07-03 04:13:09 CEST
# mode: julia
	tt.parent.rtrian
# time: 2023-07-03 04:16:00 CEST
# mode: julia
	move_contributions(tt,arr)
# time: 2023-07-03 04:16:08 CEST
# mode: julia
	move_contributions(arr,tt)
# time: 2023-07-03 04:16:20 CEST
# mode: julia
	move_contributions(arr,tt)[1]
# time: 2023-07-03 04:16:36 CEST
# mode: julia
	typeof(move_contributions(arr,tt))
# time: 2023-07-03 04:16:49 CEST
# mode: julia
	move_contributions(arr,tt)[1][1]
# time: 2023-07-03 04:17:25 CEST
# mode: julia
	length(move_contributions(arr,tt)[1])
# time: 2023-07-03 04:18:13 CEST
# mode: julia
	tt.parent
# time: 2023-07-03 04:18:18 CEST
# mode: julia
	tt.parent.dtrian
# time: 2023-07-03 04:19:46 CEST
# mode: julia
	skeleton_glue = get_glue(tt.parent.dtrian,Val(2))
# time: 2023-07-03 04:20:35 CEST
# mode: julia
	skeleton_glue.minus
# time: 2023-07-03 04:20:50 CEST
# mode: julia
	skeleton_glue.minus.tface_to_mface
# time: 2023-07-03 04:21:13 CEST
# mode: julia
	tt.cell_to_parent_cell
# time: 2023-07-03 04:21:59 CEST
# mode: julia
	skeleton_glue.minus.tface_to_mface[tt.cell_to_parent_cell]
# time: 2023-07-03 04:22:27 CEST
# mode: julia
	num_cells(Γ₉)
# time: 2023-07-03 04:24:16 CEST
# mode: julia
	tt.parent.rtrian.trian.glue.face_to_bgface
# time: 2023-07-03 04:25:04 CEST
# mode: julia
	num_faces(model)
# time: 2023-07-03 04:25:09 CEST
# mode: julia
	num_faces(model,0)
# time: 2023-07-03 04:25:12 CEST
# mode: julia
	num_faces(model,1)
# time: 2023-07-03 04:25:15 CEST
# mode: julia
	num_faces(model,2)
# time: 2023-07-03 04:25:22 CEST
# mode: julia
	num_faces(model,3)
# time: 2023-07-03 04:27:31 CEST
# mode: julia
	num_cells(Γf_minus)
# time: 2023-07-03 04:27:39 CEST
# mode: julia
	num_cells(Γf_plus)
# time: 2023-07-03 04:32:52 CEST
# mode: julia
	using Gridap
# time: 2023-07-03 04:34:13 CEST
# mode: julia
	exit()
# time: 2023-07-03 04:47:03 CEST
# mode: julia
	typeof(arr_plus)
# time: 2023-07-03 04:47:14 CEST
# mode: julia
	typeof(arr_minus)
# time: 2023-07-03 04:47:48 CEST
# mode: julia
	arr
# time: 2023-07-03 04:47:55 CEST
# mode: julia
	arr[1]
# time: 2023-07-03 04:50:05 CEST
# mode: julia
	arr[1][2]
# time: 2023-07-03 04:56:35 CEST
# mode: julia
	arr_plus[1]
# time: 2023-07-03 06:16:47 CEST
# mode: julia
	cf_f
# time: 2023-07-03 06:19:09 CEST
# mode: julia
	μ_f.cell_basis
# time: 2023-07-03 06:19:29 CEST
# mode: julia
	μ_f.cell_basis[1]
# time: 2023-07-03 06:19:48 CEST
# mode: julia
	ue_f.cell_field
# time: 2023-07-03 06:22:12 CEST
# mode: julia
	GenericField
# time: 2023-07-03 06:22:18 CEST
# mode: julia
	Gridap.Fields.GenericField
# time: 2023-07-03 06:22:25 CEST
# mode: julia
	Gridap.Fields.GenericField(_get_y)
# time: 2023-07-03 06:22:41 CEST
# mode: julia
	Gridap.Fields.GenericField(_get_y)(Point(0.0))
# time: 2023-07-03 06:22:48 CEST
# mode: julia
	Gridap.Fields.GenericField(_get_y)(Point(0.0,1.0))
# time: 2023-07-03 06:24:36 CEST
# mode: julia
	ue_c.cell_field
# time: 2023-07-03 06:28:58 CEST
# mode: julia
	Array(undef,(2,3))
# time: 2023-07-03 06:29:08 CEST
# mode: julia
	Matrix(undef,(2,3))
# time: 2023-07-03 06:29:16 CEST
# mode: julia
	Array{Float64}(undef,(2,3))
# time: 2023-07-03 06:29:46 CEST
# mode: julia
	x = MatrixBlock(zeros(3,3),fill(true,(3,3)))
# time: 2023-07-03 06:29:55 CEST
# mode: julia
	x = Gridap.Fields.MatrixBlock(zeros(3,3),fill(true,(3,3)))
# time: 2023-07-03 06:30:07 CEST
# mode: julia
	x = Gridap.Fields.ArrayBlock(zeros(3,3),fill(true,(3,3)))
# time: 2023-07-03 06:30:14 CEST
# mode: julia
	eltype(x)
# time: 2023-07-03 06:33:45 CEST
# mode: julia
	x
# time: 2023-07-03 06:33:49 CEST
# mode: julia
	num_dims(x)
# time: 2023-07-03 06:33:57 CEST
# mode: julia
	dims(x)
# time: 2023-07-03 07:15:07 CEST
# mode: julia
	Γ
# time: 2023-07-03 07:15:27 CEST
# mode: julia
	boundary_faces
# time: 2023-07-03 07:15:35 CEST
# mode: julia
	map(num_cells,Γ)
# time: 2023-07-03 07:16:16 CEST
# mode: julia
	Γ[1].glue.face_to_bgface
# time: 2023-07-03 07:18:31 CEST
# mode: julia
	num_faces(model,1)
# time: 2023-07-03 07:18:35 CEST
# mode: julia
	num_faces(model,2)
# time: 2023-07-03 07:40:46 CEST
# mode: julia
	partial_glues
# time: 2023-07-03 07:42:24 CEST
# mode: julia
	zip(Γ,pglues)
# time: 2023-07-03 07:42:28 CEST
# mode: julia
	zip(Γ,pglues)[1]
# time: 2023-07-03 07:42:45 CEST
# mode: julia
	first(zip(Γ,pglues))
# time: 2023-07-03 07:45:32 CEST
# mode: julia
	Γ
# time: 2023-07-03 07:45:37 CEST
# mode: julia
	Γ[1]
# time: 2023-07-03 07:45:45 CEST
# mode: julia
	Γ[1].glue.bgface_to_lcell
# time: 2023-07-03 07:45:52 CEST
# mode: julia
	Γ[1].glue.face_to_bgface
# time: 2023-07-03 07:46:38 CEST
# mode: julia
	Γ[2].glue.face_to_bgface
# time: 2023-07-03 07:47:28 CEST
# mode: julia
	pglues
# time: 2023-07-03 07:47:32 CEST
# mode: julia
	pglues[1]
# time: 2023-07-03 07:47:47 CEST
# mode: julia
	totalColumns
# time: 2023-07-04 03:08:28 CEST
# mode: pkg
	add GridapPardiso
# time: 2023-07-04 03:09:06 CEST
# mode: pkg
	build GridPardiso
# time: 2023-07-04 03:09:15 CEST
# mode: pkg
	status
# time: 2023-07-04 03:09:33 CEST
# mode: julia
	exit()
# time: 2023-07-04 03:17:37 CEST
# mode: julia
	source /opt/intel/mkl/bin/mklvars.sh intel64
# time: 2023-07-04 03:17:43 CEST
# mode: julia
	exit()
# time: 2023-07-04 03:48:33 CEST
# mode: julia
	source /home/deck/.local/podman/lib/containers/storage/overlay/314eafab7e0b9fffc8bb3f907fbf5c59be5afcc626f6729172bb8410ccd4fa35/diff/opt/intel/oneapi/mkl/2023.1.0/env/vars.sh intel64
# time: 2023-07-04 04:02:54 CEST
# mode: julia
	include("MDC_SL_Kinematic3D_Inrtf_SuperElement_H1_ambEdge")
# time: 2023-07-04 04:04:21 CEST
# mode: julia
	include("MDC_SL_Kinematic3D_Inrtf_SuperElement_H1_ambEdge.jl")
# time: 2023-07-04 04:22:25 CEST
# mode: julia
	cd model
# time: 2023-07-04 04:22:30 CEST
# mode: julia
	exit()
# time: 2023-07-04 04:23:24 CEST
# mode: julia
	include("MDC_SL_Kinematic3D_Inrtf_SuperElement_H1_ambEdge.jl")
# time: 2023-07-04 04:24:12 CEST
# mode: julia
	push!(LOAD_PATH, pwd()*"/src")
# time: 2023-07-04 04:24:21 CEST
# mode: julia
	include("MDC_SL_Kinematic3D_Inrtf_SuperElement_H1_ambEdge.jl")
# time: 2023-07-04 04:25:05 CEST
# mode: julia
	push!(LOAD_PATH, pwd()*"/src")
# time: 2023-07-04 04:25:34 CEST
# mode: julia
	exit()
# time: 2023-07-04 04:38:49 CEST
# mode: julia
	run(`ls`)
# time: 2023-07-04 04:39:10 CEST
# mode: julia
	run(`source /home/deck/vars.sh intel64`)
# time: 2023-07-04 04:39:25 CEST
# mode: julia
	run(`source ~/home/deck/vars.sh intel64`)
# time: 2023-07-04 04:40:29 CEST
# mode: julia
	include("MDC_SL_Kinematic3D_Inrtf_SuperElement_H1_ambEdge.jl")
# time: 2023-07-04 04:40:46 CEST
# mode: julia
	exit()
# time: 2023-07-04 04:41:46 CEST
# mode: pkg
	up
# time: 2023-07-04 04:43:20 CEST
# mode: pkg
	]build
# time: 2023-07-04 04:43:27 CEST
# mode: pkg
	build
# time: 2023-07-04 05:02:17 CEST
# mode: julia
	λh
# time: 2023-07-04 05:09:07 CEST
# mode: julia
	λh_Λc_H1
# time: 2023-07-04 05:09:13 CEST
# mode: julia
	λh_Λc_H1.free_values
