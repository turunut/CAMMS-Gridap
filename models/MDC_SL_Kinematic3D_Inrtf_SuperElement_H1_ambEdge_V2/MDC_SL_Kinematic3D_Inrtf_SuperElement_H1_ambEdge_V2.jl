push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/Materials")

using Gridap
using Gridap.Geometry
using Gridap.Adaptivity
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.CellData
using modCT
using modModel
using modSubroutines
using modInterface

using GridapGmsh

using FillArrays

prblName = "MDC_SL_Kinematic3D_Inrtf_SuperElement_H1_ambEdge_V2"
projFldr = pwd()

order  = 2
degree = 2*order

#--------------------------------------------------
# Definim el model volumetric

domain = (0,4,0,4,-0.5,0.5)
partition = (8,8,4)
model = CartesianDiscreteModel(domain,partition)

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"faceY0",[23])
add_tag_from_tags!(labels,"faceY1",[24])
add_tag_from_tags!(labels,"faceX0",[25])
add_tag_from_tags!(labels,"faceX1",[26])


#--------------------------------------------------
# Definim les propietats del model solid 

modlType = Solid()

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(72000, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------
# Definim la triangulacio del model solid

Ω  = Triangulation(model)
dΩ = Measure(Ω,degree)


#--------------------------------------------------
# Definim el Cell Field de CT del model solid

dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1)

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#---------------------------
# Definim les propietats de cada layer per fer un model Shell equivalent

modlType_2D = PlaneStress()

ct1_2D = modModel.computeCT(modlType_2D, CT1, true)
ct2_2D = modModel.computeCT(modlType_2D, CT2, true)

CTs_2D = hcat(ct1_2D)

CTf_2D = get_CT_CellField(modlType_2D, CTs_2D, tags, Ω)


#--------------------------------------------------
# Acotem les triangulacions de les interficies

# Totes juntes
Γ₉ = BoundaryTriangulation(model,tags=["tag_23","tag_26","tag_24","tag_25"])

# Cada una de les cares
Γ₀ = BoundaryTriangulation(model,tags=["tag_23"])
Γ₁ = BoundaryTriangulation(model,tags=["tag_26"])
Γ₂ = BoundaryTriangulation(model,tags=["tag_24"])
Γ₃ = BoundaryTriangulation(model,tags=["tag_25"])


#--------------------------------------------------
# Definim els models 1x1 fine i coarse del conjunt de les interficies
# Creem una llista de coordenades en forma de integer dels CG dels elements de les cares per
# despres poder poder filtrarlos per la seva posicio 

Dc      = num_cell_dims(model)
topo    = get_grid_topology(model)
f2n_map = Geometry.get_faces(topo,Dc-1,0)
coords  = Geometry.get_vertex_coordinates(topo)
face_coords = map(Reindex(coords),f2n_map)

tol = 1.e-3
com_coords = lazy_map(N->sum(N)/length(N),face_coords) # center of mass
int_coords = map(N->VectorValue(Int(floor(N[1]/tol)),Int(floor(N[2]/tol)),Int(floor(N[3]/tol))),com_coords)

# Definim leix i posicio en aquest de cada una de les interficies

z_coord(x) = x[end]; zf = CellField(z_coord,Ω)

axis_id₀ = 2; face_pos₀ = minimum(lazy_map(c->c[axis_id₀],int_coords))

axis_id₁ = 1; face_pos₁ = maximum(lazy_map(c->c[axis_id₁],int_coords))

axis_id₂ = 2; face_pos₂ = maximum(lazy_map(c->c[axis_id₂],int_coords))

axis_id₃ = 1; face_pos₃ = minimum(lazy_map(c->c[axis_id₃],int_coords))

#glue₉ = create_interface_global([Γ₀,Γ₁,Γ₂,Γ₃],
#                                int_coords, 
#                                [axis_id₀, axis_id₁, axis_id₂,axis_id₃],
#                                [face_pos₀,face_pos₁,face_pos₂,face_pos₃],
#                                [false,false,true,true])

# Definim el model global 1x1 i la glue corresponent de totes les interficies

Γ = [Γ₀,Γ₁,Γ₂,Γ₃]
axis_id = [axis_id₀, axis_id₁, axis_id₂,axis_id₃]
axis_int_coord = [face_pos₀,face_pos₁,face_pos₂,face_pos₃]
inv_list = [false,false,true,true]

global partial_glues = AdaptivityGlue[]

global totalColumns = 0
global n_refs = 0 # Divisio vertical
for iintrf in eachindex(Γ)
  n2o_faces, child_ids, o2n_faces = get_comp_glue(Γ[iintrf], int_coords, axis_id[iintrf], axis_int_coord[iintrf], inv_list[iintrf])
  n_coarse = length(o2n_faces)
  global n_refs = length(o2n_faces[1])

  n2o_faces[3] .+= totalColumns

  ref_grid_domain = (!inv_list[iintrf]) ? (0,1,0,1) : (1,0,0,1)
  ref_grid = Geometry.UnstructuredDiscreteModel(CartesianDiscreteModel(ref_grid_domain,(1,n_refs)))
  rr = Adaptivity.RefinementRule(Adaptivity.GenericRefinement(),QUAD,ref_grid)
  rrules = Fill(rr,n_coarse)
  glue = AdaptivityGlue(n2o_faces,child_ids,rrules) # From coarse to fine
  push!(partial_glues, glue)
  global totalColumns += n_coarse
end

nF = num_cells(Γ₉)
global_n2o_faces = zeros(Int,nF)
global_child_ids = zeros(Int,nF)
global_rrules_data = Vector{typeof(first(first(partial_glues).refinement_rules))}(undef,length(partial_glues))
global_rrules_ptrs = zeros(Int32,totalColumns)

global nCoarse = 1
boundary_faces = Γ₉.glue.face_to_bgface
for (i,glue,Γi) in zip(1:length(partial_glues),partial_glues,Γ)
  interface_faces = Γi.glue.face_to_bgface
  face_map = map(face->findfirst(x->x==face,boundary_faces),interface_faces)
  global_n2o_faces[face_map] .= glue.n2o_faces_map[3]
  global_child_ids[face_map] .= glue.n2o_cell_to_child_id
  global_rrules_data[i] = glue.refinement_rules.value

  n_coarse_i = num_cells(Γi) ÷ n_refs
  global_rrules_ptrs[nCoarse:nCoarse+n_coarse_i-1] .= i
  global nCoarse += n_coarse_i
end
@assert all(global_n2o_faces .!= 0)
@assert all(global_child_ids .!= 0)

global_rrules = CompressedArray(global_rrules_data,global_rrules_ptrs)
glue = AdaptivityGlue([Int[],Int[],global_n2o_faces],global_child_ids,global_rrules) # From coarse to fine

cface_model = CartesianDiscreteModel((0,1,0,1),(totalColumns,1),isperiodic=(true,false))
Γc  = Triangulation(cface_model)
Γf  = Gridap.Adaptivity.GluedTriangulation(Γ₉,Γc,glue)
Ψf  = EdgeTriangulation(Γf,["tag_17","tag_18","tag_19","tag_20"])

# Definim 

Γ₀_glue = Gridap.Adaptivity.GluedTriangulation(Γ₀,Γc,partial_glues[1])
Γ₁_glue = Gridap.Adaptivity.GluedTriangulation(Γ₁,Γc,partial_glues[2])
Γ₂_glue = Gridap.Adaptivity.GluedTriangulation(Γ₂,Γc,partial_glues[3])
Γ₃_glue = Gridap.Adaptivity.GluedTriangulation(Γ₃,Γc,partial_glues[4])

intrf₀  = Intrf_Kinematic3DV2(Ω, Γf, Ψ₀, Γ₀_glue, CTf_2D[1], degree)

intrf₁  = Intrf_Kinematic3DV2(Ω, Γf, Ψ₁, Γ₁_glue, CTf_2D[1], degree)

intrf₂  = Intrf_Kinematic3DV2(Ω, Γf, Ψ₂, Γ₂_glue, CTf_2D[1], degree)

intrf₃  = Intrf_Kinematic3DV2(Ω, Γf, Ψ₃, Γ₃_glue, CTf_2D[1], degree)


#--------------------------------------------------
# Definim els FE del model volumetric i dels multiplicadors de lagrange (un per cada columna delements) 
# Definim el model 3D
reffe_u  = ReferenceFE(lagrangian,VectorValue{3,Float64},order)#;space=:S)

Vu = TestFESpace(model,reffe_u;conformity=:H1)
Uu = TrialFESpace(Vu)

# Definim l'espai de l'acoplament
dofs  = 3
reffe = ReferenceFE(lagrangian,VectorValue{dofs,Float64},(1,0))
Vλ = FESpace(Γc,reffe)
#reffe = ReferenceFE(lagrangian,VectorValue{dofs,Float64},0)
#Vλ = FESpace(Γc,reffe;conformity=:L2)
Uλ = TrialFESpace(Vλ)

# Ajuntem els dos espais anteriors
V = MultiFieldFESpace([Vu,Vλ])
U = MultiFieldFESpace([Uu,Uλ])


#--------------------------------------------------
# Model linea fine amb periodicitat als seus extrems

modelΛf = CartesianDiscreteModel((0,1),totalColumns;isperiodic=Tuple(true))
Λf = Triangulation(modelΛf)

dofs = 3
reffeΛ = ReferenceFE(lagrangian,VectorValue{dofs,Float64},order)
VΛf = FESpace(Λf,reffeΛ,conformity=:H1)
UΛf = TrialFESpace(VΛf)


#--------------------------------------------------
interfaces = [intrf₀, intrf₁, intrf₂, intrf₃]
crs_dsc = [1,1,1,1]

# Defineixo el model lineal coarse amb peridicitat als extrems
domainΛc = (0, 1); partitionΛc = (sum(crs_dsc))
modelΛc  = CartesianDiscreteModel(domainΛc, partitionΛc; isperiodic=Tuple(true))
Λc  = Triangulation(modelΛc)
# i el seu FEspace
VΛc = FESpace(modelΛc,reffeΛ,conformity=:H1); UΛc = TrialFESpace(VΛc)


#--------------------------------------------------
# Dic quin DOF vui activar
active_DOF = 1

# Projeccio de la funcio en el model linea al llarg de leix Y del model 2D coarse de les interfaces
_get_y(x) = VectorValue(x[2])
function π_Λe_Γc(f::CellField, Γc::Triangulation)
  _data = CellData.get_data(f)
  _cellmap = Fill(Broadcasting(_get_y),length(_data))
  data = lazy_map(∘,_data,_cellmap) # data = lazy_map(Broadcasting(∘),_data,_cellmap) 
  return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
end

# Defineixo la funcio activant el DOF
xC = zero_free_values(UΛc);
if active_DOF != 0; xC[active_DOF] = 1.0; end
fun_uC = FEFunction(UΛc,xC)

# Interpolo la funcio del model lineal coarse al model lineal fine
uC_intrp = Interpolable(fun_uC)
# Portem la funcio de la linea coarse al model fine
fun_ue = interpolate(uC_intrp,UΛf)
fun_ue(Point(0.0))

ue_c = π_Λe_Γc(fun_ue,Γc)



struct InterfaceProblem
  Ω  ::Triangulation
  Γf ::Triangulation
  Γc ::Triangulation
  Λ  ::Triangulation
  Ψ  ::Triangulation
  dΨ ::Measure
  intrf :: Vector{Intrf_Kinematic3DV2}
end

dΨ = Measure(Ψf,degree)
prob = InterfaceProblem(Ω, Γf, Γc, Λf, Ψf, dΨ, [intrf₀, intrf₁, intrf₂, intrf₃])

function func(problem::InterfaceProblem,u,v,λ,μ)

  λ_f = change_domain(λ,problem.Γf,ReferenceDomain())
  μ_f = change_domain(μ,problem.Γf,ReferenceDomain())

  contr = DomainContribution()
  for intrf in problem.intrf
    dΓfi = intrf.dΓi
    ci = ∫( (λ_f⋅v) + (μ_f⋅u) )*dΓfi
    contr += ci
  end
  return contr
end

function get_plus_component(x::Gridap.Fields.MatrixBlock)
  array = Matrix{Matrix{Float64}}(undef,size(x))
  for i in eachindex(x.touched)
    if x.touched[i]
      array[i] = x.array[i].array[1,1]
    end
  end
  return Gridap.Fields.ArrayBlock(array,x.touched)
end

function get_minus_component(x::Gridap.Fields.MatrixBlock)
  array = Matrix{Matrix{Float64}}(undef,size(x))
  for i in eachindex(x.touched)
    if x.touched[i]
      array[i] = x.array[i].array[2,2]
    end
  end
  return Gridap.Fields.ArrayBlock(array,x.touched)
end

function FESpaces.get_cell_fe_data(fun,f,ttrian::Geometry.TriangulationView{Dc,Dp,<:Adaptivity.GluedTriangulation}) where {Dc,Dp}
  fe_data = FESpaces.get_cell_fe_data(fun,f,ttrian.parent)
  return lazy_map(Reindex(fe_data),ttrian.cell_to_parent_cell)
end

function func_edge(problem::InterfaceProblem,u,v,λ,μ)
  Γf = problem.Γf
  λ_f = change_domain(λ,Γf,ReferenceDomain())
  μ_f = change_domain(μ,Γf,ReferenceDomain())

  dΨ = problem.dΨ
  Ψ  = problem.Ψ
  arr = get_array(∫( mean((λ_f⋅v) + (μ_f⋅u)) )*dΨ)

  skeleton_glue  = get_glue(Ψ.parent.dtrian,Val(2))
  minus_Γf_faces = skeleton_glue.minus.tface_to_mface[Ψ.cell_to_parent_cell]
  plus_Γf_faces  = skeleton_glue.plus.tface_to_mface[Ψ.cell_to_parent_cell]

  Γf_plus   = Geometry.TriangulationView(Γf,plus_Γf_faces)
  Γf_minus  = Geometry.TriangulationView(Γf,minus_Γf_faces)
  arr_plus  = lazy_map(get_plus_component,arr)
  arr_minus = lazy_map(get_minus_component,arr)

  contr = DomainContribution()
  CellData.add_contribution!(contr,Γf_plus,arr_plus)
  CellData.add_contribution!(contr,Γf_minus,arr_minus)

  return contr
end

function func_edge_vector(problem::InterfaceProblem,v,μ,ue_c)
  Γf = problem.Γf
  cf = μ⋅ue_c
  cf_f = change_domain(cf,Γf,ReferenceDomain())

  dΨ = problem.dΨ
  Ψ  = problem.Ψ
  arr = get_array(∫( mean(cf_f) )*dΨ)

  skeleton_glue  = get_glue(Ψ.parent.dtrian,Val(2))
  minus_Γf_faces = skeleton_glue.minus.tface_to_mface[Ψ.cell_to_parent_cell]
  plus_Γf_faces  = skeleton_glue.plus.tface_to_mface[Ψ.cell_to_parent_cell]

  Γf_plus   = Geometry.TriangulationView(Γf,plus_Γf_faces)
  Γf_minus  = Geometry.TriangulationView(Γf,minus_Γf_faces)
  arr_plus  = lazy_map(get_plus_component,arr)
  arr_minus = lazy_map(get_minus_component,arr)

  contr = DomainContribution()
  CellData.add_contribution!(contr,Γf_plus,arr_plus)
  CellData.add_contribution!(contr,Γf_minus,arr_minus)

  return contr
end


aΩ((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ
aΓ((u,λ),(v,μ)) = func(prob,u,v,λ,μ)
aΨ((u,λ),(v,μ)) = func_edge(prob,u,v,λ,μ)

a((u,λ),(v,μ)) =  aΩ((u,λ),(v,μ)) + 
                   aΓ((u,λ),(v,μ)) + 
                   aΨ((u,λ),(v,μ))                                       

A = assemble_matrix(a,U,V)

line = A.m-2
for i in 1:7904
  if abs( A[line,i] ) > 0.0
    println( i,"   ",A[line,i] )
  end
end

f = VectorValue(0.0,0.0,0.0)
g₀ = VectorValue(1.0,0.0,0.0)
g₁ = VectorValue(0.0,0.0,0.0)
g₂ = VectorValue(0.0,0.0,0.0)
g₃ = VectorValue(0.0,0.0,0.0)

#dΓf = Measure(Γf, degree)
#l((v,μ,β₀,β₁,β₂,β₃)) = ∫(v⋅f)*dΩ + 
#                       ∫(μ⋅ue_c)*dΓf + 
#                       ∫(β₀⋅g₀)*dΨ₀ + 
#                       ∫(β₁⋅g₁)*dΨ₁ + 
#                       ∫(β₂⋅g₂)*dΨ₂ + 
#                       ∫(β₃⋅g₃)*dΨ₃

dΓf = Measure(Γf, degree)
l((v,μ)) = ∫(v⋅f)*dΩ + 
           ∫(μ⋅ue_c)*dΓf + 
           func_edge_vector(prob,v,μ,ue_c)
             
b = assemble_vector(l,V)

b[(Uu.nfree+1):end]
b[end-2:end]
b[end-11:end-9]

v,μ = get_fe_basis(V);

cf = μ⋅ue_c
cf_f = mean(change_domain(cf,Γf,ReferenceDomain()))

pts = get_cell_points(dΨ)
cf_f(pts)

#--------------------------------------------------
# Multiplicadors

#op = AffineFEOperator(a,l,U,V)
op = AffineFEOperator(U,V,A,b)

ls = LUSolver()
solver = LinearFESolver(ls)

xh = solve(op);
uh, λh, αh = xh;

println(get_free_dof_values(αh)[1])
for i in 1:3:(3*partition[1]+1)
  println(get_free_dof_values(λh)[i])
end

#x = get_free_dof_values(xh)
#xλ = Gridap.MultiField.restrict_to_field(U,x,2)

#λ₀ = xh.single_fe_functions[2]
#λ₁ = xh.single_fe_functions[3]
#λ₂ = xh.single_fe_functions[4]
#λ₃ = xh.single_fe_functions[5]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])
