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

prblName = "MDC_SL_Kinematic3D_Inrtf_SuperElement_H1_verificacio"
projFldr = pwd()

order  = 1
degree = 2*order

#--------------------------------------------------
# Definim el model volumetric

domain = (0,1,0,1,0,1)
partition = (2,2,1)
model = CartesianDiscreteModel(domain,partition)

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"faceY0",[23])
add_tag_from_tags!(labels,"faceY1",[24])
add_tag_from_tags!(labels,"faceX0",[25])
add_tag_from_tags!(labels,"faceX1",[26])

#--------------------------------------------------
# Definim la triangulacio del model solid

Ω  = Triangulation(model)
dΩ = Measure(Ω,degree)

#--------------------------------------------------
# Acotem les triangulacions de les interficies

# Totes juntes
Γ₉ = BoundaryTriangulation(model,tags=["tag_23","tag_26","tag_24","tag_25"])

# Cada una de les cares
Γ₀ = BoundaryTriangulation(model,tags=["tag_23"])
Γ₁ = BoundaryTriangulation(model,tags=["tag_26"])
Γ₂ = BoundaryTriangulation(model,tags=["tag_24"])
Γ₃ = BoundaryTriangulation(model,tags=["tag_25"])

# Model lineal al llarg de lespessor
Ψ = EdgeTriangulation(model,["tag_17"])


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

# Definim el model global 1x1 i la glue corresponent de totes les interficies

Γ = [Γ₀,Γ₁,Γ₂,Γ₃]
axis_id = [axis_id₀, axis_id₁, axis_id₂,axis_id₃]
axis_int_coord = [face_pos₀,face_pos₁,face_pos₂,face_pos₃]
inv_list = [false,false,true,true]


rr = Adaptivity.RefinementRule(Adaptivity.GenericRefinement(),QUAD,Geometry.UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(1,n_refs))))
rri = Adaptivity.RefinementRule(Adaptivity.GenericRefinement(),QUAD,Geometry.UnstructuredDiscreteModel(CartesianDiscreteModel((1,0,0,1),(1,n_refs))))

global partial_glues = Vector{AdaptivityGlue}(undef,length(Γ))

global totalColumns = 0
global n_refs = 0 # Divisio vertical
for iintrf in eachindex(Γ)
  n2o_faces, child_ids, o2n_faces = get_comp_glue(Γ[iintrf], int_coords, axis_id[iintrf], axis_int_coord[iintrf], inv_list[iintrf])

  n_coarse = length(o2n_faces)
  global n_refs = length(o2n_faces[1])

  n2o_faces[3] .+= totalColumns

  #ref_grid_domain = (!inv_list[iintrf]) ? (0,1,0,1) : (1,0,0,1)
  #ref_grid = Geometry.UnstructuredDiscreteModel(CartesianDiscreteModel(ref_grid_domain,(1,n_refs)))
  #rr = Adaptivity.RefinementRule(Adaptivity.GenericRefinement(),QUAD,ref_grid)

  _rr = (!inv_list[iintrf]) ? rr : rri
  rrules = Fill(_rr,n_coarse)
  glue = AdaptivityGlue(n2o_faces,child_ids,rrules) # From coarse to fine
  partial_glues[iintrf] = glue
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

  n_coarse_i = num_cells(Γi)÷n_refs
  global_rrules_ptrs[nCoarse:nCoarse+n_coarse_i-1] .= i
  global nCoarse += n_coarse_i
end
@assert all(global_n2o_faces .!= 0)
@assert all(global_child_ids .!= 0)

global_rrules = CompressedArray(global_rrules_data,global_rrules_ptrs)
glue = AdaptivityGlue([Int[],Int[],global_n2o_faces],global_child_ids,global_rrules) # From coarse to fine

Geometry.get_cell_coordinates(Γ₉)
glue.n2o_faces_map[3]

rr_fine = Adaptivity.get_new_cell_refinement_rules(glue)
rr_coarse = glue.refinement_rules


cface_model = CartesianDiscreteModel((0,1,0,1),(totalColumns,1),isperiodic=(true,false))
Γc  = Triangulation(cface_model)
Γf  = Gridap.Adaptivity.GluedTriangulation(Γ₉,Γc,glue)

# Definim 

Γ₀_glue = Gridap.Adaptivity.GluedTriangulation(Γ₀,Γc,partial_glues[1])
Γ₁_glue = Gridap.Adaptivity.GluedTriangulation(Γ₁,Γc,partial_glues[2])
Γ₂_glue = Gridap.Adaptivity.GluedTriangulation(Γ₂,Γc,partial_glues[3])
Γ₃_glue = Gridap.Adaptivity.GluedTriangulation(Γ₃,Γc,partial_glues[4])

intrf₀  = Intrf_Kinematic3DV2(Ω, Γf, Ψ, Γ₀_glue, CTf_2D[1], degree)

intrf₁  = Intrf_Kinematic3DV2(Ω, Γf, Ψ, Γ₁_glue, CTf_2D[1], degree)

intrf₂  = Intrf_Kinematic3DV2(Ω, Γf, Ψ, Γ₂_glue, CTf_2D[1], degree)

intrf₃  = Intrf_Kinematic3DV2(Ω, Γf, Ψ, Γ₃_glue, CTf_2D[1], degree)


#--------------------------------------------------
# Definim els FE del model volumetric i dels multiplicadors de lagrange (un per cada columna delements) 
# Definim el model 3D
reffe_u  = ReferenceFE(lagrangian,Float64,order)

Vu = TestFESpace(Γf ,reffe_u;conformity=:H1)
Uu = TrialFESpace(Vu)

# Definim l'espai de l'acoplament
dofs  = 1
reffe = ReferenceFE(lagrangian,Float64,(1,1))
Vλ = FESpace(Γc,reffe)
#reffe = ReferenceFE(lagrangian,VectorValue{dofs,Float64},0)
#Vλ = FESpace(Γc,reffe;conformity=:L2)
Uλ = TrialFESpace(Vλ)

# Ajuntem els dos espais anteriors
V = MultiFieldFESpace([Vu,Vλ])
U = MultiFieldFESpace([Uu,Uλ])


uf = get_fe_basis(Uu)
uc = get_fe_basis(Uλ)
ucf = change_domain(uc,Γf,ReferenceDomain())

pts = get_cell_points(Γf)

yf = lazy_map(Reindex(uf(pts)),f2c_cell_map)
yc = lazy_map(Reindex(ucf(pts)),f2c_cell_map)



f2c_cell_map = vcat(glue.o2n_faces_map...)

f_dofs = lazy_map(Reindex(get_cell_dof_ids(Uu)),f2c_cell_map)
c_dofs = get_cell_dof_ids(Uλ)

#--------------------------------------------------
# Model linea fine amb periodicitat als seus extrems

modelΛf = CartesianDiscreteModel((0,1),totalColumns;isperiodic=Tuple(true))
Λf = Triangulation(modelΛf)

dofs = 1
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
if active_DOF != 0; xC[active_DOF] = -1.0; end
fun_uC = FEFunction(UΛc,xC)

# Interpolo la funcio del model lineal coarse al model lineal fine
uC_intrp = Interpolable(fun_uC)
# Portem la funcio de la linea coarse al model fine
fun_ue = interpolate(uC_intrp,UΛf)
fun_ue(Point(0.0))

ue_c = π_Λe_Γc(fun_ue,Γc)



struct InterfaceProblem
  Ω::Triangulation
  Γf::Triangulation
  Γc::Triangulation
  Λ::Triangulation
  intrf :: Vector{Intrf_Kinematic3DV2}
end

prob = InterfaceProblem(Ω, Γf, Γc, Λf, [intrf₀, intrf₁, intrf₂, intrf₃])

function func(problem::InterfaceProblem,u,v,λ,μ)

  λ_f = change_domain(λ,problem.Γf,ReferenceDomain())
  μ_f = change_domain(μ,problem.Γf,ReferenceDomain())

  contr = DomainContribution()
  for intrf in problem.intrf
    dΓfi = intrf.dΓi
    ci = ∫( (λ_f⋅v) + (μ_f⋅u) )*dΓfi
    #println(map(c->count(abs.(c[2,1]) .> tol),get_array(ci)))
    contr += ci
  end
  return contr
end

#aΩ((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ
aΩ((u,λ),(v,μ)) = ∫( ∇(v)⊙∇(u) )*dΩ
aΓ((u,λ),(v,μ)) = func(prob,u,v,λ,μ)

a((u,λ),(v,μ)) =  aΩ((u,λ),(v,μ)) + aΓ((u,λ),(v,μ))
#a((u,λ),(v,μ)) = aΩ((u,λ),(v,μ))

######fun_mult(u,v,λ,μ) = func(prob,u,v,λ,μ)
######UU = get_trial_fe_basis(U)
######VV = get_fe_basis(V)
######contrA = fun_mult(UU[1],VV[1],UU[2],VV[2])
######elementA = first(contrA.dict).second[1]

#intrf = intrf₀
#intrf = intrf₃
#Γi = intrf.Γi.trian
#dΓfi = Measure(Γi,4)
#u,λ = get_trial_fe_basis(U);
#v,μ = get_fe_basis(U);
#
#pts = get_cell_points(Γi)
#λf = change_domain(λ,Γf,ReferenceDomain())
#μf = change_domain(μ,Γf,ReferenceDomain())
#
#arr = get_array(∫( (λf⋅v) + (μf⋅u) )*dΓfi)
#
##for (i,intrf) in enumerate([intrf₀,intrf₃])
#  intrf = intrf₃
#  println(i)
#  Γi = intrf.Γi
#  pts = get_cell_points(Γi)
#  ux = u(pts)
#  λx = λf(pts)


#end



#ui = change_domain(u,prob.intrf[1].Γi,ReferenceDomain())

A = assemble_matrix(a,U,V)

#for j in [1090,1093,1096,1099]
#  i = A.rowval[A.colptr[j]:A.colptr[j+1]-1]
#  vals = A.nzval[A.colptr[j]:A.colptr[j+1]-1]
#  println(vals)
#end
#
#
#tol = 1.e-6
#for i in 1:1149
#  if abs(A[3970,i]) > tol; 
#    println(i,"  ->  ",A[3970,i])
#  end
#end
#for i in 1:1149
#  if abs(A[3973,i]) > tol; 
#    println(i,"  ->  ",A[3973,i])
#  end
#end
#for i in 1:1149
#  if abs(A[3976,i]) > tol
#    println(i,"  ->  ",A[3976,i])
#  end
#end

###for i in 1:1149
###  if abs(A[1099,i]) > tol; println(A[1099,i]); end
###end
###
###
###_A = Matrix(A)
###j = 1090
###while j <= 1149
###  idx = findall(x->abs(x)>tol,_A[1:1149,j])
###  println(j,",",length(idx))#," -> ",A[j,idx])
###  j = j+3
###end
###
###
###dof_ids = get_cell_dof_ids(Uu,Γ₀_glue)
###cells = [[1],[1,2],[2,3],[3,4]]
###for (cell,j) in enumerate([1090,1093,1096,1099])
###  idx = findall(x->x > tol , A[j,1:1149])
###  dof_ids_i = vcat(view(dof_ids,cells[cell])...)
###  dof_pos = map(dof->findlast(x->x==dof,dof_ids_i),idx)
###  dof_pos = map(dof->(dof-1)%81+1,dof_pos)
###
###  println(j," -> ",dof_pos)
###end
###
###A[1090,1]   = A[1090,1]*2
###A[1090,19]  = 0.0
###A[1090,109] = A[1090,109]*2
###A[1090,127] = 0.0
###A[1090,217] = A[1090,217]
###A[1090,223] = A[1090,223]
###A[1090,229] = A[1090,229]*4
###A[1090,235] = 0.0
###A[1090,241] = A[1090,241]/2
###A[1090,247] = A[1090,247]/2
###A[1090,691] = A[1090,691]
###A[1090,697] = A[1090,697]

#for i in 82:101
#  println( maximum(A[i,:]) )
#end

f = VectorValue(0.0)
dΓf = Measure(Γf, degree)
l((v,μ)) = ∫(v⋅f)*dΩ + ∫(μ⋅ue_c)*dΓf

b = assemble_vector(l,V)

b[(Uu.nfree+1):end]

b .= 0.0
b[9] = 1.0


#--------------------------------------------------
# Multiplicadors

#op = AffineFEOperator(a,l,U,V)
op = AffineFEOperator(U,V,A,b)

ls = LUSolver()
solver = LinearFESolver(ls)

xh = solve(op);
uh, λh = xh;

#for i in 1:3:(3*partition[1]+1)
#  println(get_free_dof_values(λh)[i])
#end

#x = get_free_dof_values(xh)
#xλ = Gridap.MultiField.restrict_to_field(U,x,2)

#λ₀ = xh.single_fe_functions[2]
#λ₁ = xh.single_fe_functions[3]
#λ₂ = xh.single_fe_functions[4]
#λ₃ = xh.single_fe_functions[5]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh])
