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

prblName = "MDC_SL_Reissner_Inrtf_SuperElement"
projFldr = pwd()

order = 2
degree = 2*order

############################################################################################
# Fine model
domain = (0,4,0,4,-0.5,0.5)
partition = (4,4,2)
model = CartesianDiscreteModel(domain,partition)

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"faceY0",[23])
add_tag_from_tags!(labels,"faceY1",[24])
add_tag_from_tags!(labels,"faceX0",[25])
add_tag_from_tags!(labels,"faceX1",[26])


#--------------------------------------------------


modlType = Solid()

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------


Ω  = Triangulation(model)
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1)

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#---------------------------


modlType_2D = PlaneStress()

ct1_2D = modModel.computeCT(modlType_2D, CT1, true)
ct2_2D = modModel.computeCT(modlType_2D, CT2, true)

CTs_2D = hcat(ct1_2D)

CTf_2D = get_CT_CellField(modlType_2D, CTs_2D, tags, Ω)


#--------------------------------------------------

Γ₉  = BoundaryTriangulation(model,tags=["tag_23","tag_26","tag_24","tag_25"])

Γ₀ = BoundaryTriangulation(model,tags=["tag_23"])
Γ₁ = BoundaryTriangulation(model,tags=["tag_26"])
Γ₂ = BoundaryTriangulation(model,tags=["tag_24"])
Γ₃ = BoundaryTriangulation(model,tags=["tag_25"])

Ψ = EdgeTriangulation(model,["tag_17"])


############################################################################################


# Coarse model

Dc      = num_cell_dims(model)
topo    = get_grid_topology(model)
f2n_map = Geometry.get_faces(topo,Dc-1,0)
coords  = Geometry.get_vertex_coordinates(topo)
face_coords = map(Reindex(coords),f2n_map)

tol = 1.e-3
com_coords = lazy_map(N->sum(N)/length(N),face_coords) # center of mass
int_coords = map(N->VectorValue(Int(floor(N[1]/tol)),Int(floor(N[2]/tol)),Int(floor(N[3]/tol))),com_coords)


############################################################################################


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

Γ = [Γ₀,Γ₁,Γ₂,Γ₃]
axis_id = [axis_id₀, axis_id₁, axis_id₂,axis_id₃]
axis_int_coord = [face_pos₀,face_pos₁,face_pos₂,face_pos₃]
inv_list = [false,false,true,true]


c2f_faces_list = []
n2o_faces_g = Vector{Vector{Int}}(undef, 0)
child_ids_g = Vector{Int}(undef, 0)

append!(n2o_faces_g, [[],[],[]])

global totalColumns = 0
global coun_max_n2o = 0
for iintrf in eachindex(Γ)
  n2o_faces, child_ids, c2f_faces = get_comp_glue(Γ[iintrf], int_coords, axis_id[iintrf], axis_int_coord[iintrf], inv_list[iintrf]) 
  n2o_faces[3] .+= coun_max_n2o
  append!(n2o_faces_g[3], n2o_faces[3])
  append!(child_ids_g, child_ids)
  global totalColumns += length(c2f_faces)
  push!(c2f_faces_list, c2f_faces)
  global coun_max_n2o = maximum(n2o_faces[3])
end

c2f_faces_glo = append_tables_globally(c2f_faces_list[1], c2f_faces_list[2], c2f_faces_list[3], c2f_faces_list[4])

rrules_g = Fill(RefinementRule(QUAD,(1,length(c2f_faces_glo[1]))),length(c2f_faces_glo))
glue = AdaptivityGlue(n2o_faces_g,child_ids_g,rrules_g) # From coarse to fine

cface_model = CartesianDiscreteModel((0,1,0,1),(totalColumns,1),isperiodic=(true,false))
Γc  = Triangulation(cface_model)
Γf  = Gridap.Adaptivity.GluedTriangulation(Γ₉,Γc,glue)


intrf₀  = Intrf_ReissnerV2(Ω, Γf, Ψ, Γ₀, CTf_2D[1], degree)

intrf₁  = Intrf_ReissnerV2(Ω, Γ₁, Ψ, Γ₁, CTf_2D[1], degree)

intrf₂  = Intrf_ReissnerV2(Ω, Γ₂, Ψ, Γ₂, CTf_2D[1], degree)

intrf₃  = Intrf_ReissnerV2(Ω, Γ₃, Ψ, Γ₃, CTf_2D[1], degree)


############################################################################################
# FESpaces 
# Model 3D
reffe_u  = ReferenceFE(lagrangian,VectorValue{3,Float64},order)

Vu = TestFESpace(model,reffe_u;conformity=:H1)
Uu = TrialFESpace(Vu)

dofs  = 5
reffe = ReferenceFE(lagrangian,VectorValue{dofs,Float64},0)
Vλ = FESpace(Γc,reffe,conformity=:L2)
Uλ = TrialFESpace(Vλ)

V = MultiFieldFESpace([Vu,Vλ])
U = MultiFieldFESpace([Uu,Uλ])

## Model linea
modelΛf = CartesianDiscreteModel((0,1),totalColumns;isperiodic=Tuple(true))
Λf = Triangulation(modelΛf)

dofs = 5
reffeΛ = ReferenceFE(lagrangian,VectorValue{dofs,Float64},order)
VΛf = FESpace(Λf,reffeΛ,conformity=:H1)
UΛf = TrialFESpace(VΛf)


#--------------------------------------------------


f = VectorValue(0.0,0.0,0.0)

interfaces = [intrf₀, intrf₁, intrf₂, intrf₃]
crs_dsc = [1,1,1,1]

domainΛc = (0, 1); partitionΛc = (totalColumns)
modelΛc  = CartesianDiscreteModel(domainΛc, partitionΛc)
Λc  = Triangulation(modelΛc)
VΛc = FESpace(Λc,reffeΛ,conformity=:H1); UΛc = TrialFESpace(VΛc)

active_DOF = 1

_get_y(x) = VectorValue(x[2])
function π_Λe_Γc(f::CellField, Γc::Triangulation)
  _data = CellData.get_data(f)
  _cellmap = Fill(Broadcasting(_get_y),length(_data))
  data = lazy_map(∘,_data,_cellmap)
  return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
end

xC = zero_free_values(UΛc);
if active_DOF != 0; xC[active_DOF] = 1.0; end
fun_uC = FEFunction(UΛc,xC)
# -----------------
uC_intrp = Interpolable(fun_uC)

fun_ue = interpolate(uC_intrp,UΛf)
ue_c = π_Λe_Γc(fun_ue,Γc)





struct InterfaceProblem
  Ω::Triangulation
  Γf::Triangulation
  Γc::Triangulation
  Λ::Triangulation
  intrf :: Vector{Intrf_ReissnerV2}
end

prob = InterfaceProblem(Ω, Γf, Γc, Λf, [intrf₀, intrf₁, intrf₂, intrf₃])

function func(problem::InterfaceProblem,u,v,λ,μ)

  λ_f = change_domain(λ,problem.Γf,ReferenceDomain())
  μ_f = change_domain(μ,problem.Γf,ReferenceDomain())

  contr = DomainContribution()
  for intrf in problem.intrf
    Γfi  = intrf.Γi
    dΓfi = intrf.dΓi

    invD = intrf.invD
    invA = intrf.invA
  
    da_fun(CT,zf,z_val) = sum(∫(    step_field(zf,z_val,intrf.Γf)*CT )*intrf.dΓf)
    db_fun(CT,zf,z_val) = sum(∫( zf*step_field(zf,z_val,intrf.Γf)*CT )*intrf.dΓf)
    da(z_val) = da_fun(intrf.CTf_2D,intrf.zf,z_val)
    db(z_val) = db_fun(intrf.CTf_2D,intrf.zf,z_val)
  
    function f_da(x); z_val = x[end]; return sum( ∫(          step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end
    function f_db(x); z_val = x[end]; return sum( ∫( intrf.zf*step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end
  
    function _my_tensor(z,CT_2D)
      da_db_arr = zeros(4,2)
  
      α_arr = invD ⋅ ( TensorValue{6,3}( 1.0, 0.0, 0.0,   z, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0,   z, 0.0,
                                         0.0, 0.0, 1.0, 0.0, 0.0,   z) ⋅ CT_2D )
  
      AAA = f_da(z)
      BBB = f_db(z)
      
      da_db_arr[1:2,1:2] = get_array( AAA )[1:2,1:2]
      da_db_arr[3:4,1:2] = get_array( BBB )[1:2,1:2]
      
      β_arr = invA ⋅ ( TensorValue{4,2}( da_db_arr ) )
      
      A_arr = zeros(5,3)
      A_arr[1:4,1:2] .= get_array( α_arr )[[1,3,4,6],[1,3]]
      A_arr[5,3]      = get_array( β_arr )[1,1]
      return TensorValue{5,3}(A_arr)
    end
  
    T = _my_tensor∘(intrf.zf, intrf.CTf_2D)

    pts = get_cell_points(dΓfi)
    T(pts)

    contr += ∫(λ_f⋅T⋅v + μ_f⋅T⋅u)*dΓfi
  end
  return contr
end

aΓ((u,λ),(v,μ)) = func(prob,u,v,λ,μ)
aΩ((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

a((u,λ),(v,μ)) =  aΩ((u,λ),(v,μ)) + aΓ((u,λ),(v,μ))

u,λ = get_trial_fe_basis(U)
v,μ = get_fe_basis(U)


#Γi = BoundaryTriangulation(Γf,tags=["tag_23"])

Γi = view(Γf,[1,2,3,4,5,6,7,8])

pts = get_cell_points(Γi)
λf = change_domain(λ,Γf,ReferenceDomain())
λfi = change_domain(λf,Γi,ReferenceDomain())


ui = change_domain(u,prob.intrf[1].Γi,ReferenceDomain())

A = assemble_matrix(a,U,V)



l((v,μ₀,μ₁,μ₂,μ₃)) = ∫(v⋅f)*dΩ + ∫( μ_glo⋅ue_c )*dΓf

b = assemble_vector(l,V)

##--------------------------------------------------
#Multiplicadors

#op = AffineFEOperator(a,l,U,V)
op = AffineFEOperator(U,V,A,b)

ls = LUSolver()
solver = LinearFESolver(ls)

xh = solve(op);
uh, λh = xh;

λ₀ = xh.single_fe_functions[2]
λ₁ = xh.single_fe_functions[3]
λ₂ = xh.single_fe_functions[4]
λ₃ = xh.single_fe_functions[5]


writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])