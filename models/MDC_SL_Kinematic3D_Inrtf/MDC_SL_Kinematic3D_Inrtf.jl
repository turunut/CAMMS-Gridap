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

prblName = "MDC_SL_Kinematic3D_Inrtf"
projFldr = pwd()

order = 1
degree = 2*order

############################################################################################
# Fine model
domain = (0,4,0,4,-0.5,0.5)
partition = (4,4,2)
model = CartesianDiscreteModel(domain,partition)

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"wall"    ,[1,3,5,7,13,15,17,19,25])
add_tag_from_tags!(labels,"laterals",[9,11,23])
add_tag_from_tags!(labels,"faceY0"  ,[23])
add_tag_from_tags!(labels,"faceY1"  ,[24])
add_tag_from_tags!(labels,"faceX0"  ,[25])
add_tag_from_tags!(labels,"faceX1"  ,[26])
add_tag_from_tags!(labels,"line"    ,[14])
add_tag_from_tags!(labels,"node"    ,[3])


#--------------------------------------------------


modlType = Solid()

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------

Ω  = Triangulation(model)
dΩ = Measure(Ω,degree)

Γ₁ = BoundaryTriangulation(model,tags=["tag_25"]) # y-z, x = 0
Γ₂ = BoundaryTriangulation(model,tags=["tag_26"]) # y-z, x = 1

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


axis_id = 1; face_B1_pos = minimum(lazy_map(c->c[axis_id],int_coords))
intrf₁  = Intrf_Kinematic3D(Γ₁, int_coords, axis_id, face_B1_pos, degree)

axis_id = 1; face_B2_pos = maximum(lazy_map(c->c[axis_id],int_coords))
intrf₂  = Intrf_Kinematic3D(Γ₂, int_coords, axis_id, face_B2_pos, degree)


############################################################################################

model_line₁, Λe₁ = get_line_model(intrf₁)
model_line₂, Λe₂ = get_line_model(intrf₂)

############################################################################################
# FESpaces 
# Model 3D
reffe_u  = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
# Cares reduides
reffe_λ₁ = ReferenceFE(lagrangian,VectorValue{3,Float64},0)
reffe_λ₂ = ReferenceFE(lagrangian,VectorValue{3,Float64},0)
# Model lneal
reffe_e₁ = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
reffe_e₂ = ReferenceFE(lagrangian,VectorValue{3,Float64},order)


Vu = TestFESpace(model,reffe_u;conformity=:H1)
Uu = TrialFESpace(Vu)

#Vλ_X0, Uλ_X0 = get_test_trial_spaces(intrf_X0, model)
Vλ₁ = FESpace(intrf₁.Γc,reffe_λ₁,conformity=:L2)
Vλ₂ = FESpace(intrf₂.Γc,reffe_λ₂,conformity=:L2)

Uλ₁ = TrialFESpace(Vλ₁)
Uλ₂ = TrialFESpace(Vλ₂)

V = MultiFieldFESpace([Vu,Vλ₁,Vλ₂])
U = MultiFieldFESpace([Uu,Uλ₁,Uλ₂])

## Models linea
#
Ve₁ = FESpace(Λe₁,reffe_e₁,conformity=:H1)
Ve₂ = FESpace(Λe₂,reffe_e₂,conformity=:H1)
Ue₁ = TrialFESpace(Ve₁)
Ue₂ = TrialFESpace(Ve₂)

# Dominis

#dΓ_X0 = Measure(Γf_X0, degree)
#dΓ_X1 = Measure(Γf_X1, degree)
#
#n_Γ_X0 = get_normal_vector(Γf_X0)
#n_Γ_X1 = get_normal_vector(Γf_X1)

#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[CT₁],
                #                            [CT₂],
                #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


tr_Γf₁(λ) = change_domain(λ,intrf₁.Γf,DomainStyle(λ))
tr_Γf₂(λ) = change_domain(λ,intrf₂.Γf,DomainStyle(λ))

_get_y(x) = VectorValue(x[2])
function π_Λe_Γc(f::CellField, Γc::Triangulation)
    _data = CellData.get_data(f)
    _cellmap = Fill(Broadcasting(_get_y),length(_data))
    data = lazy_map(∘,_data,_cellmap)
    return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
end


f = VectorValue(0.0,0.0,0.0)

xe₁ = zero_free_values(Ue₁); xe₁[3] = 1.0
ue₁ = FEFunction(Ue₁,xe₁)
ue_c₁ = π_Λe_Γc(ue₁, intrf₁.Γc)

xe₂ = zero_free_values(Ue₂)
ue₂ = FEFunction(Ue₂,xe₂)
ue_c₂ = π_Λe_Γc(ue₂, intrf₂.Γc)

z_coord(x) = x[3]
z_cf = CellField(z_coord,Ω)

aΩ((u,λ,α),(v,μ,β)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

## Ordre 0
aΓ₁((u,λ,α),(v,μ,β)) = ∫( tr_Γf₁(λ)⋅v )*intrf₁.dΓ + ∫( tr_Γf₁(μ)⋅u )*intrf₁.dΓ
aΓ₂((u,λ,α),(v,μ,β)) = ∫( tr_Γf₂(α)⋅v )*intrf₂.dΓ + ∫( tr_Γf₂(β)⋅u )*intrf₂.dΓ

a((u,λ,α),(v,μ,β))  =  aΩ((u,λ,α),(v,μ,β)) + 
                      aΓ₁((u,λ,α),(v,μ,β)) + 
                      aΓ₂((u,λ,α),(v,μ,β))

l((v,μ,β)) = ∫(v⋅f)*dΩ + 
             ∫(tr_Γf₁(μ)⋅ue_c₁)*intrf₁.dΓ + 
             ∫(tr_Γf₂(β)⋅ue_c₂)*intrf₂.dΓ

A = assemble_matrix(a,U,V)
b = assemble_vector(l,V)

##--------------------------------------------------

op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

xh = solve(op);
uh, λh = xh;

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])