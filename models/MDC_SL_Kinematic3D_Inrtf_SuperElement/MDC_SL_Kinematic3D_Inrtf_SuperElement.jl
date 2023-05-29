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

prblName = "MDC_SL_Kinematic3D_Inrtf_SuperElement"
projFldr = pwd()

order = 2
degree = 2*order

############################################################################################
# Fine model
domain = (0,4,0,4,-0.5,0.5)
partition = (12,12,4)
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

Γ₁ = BoundaryTriangulation(model,tags=["tag_24"]) # y-z, x = 0
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

axis_id = 2; face_B1_pos = maximum(lazy_map(c->c[axis_id],int_coords))
intrf₁  = Intrf_Kinematic3D(Γ₁, int_coords, axis_id, face_B1_pos, degree, true)

axis_id = 1; face_B2_pos = maximum(lazy_map(c->c[axis_id],int_coords))
intrf₂  = Intrf_Kinematic3D(Γ₂, int_coords, axis_id, face_B2_pos, degree, false)

############################################################################################
# FESpaces 
# Model 3D
reffe_u  = ReferenceFE(lagrangian,VectorValue{3,Float64},order)


Vu = TestFESpace(model,reffe_u;conformity=:H1)
Uu = TrialFESpace(Vu)

#Vλ_X0, Uλ_X0 = get_test_trial_spaces(intrf_X0, model)
Vλ₁, Uλ₁ = get_test_trial_spaces(intrf₁)
Vλ₂, Uλ₂ = get_test_trial_spaces(intrf₂)

V = MultiFieldFESpace([Vu,Vλ₁,Vλ₂])
U = MultiFieldFESpace([Uu,Uλ₁,Uλ₂])

## Models linea
Ve₁, Ue₁, reffe_e₁ = get_line_test_trial_spaces(intrf₁,order)
Ve₂, Ue₂, reffe_e₂ = get_line_test_trial_spaces(intrf₂,order)

#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1)

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


_get_y(x) = VectorValue(x[2])
function π_Λe_Γc(f::CellField, Γc::Triangulation)
    _data = CellData.get_data(f)
    _cellmap = Fill(Broadcasting(_get_y),length(_data))
    data = lazy_map(∘,_data,_cellmap)
    return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
end


f = VectorValue(0.0,0.0,0.0)


domainC = (0, 1); partitionC = (3)
modelC  = CartesianDiscreteModel(domainC, partitionC)
ΩC = Triangulation(modelC)
VC = FESpace(ΩC,reffe_e₁,conformity=:H1); UC = TrialFESpace(VC)
# Definim la funcio a partir del DOFs
xC = zero_free_values(UC); xC[5] = 1.0
uC = FEFunction(UC,xC)
# -----------------
uC_intrp = Interpolable(uC)

ue₁ = interpolate(uC_intrp,Ue₁)
ue_c₁ = π_Λe_Γc(ue₁,intrf₁.Γc)

#for i in 0.0:0.025:1.0
#    println( i, "  ", ue₁( Point(i) ) )
#end
#
#for i in 0.0:0.025:1.0
#    println( i, "  ", uC( Point(i) ) )
#end

xe₂ = zero_free_values(Ue₂); xe₂[13] = 0.0
ue₂ = FEFunction(Ue₂,xe₂)
ue_c₂ = π_Λe_Γc(ue₂,intrf₂.Γc)


z_coord(x) = x[3]
z_cf = CellField(z_coord,Ω)

aΩ((u,λ,α),(v,μ,β)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ



aΓ₁((u,λ,α),(v,μ,β)) = contribute_matrix(intrf₁, (u,λ,α), (v,μ,β), 1, 2)
aΓ₂((u,λ,α),(v,μ,β)) = contribute_matrix(intrf₂, (u,λ,α), (v,μ,β), 1, 3)

la₁((v,μ,β)) = contribute_vector(intrf₁, (v,μ,β), 2, ue_c₁)
la₂((v,μ,β)) = contribute_vector(intrf₂, (v,μ,β), 3, ue_c₂)

a((u,λ,α),(v,μ,β)) =  aΩ((u,λ,α),(v,μ,β)) + 
                     aΓ₁((u,λ,α),(v,μ,β)) + 
                     aΓ₂((u,λ,α),(v,μ,β))

l((v,μ,β)) = ∫(v⋅f)*dΩ + la₁((v,μ,β)) + la₂((v,μ,β))

A = assemble_matrix(a,U,V)
b = assemble_vector(l,V)

##--------------------------------------------------

#op = AffineFEOperator(a,l,U,V)
op = AffineFEOperator(U,V,A,b)

ls = LUSolver()
solver = LinearFESolver(ls)

xh = solve(op);
uh, λh = xh;

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])