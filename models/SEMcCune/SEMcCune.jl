
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

prblName = "testGluedTriang"
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
add_tag_from_tags!(labels,"wall",[1,3,5,7,13,15,17,19,25])
add_tag_from_tags!(labels,"laterals",[9,11,23])
add_tag_from_tags!(labels,"faceY0",[23])
add_tag_from_tags!(labels,"faceY1",[24])
add_tag_from_tags!(labels,"faceX0",[25])
add_tag_from_tags!(labels,"faceX1",[26])
add_tag_from_tags!(labels,"line",[14])
add_tag_from_tags!(labels,"node",[3])


#--------------------------------------------------


modlType = Solid()

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------

Ω  = Triangulation(model)

#boundary_tags = ["faceX1"]
#Γ = BoundaryTriangulation(model,tags=boundary_tags)

#ΓY0  = BoundaryTriangulation(model,tags=["tag_23"]) # x-z, y = 0
#ΓY1  = BoundaryTriangulation(model,tags=["tag_24"]) # x-z, y = 1

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

axis_id = 1; face_B2_pos = maximum(lazy_map(c->c[axis_id],int_coords))
Γ_X1  = BoundaryTriangulation(model,tags=["tag_26"]) # y-z, x = 1
intrf_X1 = McCune(Γ_X1, int_coords, axis_id, face_B2_pos)
c2f_faces, cface_model_X1, Γc_X1, Γf_X1 = define_corse_fine_Γ(Γ_X1, int_coords, axis_id, face_B2_pos)

axis_id = 1; face_B2_pos = minimum(lazy_map(c->c[axis_id],int_coords))
Γ_X0  = BoundaryTriangulation(model,tags=["tag_25"]) # y-z, x = 0
intrf_X0 = McCune(Γ_X0, int_coords, axis_id, face_B2_pos)
c2f_faces, cface_model_X0, Γc_X0, Γf_X0 = define_corse_fine_Γ(Γ_X0, int_coords, axis_id, face_B1_pos)

############################################################################################

line_model_X1, Λe_X1 = get_line_model(intrf_X1)

line_model_X0, Λe_X0 = get_line_model(intrf_X0)

############################################################################################
# FESpaces 
# Model 3D
reffe_u = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
# Cares reduides
reffe_λ_X1 = ReferenceFE(lagrangian,Float64,0)
reffe_λ_X0 = ReferenceFE(lagrangian,Float64,0)

Vu = TestFESpace(Ω,reffe_u;conformity=:H1)
Vλ_X1 = FESpace(Γc_X1,reffe_λ_X1,conformity=:L2)
Vλ_X0 = FESpace(Γc_X0,reffe_λ_X0,conformity=:L2)

Uu = TrialFESpace(Vu)
Uλ_X1 = TrialFESpace(Vλ_X1)
Uλ_X0 = TrialFESpace(Vλ_X0)

V = MultiFieldFESpace([Vu,Vλ_X1,Vλ_X0])
U = MultiFieldFESpace([Uu,Uλ_X1,Uλ_X0])

# -----------------------------------------------
# Models linea
reffe_e_X1 = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
reffe_e_X0 = ReferenceFE(lagrangian,VectorValue{3,Float64},order)

Ve_X1 = FESpace(Λe_X1,reffe_e_X1,conformity=:H1)
Ve_X0 = FESpace(Λe_X0,reffe_e_X0,conformity=:H1)

Ue_X1 = TrialFESpace(Ve_X1)
Ue_X0 = TrialFESpace(Ve_X0)


# Dominis

dΩ = Measure(Ω,  degree)
dΓ_X1 = Measure(Γf_X1, degree)
dΓ_X0 = Measure(Γf_X0, degree)

n_Γ_X1 = get_normal_vector(Γf_X1)
n_Γ_X0 = get_normal_vector(Γf_X0)

#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[CT₁],
                #                            [CT₂],
                #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


tr_Γf_X1(λ) = change_domain(λ,Γf_X1,DomainStyle(λ))
tr_Γf_X0(λ) = change_domain(λ,Γf_X0,DomainStyle(λ))

_get_y(x) = VectorValue(x[2])
function π_Λe_Γc(f::CellField, Γc::Triangulation)
    _data = CellData.get_data(f)
    _cellmap = Fill(Broadcasting(_get_y),length(_data))
    data = lazy_map(∘,_data,_cellmap)
    return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
end


f = VectorValue(0.0,0.0,0.0)

xe_X1 = zero_free_values(Ue_X1); xe_X1[3] = 1.0
ue_X1 = FEFunction(Ue_X1,xe_X1)
ue_c_X1 = π_Λe_Γc(ue_X1, Γc_X1)

xe_X0 = zero_free_values(Ue_X0); xe_X0[3] = 0.0
ue_X0 = FEFunction(Ue_X0,xe_X0)
ue_c_X0 = π_Λe_Γc(ue_X0, Γc_X0)

z_coord(x) = x[3]
z_cf = CellField(z_coord,Ω)

aΩ((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

## Ordre 0
aΓ((u,λ),(v,μ)) = ∫( tr_Γf(λ)*(v⋅n_Γ) )*dΓ + ∫( tr_Γf(μ)*(u⋅n_Γ) )*dΓ

a((u,λ),(v,μ))  = aΩ((u,λ),(v,μ)) + aΓ((u,λ),(v,μ))

l((v,μ)) = ∫(v⋅f)*dΩ + 
           ∫(tr_Γf(μ*ue_c))*dΓ

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