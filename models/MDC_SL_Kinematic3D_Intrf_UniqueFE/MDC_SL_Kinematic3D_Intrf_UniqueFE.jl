
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

using FillArrays

prblName = "MDC_SL_Kinematic3D_Intrf_UniqueFE"
projFldr = pwd()

order = 1
degree = 2*order

############################################################################################
# Fine model
domain = (0,4,0,4,-0.5,0.5)
partition = (16,16,4)
model = CartesianDiscreteModel(domain,partition)

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"wall",[1,3,5,7,13,15,17,19,25])
add_tag_from_tags!(labels,"laterals",[9,11,23])
add_tag_from_tags!(labels,"face",[26])
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
dΩ = Measure(Ω,degree)

#Γ  = BoundaryTriangulation(model,tags=["tag_26"]) # y-z, x = 1
boundary_tags = ["tag_26"]
Γ = BoundaryTriangulation(model,tags=boundary_tags)

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

axis_id = 1; face_B1_pos = maximum(lazy_map(c->c[axis_id],int_coords))
intrf  = Intrf_Kinematic3D(Γ, int_coords, axis_id, face_B1_pos, degree)

############################################################################################
# FESpaces 

reffe_u = ReferenceFE(lagrangian,VectorValue{3,Float64},order)

Vu = TestFESpace(Ω,reffe_u;
                 conformity=:H1,
                 dirichlet_tags=["wall"],
                 dirichlet_masks=[(true,true,true)])
g1(x) = VectorValue(0.0,0.0,0.0)
Uu = TrialFESpace(Vu,[g1])

Vλ, Uλ = get_test_trial_spaces(intrf)

U = MultiFieldFESpace([Uu,Uλ])
V = MultiFieldFESpace([Vu,Vλ])

Ve, Ue = get_line_test_trial_spaces(intrf,order)


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

xe = zero_free_values(Ue); xe[16] = 1.0
ue = FEFunction(Ue,xe)
ue_c = π_Λe_Γc(ue,intrf.Γc)

z_coord(x) = x[3]
z_cf = CellField(z_coord,Ω)


aΩ((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

aΓ((u,λ),(v,μ)) = contribute_matrix(intrf, (u,λ),(v,μ), 1, 2)

a((u,λ),(v,μ))  = aΩ((u,λ),(v,μ)) + aΓ((u,λ),(v,μ))


la((v,μ)) = contribute_vector(intrf, (v,μ), 2, ue_c)

l((v,μ)) = ∫(v⋅f)*dΩ + la((v,μ))


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