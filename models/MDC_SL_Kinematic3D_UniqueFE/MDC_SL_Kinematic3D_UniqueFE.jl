
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

using FillArrays

prblName = "MDC_SL_Kinematic3D_UniqueFE"
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

#Γ  = BoundaryTriangulation(model,tags=["tag_26"]) # y-z, x = 1
boundary_tags = ["face"]
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

face_x = maximum(lazy_map(c->c[1],int_coords))
interface_faces = findall(c->c[1]==face_x,int_coords)
interface_coords = view(int_coords,interface_faces)

y_coords = map(c->c[2],interface_coords)
y_unique = sort(unique(map(c->c[2],interface_coords)))
y_counts = [count(==(y),y_coords) for y in sort(unique(y_coords))]
y_ptrs = Gridap.Adaptivity.counts_to_ptrs(y_counts)

perm = sortperm(interface_coords,by=x->x[2])
data = lazy_map(Reindex(interface_faces),perm)
c2f_faces = Table(data,y_ptrs)

n2o_cells = zeros(Int, length(Γ.glue.face_to_bgface))
child_ids = zeros(Int, length(Γ.glue.face_to_bgface))
for (islide, slide_list) in enumerate(c2f_faces)
    for (izpos, num) in enumerate(slide_list)
        ind = findall(x->x==num, Γ.glue.face_to_bgface)[1]
        n2o_cells[ind] = islide
        child_ids[ind] = izpos
    end  
end

n2o_faces = [Int[],Int[],n2o_cells]
rrules = Fill(RefinementRule(QUAD,(1,length(c2f_faces[1]))),length(c2f_faces))
glue = AdaptivityGlue(n2o_faces,child_ids,rrules) # From coarse to fine

cface_model = CartesianDiscreteModel((0,1,0,1),(length(c2f_faces),1))
Γc = Triangulation(cface_model)
Γf = Adaptivity.GluedTriangulation(Γ,Γc,glue)

############################################################################################

line_model = CartesianDiscreteModel((0,1),(length(c2f_faces)))
Λe = Triangulation(line_model)

############################################################################################
# FESpaces 

reffe_u = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
reffe_λ = ReferenceFE(lagrangian,VectorValue{3,Float64},0)
reffe_e = ReferenceFE(lagrangian,VectorValue{3,Float64},order)

Vu = TestFESpace(Ω,reffe_u;
                 conformity=:H1,
                 dirichlet_tags=["wall"],
                 dirichlet_masks=[(true,true,true)])
g1(x) = VectorValue(0.0,0.0,0.0)
Uu = TrialFESpace(Vu,[g1])

Vλ = FESpace(Γc,reffe_λ,conformity=:L2)
Uλ = TrialFESpace(Vλ)

U = MultiFieldFESpace([Uu,Uλ])
V = MultiFieldFESpace([Vu,Vλ])

Ve = FESpace(Λe,reffe_e,conformity=:H1)
Ue = TrialFESpace(Ve)

dΩ = Measure(Ω,  degree)
dΓ = Measure(Γf, degree)

n_Γ = get_normal_vector(Γf)

#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[CT₁],
                #                            [CT₂],
                #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


tr_Γf(λ) = change_domain(λ,Γf,DomainStyle(λ))

_get_y(x) = VectorValue(x[2])
function π_Λe_Γc(f::CellField)
    _data = CellData.get_data(f)
    _cellmap = Fill(Broadcasting(_get_y),length(_data))
    data = lazy_map(∘,_data,_cellmap)
    return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
end


f = VectorValue(0.0,0.0,0.0)

xe = zero_free_values(Ue); xe[16] = 1.0
ue = FEFunction(Ue,xe)
ue_c = π_Λe_Γc(ue)

z_coord(x) = x[3]
z_cf = CellField(z_coord,Ω)

aΩ((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

## Ordre 0
aΓ((u,λ),(v,μ)) = ∫( tr_Γf(λ)⋅v )*dΓ + ∫( tr_Γf(μ)⋅u )*dΓ
# Ordre 1
#aΓ((u,λ),(v,μ)) = ∫( tr_Γf(λ)*(ctˣ(z_cf,1)*v⋅n_Γ) )*dΓ + ∫( tr_Γf(μ)*(ctˣ(z_cf,1)*u⋅n_Γ) )*dΓ

a((u,λ),(v,μ))  = aΩ((u,λ),(v,μ)) + aΓ((u,λ),(v,μ))

l((v,μ)) = ∫(v⋅f)*dΩ + 
           ∫(tr_Γf(μ)⋅ue_c)*dΓ

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