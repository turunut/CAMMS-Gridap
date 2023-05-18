
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
Γ_X1  = BoundaryTriangulation(model,tags=["tag_26"]) # y-z, x = 1
#ΓY1  = BoundaryTriangulation(model,tags=["tag_24"]) # x-z, y = 1
Γ_X0  = BoundaryTriangulation(model,tags=["tag_25"]) # y-z, x = 0

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

function create_Interface(Γ::Triangulation, int_coords::Vector{VectorValue{3, Int64}}, axis_id::Int64, axis_int_coord::Int64)
    axis_p = [2,1][axis_id]
    println(axis_p)
    face_x = axis_int_coord
    interface_faces = findall(c->c[axis_id]==face_x,int_coords)
    interface_coords = view(int_coords,interface_faces)
    
    y_coords = map(c->c[axis_p],interface_coords)
    y_unique = sort(unique(map(c->c[axis_p],interface_coords)))
    y_counts = [count(==(y),y_coords) for y in sort(unique(y_coords))]
    y_ptrs = Gridap.Adaptivity.counts_to_ptrs(y_counts)
    
    perm = sortperm(interface_coords,by=x->x[axis_p])
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
    return glue, c2f_faces
end

############################################################################################

axis_id = 1; face_B2_pos = maximum(lazy_map(c->c[axis_id],int_coords))
glue_X1, c2f_faces_X1 = create_Interface(Γ_X1, int_coords, axis_id, face_B2_pos)

cface_model_X1 = CartesianDiscreteModel((0,1,0,1),(length(c2f_faces_X1),1))
Γc_X1 = Triangulation(cface_model_X1)
Γf_X1 = Adaptivity.GluedTriangulation(Γ_X1,Γc_X1,glue_X1)
#-------------------------------------------------------------------------------------------
axis_id = 1; face_B2_pos = minimum(lazy_map(c->c[axis_id],int_coords))
glue_X0, c2f_faces_X0 = create_Interface(Γ_X0, int_coords, axis_id, face_B2_pos)

cface_model_X0 = CartesianDiscreteModel((0,1,0,1),(length(c2f_faces_X0),1))
Γc_X0 = Triangulation(cface_model_X0)
Γf_X0 = Adaptivity.GluedTriangulation(Γ_X0,Γc_X0,glue_X0)

############################################################################################

#line_model = CartesianDiscreteModel((0,1),(length(c2f_faces)))
#Λe = Triangulation(line_model)

line_model_X1 = CartesianDiscreteModel((0,1),(length(c2f_faces_X1)))
Λe_X1 = Triangulation(line_model_X1)

#line_model = CartesianDiscreteModel((0,1),(length(c2f_faces)))
#Λe = Triangulation(line_model)

line_model_X0 = CartesianDiscreteModel((0,1),(length(c2f_faces_X0)))
Λe_X0 = Triangulation(line_model_X0)

############################################################################################
# FESpaces 

reffe_u = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
reffe_λ_X1 = ReferenceFE(lagrangian,Float64,0)
reffe_e_X1 = ReferenceFE(lagrangian,Float64,order)
reffe_λ_X0 = ReferenceFE(lagrangian,Float64,0)
reffe_e_X0 = ReferenceFE(lagrangian,Float64,order)

Vu = TestFESpace(Ω,reffe_u;conformity=:H1)
Uu = TrialFESpace(Vu)

Vλ_X1 = FESpace(Γc_X1,reffe_λ_X1,conformity=:L2)
Uλ_X1 = TrialFESpace(Vλ_X1)

Vλ_X0 = FESpace(Γc_X0,reffe_λ_X0,conformity=:L2)
Uλ_X0 = TrialFESpace(Vλ_X0)

Ve_X1 = FESpace(Λe_X1,reffe_e_X1,conformity=:H1)
Ue_X1 = TrialFESpace(Ve_X1)

Ve_X0 = FESpace(Λe_X0,reffe_e_X0,conformity=:H1)
Ue_X0 = TrialFESpace(Ve_X0)

U = MultiFieldFESpace([Uu,Uλ_X1,Uλ_X0])
V = MultiFieldFESpace([Vu,Vλ_X1,Vλ_X0])

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
function π_Λe_Γc(f::CellField)
    _data = CellData.get_data(f)
    _cellmap = Fill(Broadcasting(_get_y),length(_data))
    data = lazy_map(∘,_data,_cellmap)
    return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
end


f = VectorValue(0.0,0.0,0.0)

xe_X1 = zero_free_values(Ue_X1); xe_X1[3] = 1.0
ue_X1 = FEFunction(Ue_X1,xe_X1)
ue_c = π_Λe_Γc(ue)

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