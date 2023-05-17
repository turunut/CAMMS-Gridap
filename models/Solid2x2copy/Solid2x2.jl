#=
Exemple LINEAR ELASTIC generant malla / Dirichelt

Exemple Linear elastic en que generem un model cuadrat de forma
procedural, el mallem i definim condicions de Dirichlet imposant
desplacaments als seus extrems esquerra i dret.
=#
push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/Materials")

using Gridap
using GridapGmsh
using Gridap.Geometry
using modCT
using modModel
using modSubroutines
using Gridap.TensorValues
using Gridap.Arrays

prblName = "Solid2x2"
projFldr = pwd()

domain = (0,10,0,10,0,10)
partition = (2,2,3)
model = CartesianDiscreteModel(domain,partition)

domain_line = (0,10)
partition_line = (2)
model_line = CartesianDiscreteModel(domain_line,partition_line)

order  = 1
degree = 2*order

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"wall",[1,3,5,7,13,15,17,19,25])
add_tag_from_tags!(labels,"face",[26])
add_tag_from_tags!(labels,"line",[14])
add_tag_from_tags!(labels,"node",[3])

labels_line = get_face_labeling(model_line)


#--------------------------------------------------


modlType = Solid()

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------

topo = get_grid_topology(model)
f2n_map = Geometry.get_faces(topo,2,0)
coords = Geometry.get_vertex_coordinates(topo)

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

p = sortperm(interface_coords,by=x->x[2])

data = lazy_map(Reindex(interface_faces),p)
tabl = Table(data,y_ptrs)

Coarse_to_fine = tabl

boundary_tags = ["face"]
Γ  = BoundaryTriangulation(model,tags=boundary_tags)

o2n_cells = zeros(Int, length(Γ.glue.face_to_bgface))
child_ids = zeros(Int, length(Γ.glue.face_to_bgface))
rrules = zeros(Int, length(tabl))

for (islide, slide_list) in enumerate(tabl)
    for (izpos, num) in enumerate(slide_list)
        ind = findall(x->x==num, Γ.glue.face_to_bgface)[1]
        o2n_cells[ind] = islide
        child_ids[ind] = izpos
    end  
end


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


topo_line = get_grid_topology(model_line)
f2n_map_line = Geometry.get_faces(topo_line,1,0)
coords_line = Geometry.get_vertex_coordinates(topo_line)

face_coords_line = map(Reindex(coords_line),f2n_map_line)

com_coords_line = lazy_map(N->sum(N)/length(N),face_coords_line)

x_coords_line = map(c->c[1],com_coords_line)

rrules = sortperm(x_coords_line)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


boundary_tags = ["face"]
Γ = BoundaryTriangulation(model,tags=boundary_tags)
dΓ = Measure(Γ,degree)
n_Γ = get_normal_vector(Γ)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



#id_face = findall(x->x=="face", model.face_labeling.tag_to_name)
#id_line = findall(x->x=="line", model.face_labeling.tag_to_name)
#id_node = findall(x->x=="node", model.face_labeling.tag_to_name)
#
## Id del node de inicial
#node_init = model.face_labeling.tag_to_entities[id_node][1]
#
## Busquem a quines linies forma part aquest node
#model.grid_topology.n_m_to_nface_to_mfaces[1,2][node_init]
#
## De les anteriors agafem la linea dintre de id_line que sera linici del contorn inferior
#
#for elem in A
#  1+1
#end


#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
VS = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["wall"],
                 dirichlet_masks=[(true,true,true)])

reffe_line = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
VL = TestFESpace(model_line,reffe_line;
                 conformity=:H1,
                 dirichlet_tags=["tag_1"],
                 dirichlet_masks=[(true,true,true)])

VcouplS = ConstantFESpace(model)
VcouplB = ConstantFESpace(model_line)

V = MultiFieldFESpace([VS,VL,VcouplS,VcouplB])

g1(x) = VectorValue(0.0,0.0,0.0)
g2(x) = VectorValue(1.0,0.0,0.0)

US = TrialFESpace(VS,[g1])
UL = TrialFESpace(VL,[g2])
UcouplS = TrialFESpace(VcouplS)
UcouplB = TrialFESpace(VcouplB)
U = MultiFieldFESpace([US,UL,UcouplS,UcouplB])


#--------------------------------------------------


# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)

# Definim l'integration mesh
L = Triangulation(model_line)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dL = Measure(L,degree)


#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[CT₁],
                #                            [CT₂],
                #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


f = VectorValue(0.0,0.0)
g = 0.0

a((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ + ∫( λ*(v⋅n_Γ) )*dΓ + ∫( λ*(v⋅n_Γ) )*dL +
                 ∫( μ*(u⋅n_Γ) )*dΓ + ∫( μ*(u⋅n_Γ) )*dΓ

l((v,μ)) = ∫(v⋅f)*dΩ +
           ∫(μ*g)*dΓ


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V0)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])
