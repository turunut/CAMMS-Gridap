#=
Exemple LINEAR ELASTIC generant malla / Dirichelt

Exemple Linear elastic en que generem un model cuadrat de forma
procedural, el mallem i definim condicions de Dirichlet imposant
desplacaments als seus extrems esquerra i dret.
=#
push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/Materials")

using Gridap
#using GridapGmsh
using GridapGiD
using Gridap.Geometry
using modCT
using modModel
using modSubroutines
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.ReferenceFEs

prblName = "Timoshenko2D"
projFldr = pwd()

#model = GiDDiscreteModel( projFldr*"/models/"*prblName*"/"*prblName )

#--------------------------------------------------


#https://github.com/gridapapps/GridapGeosciences.jl/blob/master/src/CubedSphereDiscreteModels.jl

rot = deg2rad(30.0)

nodes = [VectorValue( 00.0*1.25*cos(rot), 00.0*1.25*sin(rot)),
         VectorValue( 10.0*1.25*cos(rot), 10.0*1.25*sin(rot)),
         VectorValue( 20.0*1.25*cos(rot), 20.0*1.25*sin(rot)),
         VectorValue( 30.0*1.25*cos(rot), 30.0*1.25*sin(rot)),
         VectorValue( 40.0*1.25*cos(rot), 40.0*1.25*sin(rot)),
         VectorValue( 50.0*1.25*cos(rot), 50.0*1.25*sin(rot)),
         VectorValue( 60.0*1.25*cos(rot), 60.0*1.25*sin(rot)),
         VectorValue( 70.0*1.25*cos(rot), 70.0*1.25*sin(rot)),
         VectorValue( 80.0*1.25*cos(rot), 80.0*1.25*sin(rot))]
c2n_map = Table([1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9],
                [1,   3,   5,   7,   9,   11,  13,  15,17])

cell_type = Int8[1,1,1,1,1,1,1,1]
polys = [SEGMENT]
reffes = map(p->LagrangianRefFE(Float64,p,1),polys)
orientation = NonOriented()

topo = UnstructuredGridTopology(nodes,c2n_map,cell_type,polys,orientation)
grid = UnstructuredGrid(nodes,c2n_map,reffes,cell_type,orientation)

#d_to_num_dfaces_to_entity = [[1,3,3,3,2],[3,3,3,3,3]]
#tag_to_name   = ["interior","boundary","left","right"]
#tag_to_entity = [[3],[1,2],[1],[2]]
#tag_to_name     = ["interior", "boundary", "left", "right"]
#tag_to_entities = [[3],        [1,2],      [1],   [2]  ]
#d_to_dface_to_entity = [[1,3,3,3,3,3,3,3,2],[3,3,3,3,3,3,3,3]] # nodes edges
tag_to_name     = ["interior", "left", "right"]
tag_to_entities = [[3],        [1],   [2]  ]
d_to_dface_to_entity = [[1,3,3,3,3,3,3,3,2],[3,3,3,3,3,3,3,3]] # nodes edges

#labels = FaceLabeling(topo)
labels = FaceLabeling(d_to_dface_to_entity, tag_to_entities, tag_to_name)

model = UnstructuredDiscreteModel(grid,topo,labels)

Dc = num_cell_dims(model)
Dp = num_point_dims(model)

Ω = Triangulation(model)

pts = get_cell_points(Ω)


#--------------------------------------------------


order  = 1
degree = 2*order

#writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)


#--------------------------------------------------


MT = Timoshenko(1.0, 5.0, -5.0)

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(MT, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(MT, CT2)


#--------------------------------------------------


reffe2 = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
reffe1 = ReferenceFE(lagrangian,VectorValue{1,Float64},order)
Vv = TestFESpace(model,reffe2;
                 conformity=:H1,
                 dirichlet_tags=["left", "right"],
                 dirichlet_masks=[(true,true), (true,true)])
Vt = TestFESpace(model,reffe1;
                 conformity=:H1,
                 dirichlet_tags=["left", "right"],
                 dirichlet_masks=[(true), (false)])
V = MultiFieldFESpace([Vv,Vt])

g0_1(x) = VectorValue(  0.0)
g1_1(x) = VectorValue(-10.0)
g0_2(x) = VectorValue(0.0, 0.0)
g1_2(x) = VectorValue(0.0,-1.0)

Uv = TrialFESpace(Vv, [g0_2,g1_2])
Ut = TrialFESpace(Vt, [g0_1,g0_1])
U = MultiFieldFESpace([Uv,Ut])


#--------------------------------------------------


# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 1 # linies
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1)

CTf = get_CT_CellField(MT, CTs, tags, Ω)


#--------------------------------------------------


#function my_rotations(m)
#  t = m(VectorValue(1.0)) - m(VectorValue(0.0)) # Posicio local 0.0 i posicio local 1.0
#  t = t/norm(t)
#  L = TensorValue([ t[1] t[2];
#                   -t[2] t[1] ]) # lᵀ=l⁻¹
#  return L
#end
#getL_op(tf) = TensorValue([ tf[1] -tf[2];
#                            tf[2]  tf[1] ])
#getL(tf::CellField) = Operation(getL_op)(tf)
#function g2l(tf, u)
#    eval = getL(tf)⋅u
#    return getL(tf)⋅u # ∘(getL(tf)⋅u)
#end

function my_tangent(m)
  t = m(VectorValue(1.0)) - m(VectorValue(0.0)) # Posicio local 0.0 i posicio local 1.0
  return t/norm(t)
end
function my_normal(m)
  t = my_tangent(m)
  return TensorValue([0.0 -1.0; 1.0 0.0])⋅t
end

function get_tangent_vector(Ω::Triangulation{1})
    cmaps = get_cell_map(Ω)
    return CellField(lazy_map(my_tangent,cmaps),Ω)
end
function get_normal_vector(Ω::Triangulation{1})
    cmaps = get_cell_map(Ω)
    return CellField(lazy_map(my_normal,cmaps),Ω)
end

tf = get_tangent_vector(Ω)
nf = get_normal_vector(Ω)


#get₁(x) = VectorValue(VectorValue(1.0,0.0)⋅x)
#get₂(x) = VectorValue(VectorValue(0.0,1.0)⋅x)
getₑ(x,ê) = VectorValue( ê ⋅ x )
#getₒ(x) = VectorValue(VectorValue(1.0,1.0)⋅x)
#n0(x) = sign(sign(sign(x)+1)-0.5) # Torna 1 si el x igual o me es gran que 0, else retorna -1
#n0₂(a,b) = sign(a+b)
#signe(x,y) = (x*y)/norm(x*y)      # Torna 1 si x y tenen el mateix signe sino -1
#getₙ₁(x) = VectorValue( sign(sum(x))*norm(x) )
#getₙ₂(x) = VectorValue( sign(sum(x))*norm(x) )

getₙ(x) = VectorValue( sign(sum(x))*norm(x) )
#∂ₙ(u,ê,d) = getₑ ∘ ( ∇(u)⋅ê, d )
∂ₙ(u,ê,d) = getₑ ∘ ( ∇(u)⋅ê, d )
∂ᵥ(θ) = getₙ ∘ ∇(θ)

a((u,θ),(v,t)) = ∫( ∂ₙ(v,tf,tf)⊙σₑ(CTf[1],∂ₙ(u,tf,tf)) + ∂ᵥ(t)⊙σₑ(CTf[2],∂ᵥ(θ)) )*dΩ + # Axial         + Axial/Bending
                 ∫( ∂ₙ(v,tf,tf)⊙σₑ(CTf[2],∂ₙ(u,tf,tf)) + ∂ᵥ(t)⊙σₑ(CTf[3],∂ᵥ(θ)) )*dΩ + # Bending/Axial + Bending
                 ∫( γ(MT,∂ₙ(v,nf,tf),t) ⊙ σₑ(CTf[4], γ(MT,∂ₙ(u,nf,tf),θ)) )*dΩ 

l((v,t)) = 0


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)
vh = uh.single_fe_functions[1]
th = uh.single_fe_functions[2]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u" => vh,
                     "θ" => th,
                     "ε" => ∂ₙ(vh,tf,tf),
                     "κ" => ∂ᵥ(th),
                     "γ" => ∂ₙ(vh,nf,tf)])


#----------------------------------


A((u,θ),(v,t)) = ∫( ∂ₙ(v,tf,tf)⊙σₑ(CTf[1],∂ₙ(u,tf,tf)) + ∂ᵥ(t)⊙σₑ(CTf[2],∂ᵥ(θ)) )*dΩ
D((u,θ),(v,t)) = ∫( ∂ₙ(v,tf,tf)⊙σₑ(CTf[2],∂ₙ(u,tf,tf)) + ∂ᵥ(t)⊙σₑ(CTf[3],∂ᵥ(θ)) )*dΩ
S((u,θ),(v,t)) = ∫( γ(MT,∂ₙ(v,nf,tf),t) ⊙ σₑ(CTf[4], γ(MT,∂ₙ(u,nf,tf),θ)) )*dΩ

UU = get_trial_fe_basis(U)
VV = get_fe_basis(V)
contrA = A(UU,VV)
elementA = first(contrA.dict).second[1]
contrD = D(UU,VV)
elementD = first(contrD.dict).second[1]
contrS = S(UU,VV)
elementS = first(contrS.dict).second[1]

contr = a(UU,VV)
element = first(contr.dict).second[1]

