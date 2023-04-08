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

rot = deg2rad(0.0)

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


order = 1

#writevtk(model,"model")

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
g1_2(x) = VectorValue(1.0,-1.0)

Uv = TrialFESpace(Vv, [g0_2,g1_2])
Ut = TrialFESpace(Vt, [g0_1,g0_1])
U = MultiFieldFESpace([Uv,Ut])


#--------------------------------------------------


degree = 1*order
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


function my_rotations(m)
  t = m(VectorValue(1.0)) - m(VectorValue(0.0)) # Posicio local 0.0 i posicio local 1.0
  t = t/norm(t)
  L = TensorValue([ t[1] t[2];
                   -t[2] t[1] ]) # lᵀ=l⁻¹
  return L
end

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

#getL_op(tf) = TensorValue([ tf[1] -tf[2];
#                            tf[2]  tf[1] ])
#getL(tf::CellField) = Operation(getL_op)(tf)
#function g2l(tf, u)
#    eval = getL(tf)⋅u
#    return getL(tf)⋅u # ∘(getL(tf)⋅u)
#end

get₁(x) = VectorValue(VectorValue(1.0,0.0)⋅x)
get₂(x) = VectorValue(VectorValue(0.0,1.0)⋅x)
getₒ(x) = VectorValue(VectorValue(1.0,1.0)⋅x)
getₙ(x) = VectorValue(sign(x[1])*norm(x))

#∂₁(u,êf) = getₒ ∘ (∇(u) ⋅ êf)
∂₁(u,êf) = getₙ ∘ (∇(u) ⋅ êf)
#∂₂(u,êf) = getₒ ∘ (∇(u) ⋅ êf)
∂ₒ(θ) = VectorValue(1.0,1.0) ⋅ (∇(θ))
∂ᵥ(θ,êf) = êf ⋅ ∇(θ)

op_TB(êf, θ) = VectorValue(êf[1]*θ[1], êf[2]*θ[1])
∂ₗ(θ,êf) = getₙ ∘ Operation(op_TB)(êf, θ)

a((u,θ),(v,t)) = ∫( ∂₁(v,tf)⊙σₑ(CTf[1],∂₁(u,tf)) + ∂ᵥ(t,tf)⊙σₑ(CTf[2],∂ᵥ(θ,tf)) )*dΩ + # Axial         + Axial/Bending
                 ∫( ∂₁(v,tf)⊙σₑ(CTf[2],∂₁(u,tf)) + ∂ᵥ(t,tf)⊙σₑ(CTf[3],∂ᵥ(θ,tf)) )*dΩ + # Bending/Axial + Bending
                 ∫( γ(MT,∂₁(v,nf),t) ⊙ σₑ(CTf[4], γ(MT,∂₁(u,nf),θ)) )*dΩ 
                 #∫( γ(MT,∂₁(v,nf),t) ⊙ σₑ(CTf[4], γ(MT,∂₁(u,nf),θ)) )*dΩ 

l((v,t)) = 0

# ARREGLAR THETA AL TALLAT


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)
vh = uh.single_fe_functions[1]
th = uh.single_fe_functions[2]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u" =>vh,
                     "θ" =>th])


#----------------------------------


A((u,θ),(v,t)) = ∫( ∂₁(v,tf)⊙σₑ(CTf[1],∂₁(u,tf)) + ∂ᵥ(t,tf)⊙σₑ(CTf[2],∂ᵥ(θ,tf)) )*dΩ
D((u,θ),(v,t)) = ∫( ∂₁(v,tf)⊙σₑ(CTf[2],∂₁(u,tf)) + ∂ᵥ(t,tf)⊙σₑ(CTf[3],∂ᵥ(θ,tf)) )*dΩ
S((u,θ),(v,t)) = ∫( γ(MT,∂₁(v,nf),t) ⊙ σₑ(CTf[4], γ(MT,∂₁(u,nf),θ)) )*dΩ
#S((u,θ),(v,t)) = ∫( γ(MT,∂₁(v,nf),∂ₗ(t,tf)) ⊙ σₑ(CTf[4], γ(MT,∂₁(u,nf),∂ₗ(θ,tf))) )*dΩ

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

