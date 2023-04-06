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
using Gridap.ReferenceFEs

prblName = "Timoshenko2D"
projFldr = pwd()


#--------------------------------------------------


#https://github.com/gridapapps/GridapGeosciences.jl/blob/master/src/CubedSphereDiscreteModels.jl

nodes = [VectorValue(0.0,0.0),VectorValue(0.0,1.0),VectorValue(1.0,1.0)]
c2n_map = Table([1,2,2,3],[1,3,5])

cell_type = Int8[1,1]
polys = [SEGMENT]
reffes = map(p->LagrangianRefFE(Float64,p,1),polys)
orientation = NonOriented()

topo   = UnstructuredGridTopology(nodes,c2n_map,cell_type,polys,orientation)
grid   = UnstructuredGrid(nodes,c2n_map,reffes,cell_type,orientation)

#d_to_num_dfaces_to_entity = [[1,3,3,3,2],[3,3,3,3,3]]
#tag_to_name   = ["interior","boundary","bc1","bc2"]
#tag_to_entity = [[3],[1,2],[1],[2]]
d_to_dface_to_entity = [[1,3,2],[3,3]]
tag_to_entities = [[3],[1,2],[1],[2]]
tag_to_name     = ["interior","boundary","bc1","bc2"]

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
                 dirichlet_tags=["bc1", "bc2"],
                 dirichlet_masks=[(true,true), (true,true)])
Vt = TestFESpace(model,reffe1;
                 conformity=:H1,
                 dirichlet_tags=["bc1", "bc2"],
                 dirichlet_masks=[(true), (false)])
V = MultiFieldFESpace([Vv,Vt])

g0_1(x) = VectorValue(  0.0)
g1_1(x) = VectorValue(-10.0)
g0_2(x) = VectorValue(0.0, 0.0)
g1_2(x) = VectorValue(0.1, 0.0)

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

CTs = hcat(ct1, ct2) # Posem un sobre els altres [[A₁,B₁,D₁,S₁],
                     #                            [A₂,B₂,D₂,S₂]]

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

function get_tangent_vector(Ω::Triangulation{1})
    cmaps = get_cell_map(Ω)
    return CellField(lazy_map(my_tangent,cmaps),Ω)
end

tf = get_tangent_vector(Ω)

#getL_op(tf) = TensorValue([ tf[1] -tf[2];
#                            tf[2]  tf[1] ])
#getL(tf::CellField) = Operation(getL_op)(tf)
#function g2l(tf, u)
#    eval = getL(tf)⋅u
#    return getL(tf)⋅u # ∘(getL(tf)⋅u)
#end

#g2l(tf,u) = ∘(tf⋅u)
get₁(x) = VectorValue(VectorValue(1.0,0.0)⋅x)
get₂(x) = VectorValue(VectorValue(0.0,1.0)⋅x)
get(x)  = VectorValue(VectorValue(1.0,1.0)⋅x)

function tangential_derivative(u::CellField,tf::CellField)
  gradient(u)⋅tf
end

g2l₁(tf,u) = get₁ ∘ tangential_derivative(u,tf)
g2l₂(tf,u) = get₂ ∘ tangential_derivative(u,tf)

function tangential_derivative_1(u::CellField)
  VectorValue(1.0,1.0)⋅u #⋅TensorValue(VectorValue(1.0,1.0))
end
g2l(u) = tangential_derivative_1(gradient(u)) # get ∘ 

function ∂∂(u)
  return CellField(TensorValue(0.1), Ω)
end

#a((u,θ),(v,t)) = ∫( ∂(g2l(tf,v))⊙σ(CTf[1],∂(g2l(tf,u))) + ∂(t)⊙σ(CTf[2],∂(θ)) )*dΩ + # Axial         + Axial/Bending
#                 ∫( ∂(g2l(tf,v))⊙σ(CTf[2],∂(g2l(tf,u))) + ∂(t)⊙σ(CTf[3],∂(θ)) )*dΩ + # Bending/Axial + Bending
#                 ∫( γ(MT,g2l(tf,v),t) ⊙ σ(CTf[4], γ(MT,g2l(tf,u),θ)) )*dΩ        # Shear

a((u,θ),(v,t)) = ∫( g2l₁(tf,v)⊙σₑ(CTf[1],g2l₁(tf,u)) + g2l(t)⊙σₑ(CTf[2],g2l(θ)) )*dΩ + # Axial         + Axial/Bending
                 ∫( g2l₁(tf,v)⊙σₑ(CTf[2],g2l₁(tf,u)) + g2l(t)⊙σₑ(CTf[3],g2l(θ)) )*dΩ + # Bending/Axial + Bending
                 #∫( γ(MT,g2l₂(tf,v),t) ⊙ σₑ(CTf[4], γ(MT,g2l₂(tf,u),θ)) )*dΩ 
                 ∫( γγ(MT,g2l₂(tf,v),t) ⊙ σₑ(CTf[4], γγ(MT,g2l₂(tf,u),θ)) )*dΩ 

l((v,t)) = 0


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

