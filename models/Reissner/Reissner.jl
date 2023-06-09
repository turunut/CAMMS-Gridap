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

prblName = "Reissner"
projFldr = pwd()

n = 1
domain = (0,4,0,4)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

order = 1
degree = 2*order

#writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"diri_l",[7, 1, 3]) # left  edge
add_tag_from_tags!(labels,"diri_r",[4]) # right edge


#--------------------------------------------------


MT = Reissner(0.5, -0.5)

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(MT, CT1)


#--------------------------------------------------


reffe2 = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
reffe1 = ReferenceFE(lagrangian,VectorValue{1,Float64},order)
Vv = TestFESpace(model,reffe2;
                 conformity=:H1,
                 dirichlet_tags=["diri_l"],
                 dirichlet_masks=[(true,true)])
Vw = TestFESpace(model,reffe1;
                 conformity=:H1,
                 dirichlet_tags=["diri_l","diri_r"],
                 dirichlet_masks=[(true),(true)])
Vt = TestFESpace(model,reffe2;
                 conformity=:H1,
                 dirichlet_tags=["diri_l"],
                 dirichlet_masks=[(true,true)])
V = MultiFieldFESpace([Vv,Vw,Vt])

g0_1(x) = VectorValue( 0.0)
g1_1(x) = VectorValue(+1.0)
g0_2(x) = VectorValue( 0.0, 0.0)
g1_2(x) = VectorValue(+1.0, 0.0)

Uv = TrialFESpace(Vv, [g0_2])
Uw = TrialFESpace(Vw, [g0_1,g1_1])
Ut = TrialFESpace(Vt, [g0_2])
U = MultiFieldFESpace([Uv,Uw,Ut])


###--------------------------------------------------
##
##
##reffe2 = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
##reffe1 = ReferenceFE(lagrangian,VectorValue{1,Float64},order)
##Vv = TestFESpace(model,reffe2;
##                 conformity=:H1,
##                 dirichlet_tags=["boundary"],
##                 dirichlet_masks=[(true,true)])
##Vw = TestFESpace(model,reffe1;
##                 conformity=:H1,
##                 dirichlet_tags=["boundary"],
##                 dirichlet_masks=[(true)])
##Vt = TestFESpace(model,reffe2;
##                 conformity=:H1,
##                 dirichlet_tags=["boundary"],
##                 dirichlet_masks=[(true,true)])
##V = MultiFieldFESpace([Vv,Vw,Vt])
##
##g0_1(x) = VectorValue( 0.0)
##g1_1(x) = VectorValue(+1.0)
##g0_2(x) = VectorValue( 0.0, 0.0)
##g1_2(x) = VectorValue(+1.0, 0.0)
##
##Uv = TrialFESpace(Vv, [g0_2])
##Uw = TrialFESpace(Vw, [g0_1])
##Ut = TrialFESpace(Vt, [g0_2])
##U = MultiFieldFESpace([Uv,Uw,Ut])
##
##
###--------------------------------------------------


# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 2 # surfaces
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[A₁,B₁,D₁,S₁],
                     #                            [A₂,B₂,D₂,S₂]]

CTf = get_CT_CellField(MT, CTs, tags, Ω)


#--------------------------------------------------


q(x) = VectorValue(-1.0)

a((u,ω,θ),(v,w,t)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) + ∂(t)⊙σ(CTf[2],∂(θ)) )*dΩ + # Axial         + Axial/Bending
                     ∫( ∂(v)⊙σ(CTf[2],∂(u)) + ∂(t)⊙σ(CTf[3],∂(θ)) )*dΩ + # Bending/Axial + Bending
                     ∫( γ(MT,∇(w),t)⊙σₑ(CTf[4], γ(MT,∇(ω),θ)) )*dΩ        # Shear

l((v,w,t)) = 0.0
#l((v,w,t)) = ∫(w⋅q)*dΩ


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)
vh = uh.single_fe_functions[1]
wh = uh.single_fe_functions[2]
th = uh.single_fe_functions[3]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u" =>vh,
                     "ω" =>wh,
                     "θ" =>th])
                     #"∂ω"=>∂(wh)]) # "∂ω"=>ε(wh)⋅VectorValue(1.0),
                     #"γ" =>γ(MT,wh,th)])


#--------------------------------------------------


A((u,ω,θ),(v,w,t)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) + ∂(t)⊙σ(CTf[2],∂(θ)) )*dΩ
D((u,ω,θ),(v,w,t)) = ∫( ∂(v)⊙σ(CTf[2],∂(u)) + ∂(t)⊙σ(CTf[3],∂(θ)) )*dΩ
S((u,ω,θ),(v,w,t)) = ∫( γ(MT,∇(w),t)⊙σₑ(CTf[4], γ(MT,∇(ω),θ)) )*dΩ

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

