#=
Exemple LINEAR ELASTIC generant malla / Dirichelt

Exemple Linear elastic en que generem un model cuadrat de forma
procedural, el mallem i definim condicions de Dirichlet imposant
desplacaments als seus extrems esquerra i dret.
=#
push!(LOAD_PATH, pwd()*"\\src")
push!(LOAD_PATH, pwd()*"\\src\\Materials")

using Gridap
using GridapGmsh
using Gridap.Geometry
using modCT
using modModel
using modSubroutines
using Gridap.TensorValues
using Gridap.Arrays

prblName = "Timoshenko"
projFldr = pwd()

n = 8
domain = (0,100)
partition = (n)
model = CartesianDiscreteModel(domain,partition)

order = 1
degree = 2*order

#writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)


#--------------------------------------------------


MT = Timoshenko(1.0, 5.0, -5.0)

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(MT, CT1)

CT2 = CT_Isotrop(7200, 0.3)
ct2 = modModel.computeCT(MT, CT2)


#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{1,Float64},order)
Vv = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["tag_1", "tag_2"],
                 dirichlet_masks=[(true), (true)])
Vw = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["tag_1", "tag_2"],
                 dirichlet_masks=[(true), (true)])
Vt = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["tag_1", "tag_2"],
                 dirichlet_masks=[(true), (false)])
V = MultiFieldFESpace([Vv,Vw,Vt])

g0(x) = VectorValue( 0.0)
g1(x) = VectorValue(-1.0)

Uv = TrialFESpace(Vv, [g0,g0])
Uw = TrialFESpace(Vw, [g0,g1])
Ut = TrialFESpace(Vt, [g0,g0])
U = MultiFieldFESpace([Uv,Uw,Ut])


#--------------------------------------------------


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


a((u,ω,θ),(v,w,t)) = ∫( ∂(v) ⊙ σₑ(CTf[1],∂(u)) + ∂(t) ⊙ σₑ(CTf[2],∂(θ)) )*dΩ + # Axial         + Axial/Bending
                     ∫( ∂(v) ⊙ σₑ(CTf[2],∂(u)) + ∂(t) ⊙ σₑ(CTf[3],∂(θ)) )*dΩ + # Bending/Axial + Bending
                     ∫( γ(MT,∇(w),t) ⊙ σₑ(CTf[4],γ(MT,∇(ω),θ)) )*dΩ        # Shear

l((v,w,t)) = 0


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
                     "θ" =>th,
                     "∂ω"=>∂(wh), # "∂ω"=>ε(wh)⋅VectorValue(1.0),
                     "γ" =>γ(MT,wh,th)])


#----------------------------------


A((u,ω,θ),(v,w,t)) = ∫( ∂(v)⊙σₑ(CTf[1],∂(u)) )*dΩ
D((u,ω,θ),(v,w,t)) = ∫( ∂(t)⊙σₑ(CTf[3],∂(θ)) )*dΩ
S((u,ω,θ),(v,w,t)) = ∫( γ(MT,∇(w),t) ⊙ σₑ(CTf[4],γ(MT,∇(ω),θ)) )*dΩ 

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

