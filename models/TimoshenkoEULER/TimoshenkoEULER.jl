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

prblName = "TimoshenkoEULER"
projFldr = pwd()

n = 1
domain = (0,2000)
partition = (n)
model = CartesianDiscreteModel(domain,partition)

order = 1

#writevtk(model,"model")

labels = get_face_labeling(model)


#--------------------------------------------------

MT = Timoshenko(1.0, 50.0, -50.0)

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(MT, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(MT, CT2)


#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{1,Float64},order)
Vw = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["tag_1", "tag_2"],
                 dirichlet_masks=[(true), (true)])
Vt = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["tag_1", "tag_2"],
                 dirichlet_masks=[(true), (true)])

reffe_q = ReferenceFE(lagrangian,VectorValue{1,Float64},0)
Vq = TestFESpace(model,reffe_q,conformity=:L2)
V = MultiFieldFESpace([Vw,Vt,Vq])

g0(x) = VectorValue(  0.0)
g1(x) = VectorValue(-10.0)

Uw = TrialFESpace(Vw, [g0,g0])
Ut = TrialFESpace(Vt, [g0,g1])
Uq = TrialFESpace(Vq)
U = MultiFieldFESpace([Uw,Ut,Uq])


#--------------------------------------------------


degree = 1*order
# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 1 # linies
matFlag = ["top", "mid", "low"]

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1, ct2) # Posem un sobre els altres [[A₁,B₁,D₁,S₁],
                     #                            [A₂,B₂,D₂,S₂]]

CTf = get_CT_CellField(MT, CTs, tags, Ω)


#--------------------------------------------------


MTm = PlaneStress()

a((ω,θ,ϙ),(w,t,q)) = ∫( εₘ(t)⊙σₘ(MTm,CTf[3],εₘ(θ)) )*dΩ                                 - ∫(     t⊙σₘ(MTm,CTf[4],ϙ) )*dΩ +
                                                                                         ∫( εₘ(w)⊙σₘ(MTm,CTf[4],ϙ) )*dΩ -
                     ∫(     q⊙σₘ(MTm,CTf[4],θ)     )*dΩ + ∫( q⊙σₘ(MTm,CTf[4],εₘ(ω)) )*dΩ

l((w,t,q)) = 0


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)
wh = uh.single_fe_functions[1]
th = uh.single_fe_functions[2]
qh = uh.single_fe_functions[3]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["ω"=>wh,
                     "θ"=>th,
                     "ϙ"=>qh,
                     "γ" =>γₘ(modlType,wh,th)])

