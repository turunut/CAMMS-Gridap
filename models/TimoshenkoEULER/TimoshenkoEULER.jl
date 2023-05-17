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
                 dirichlet_masks=[(true), (false)])

reffe_q = ReferenceFE(lagrangian,VectorValue{1,Float64},0)
Vq = TestFESpace(model,reffe_q,conformity=:L2)
V = MultiFieldFESpace([Vw,Vt,Vq])

g0(x) = VectorValue( 0.0)
g1(x) = VectorValue(-1.0)

Uw = TrialFESpace(Vw, [g0,g1])
Ut = TrialFESpace(Vt, [g0,g0])
Uq = TrialFESpace(Vq)
U = MultiFieldFESpace([Uw,Ut,Uq])


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


a((ω,θ,ϙ),(w,t,q)) = ∫( ∂(t)⊙σ(CTf[3],∂(θ)) )*dΩ                                   - ∫(    t⊙τ(CTf[4],ϙ)        )*dΩ + 
                                                                                     ∫( ∂(w)⊙τ(CTf[4],toten(ϙ)) )*dΩ -
                     ∫(    q⊙τ(CTf[4],θ)    )*dΩ + ∫( q⊙τ(CTf[4],tovec(∂(ω))) )*dΩ

#a((ω,θ,ϙ),(w,t,q)) = ∫( ∂(t)⊙(M∘∂(θ)) )*dΩ                       - ∫(    t⊙(Q∘ϙ) )*dΩ +
#                                                                   ∫( ∂(w)⊙(Q∘ϙ) )*dΩ -
#                     ∫(    q⊙(Q∘θ)    )*dΩ + ∫( q⊙(Q∘∂(ω)) )*dΩ

l((w,t,q)) = 0


#vω = get_trial_fe_basis(Uw)
#vθ = get_trial_fe_basis(Ut)
#vϙ = get_trial_fe_basis(Uq)
#vw = get_fe_basis(Vw)
#vt = get_fe_basis(Vt)
#vq = get_fe_basis(Vq)
#ca  = ∂(vt)⊙σ(CTf[3],∂(vθ))
#cc1 = vt⊙τ(CTf[4],vϙ)
#cb1 = ∂(vw)⊙τ(CTf[4],toten(vϙ))
#cc2 = vq⊙τ(CTf[4],vθ)
#cb2 = vq⊙τ(CTf[4],tovec(∂(vω)))
#
#quad = dΩ.quad
#pts = get_cell_points(quad)
#ca(pts)[1]
#cc1(pts)[1]
#cb1(pts)[1]
#cc2(pts)[1]
#cb2(pts)[1]


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)
wh = uh.single_fe_functions[1]
th = uh.single_fe_functions[2]
qh = uh.single_fe_functions[3]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["ω" =>wh,
                     "∂ω"=>∂(wh),
                     "θ" =>th,
                     "ϙ" =>qh,
                     "γ" =>γ(MT,wh,th)])

