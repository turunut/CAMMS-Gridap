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

prblName = "PlaneStress1x1"
projFldr = pwd()

n = 1
d = 10
domain = (0,d,0,d)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

order = 1
degree = 2

#writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"diri_l",[7, 1, 3]) # left  edge
add_tag_from_tags!(labels,"diri_r",[4]) # right edge


#--------------------------------------------------


modlType = PlaneStress()

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(modlType, CT1)


#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
V0 = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["diri_l","diri_r" ],
                 dirichlet_masks=[(true,true), (true,true)])

g1(x) = VectorValue(0.0,0.0)
g2(x) = VectorValue(0.1,0.0)

U = TrialFESpace(V0,[g1,g2])


#--------------------------------------------------


# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 2
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[CT₁],
                          #                            [CT₂],
                          #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


a(u,v) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ
l(v)   = 0





UU = get_trial_fe_basis(U)
VV = get_fe_basis(V0)
contrA = a(UU,VV)
elementA = first(contrA.dict).second[1]

println(elementA[1,1])
println(elementA[1,2])


#--------------------------------------------------


#op = AffineFEOperator(a,l,U,V0)
#
#ls = LUSolver()
#solver = LinearFESolver(ls)
#
#uh = solve(op)
#
#writevtk(Ω,"models/"*prblName*"/"*prblName,
#         cellfields=["u"=>uh,
#                     "ε"=>∂(uh),
#                     "σ"=>σ(CTf[1],∂(uh))])
