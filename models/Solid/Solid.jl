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

prblName = "Solid"
projFldr = pwd()

n = 5
r = 12
domain = (0,r,0,1,0,1)
partition = (r*n,n,n)
model = CartesianDiscreteModel(domain,partition)

order  = 1
degree = 2*order

#writevtk(model,"model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"supportA",[1,3,5,7,13,15,17,19,25])
add_tag_from_tags!(labels,"supportB",[2,4,6,8,14,16,18,20,26])
add_tag_from_tags!(labels,"supports",["supportA","supportB"])


#--------------------------------------------------


modlType = Solid()

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
V0 = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["supportA", "supportB"],
                 dirichlet_masks=[(true,true,true), (true,true,true)])

g1(x) = VectorValue(0.0,0.0,0.0)
g2(x) = VectorValue(0.0,0.0,1.0)

U = TrialFESpace(V0,[g1,g2])


#--------------------------------------------------


# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[CT₁],
                          #                            [CT₂],
                          #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


a(u,v) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ
l(v)   = 0


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V0)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])
