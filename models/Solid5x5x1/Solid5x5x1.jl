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

prblName = "Solid5x5x1"
projFldr = pwd()

domain = (0,4,0,4,-0.5,0.5)
partition = (5,5,1)
model = CartesianDiscreteModel(domain,partition)

order  = 2
degree = 2*order

#writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"supportA",[1,3,5,7,13,15,17,19,25])
add_tag_from_tags!(labels,"supportB",[2,4,6,8,14,16,18,20,26])
add_tag_from_tags!(labels,"supports",["supportA","supportB"])


#--------------------------------------------------


modlType = Solid()

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------


# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)

boundary_tags = ["supportA"]
ΓA  = BoundaryTriangulation(model,tags=boundary_tags)
dΓA = Measure(ΓA,degree)

boundary_tags = ["supportB"]
ΓB  = BoundaryTriangulation(model,tags=boundary_tags)
dΓB = Measure(ΓB,degree)


#--------------------------------------------------


reffeu = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
V0 = TestFESpace(model,reffeu;conformity=:H1)
U0 = TrialFESpace(V0)

reffeλ = ReferenceFE(lagrangian,VectorValue{3,Float64},(1,0))
VA = FESpace(ΓA,reffeλ)
UA = TrialFESpace(VA)

VB = FESpace(ΓB,reffeλ)
UB = TrialFESpace(VB)

V = MultiFieldFESpace([V0,VA,VB])
U = MultiFieldFESpace([U0,UA,UB])


#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[CT₁],
                          #                            [CT₂],
                          #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


a((u,λ₁,λ₂),(v,μ₁,μ₂)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ  +
                         ∫( (λ₁⋅v) + (μ₁⋅u)     )*dΓA +
                         ∫( (λ₂⋅v) + (μ₂⋅u)     )*dΓB

g1 = VectorValue(0.0,0.0,0.0)
g2 = VectorValue(1.0,0.0,0.0)

l((v,μ₁,μ₂)) = ∫( (μ₁⋅g1) )*dΓA +
               ∫( (μ₂⋅g2) )*dΓB

A = assemble_matrix(a,U,V)

λa((u,λ₁,λ₂),(v,μ₁,μ₂)) = ∫( (λ₂⋅v) + (μ₂⋅u) )*dΓB
λA = assemble_matrix(λa,U,V)

for i in 1:1089
  if abs( λA[1093,i] ) > 0.000001; println(λA[1093,i]); end 
end

for i in 1:1089
  if abs( λA[1120,i] ) > 0.000001; println(λA[1120,i]); end 
end


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

xh = solve(op);
uh, λ1, λ2 = xh;

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])
