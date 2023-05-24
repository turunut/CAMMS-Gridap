
#Exemple LINEAR ELASTIC generant malla / Dirichelt
#
#Exemple Linear elastic en que generem un model cuadrat de forma
#procedural, el mallem i definim condicions de Dirichlet imposant
#desplacaments als seus extrems esquerra i dret.

push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/Materials")

using Gridap
using GridapGmsh
using Gridap.Geometry
using modCT
using modModel
using modSubroutines
using modInterface
using Gridap.TensorValues
using Gridap.Arrays

prblName = "MDC_PS-TB_Kinematic_UniqueFE"
projFldr = pwd()

model = GmshDiscreteModel( projFldr*"/models/"*prblName*"/"*prblName*".msh" )

order  = 2
degree = 2*order

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)


#--------------------------------------------------


modlType = PlaneStress()

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------


# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)

boundary_tags_a = ["left", "left_points"]
Γa = BoundaryTriangulation(model,tags=boundary_tags_a)
boundary_tags_b = ["right", "right_points"]
Γb = BoundaryTriangulation(model,tags=boundary_tags_b)

intrfA = Intrf_Kinematic2D(Γa, degree)
intrfB = Intrf_Kinematic2D(Γb, degree)


#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},order)

Vu = TestFESpace(model,reffe;
                 conformity=:H1)
Uu = TrialFESpace(Vu)

Vλa, Uλa = get_test_trial_spaces(intrfA, model)
Vλb, Uλb = get_test_trial_spaces(intrfB, model)

V = MultiFieldFESpace([Vu,Vλa,Vλb])
U = MultiFieldFESpace([Uu,Uλa,Vλb])


#--------------------------------------------------


dimens  = 2
matFlag = ["top", "mid", "low"]

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1, ct2, ct1) # Posem un sobre els altres [[CT₁],
                          #                            [CT₂],
                          #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------

# External forces
f(x) = VectorValue(0.0,0.0)

aΩ((u,λ,α),(v,μ,β)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

#aΓ((u,λ),(v,μ)) = ∫( get_x∘(λ)*(v⋅VectorValue(1.0,0.0)) + get_x∘(μ)*(u⋅VectorValue(1.0,0.0)) )*dΓ + 
#                  ∫( get_y∘(λ)*(v⋅VectorValue(0.0,1.0)) + get_y∘(μ)*(u⋅VectorValue(0.0,1.0)) )*dΓ
#aΓ((u,λ),(v,μ)) = ∫( (λ⋅v) + (μ⋅u) )*intrf.dΓ

aΓa((u,λ,α),(v,μ,β)) = contribute_matrix(intrfA, (u,λ,α),(v,μ,β), 1, 2)
aΓb((u,λ,α),(v,μ,β)) = contribute_matrix(intrfB, (u,λ,α),(v,μ,β), 1, 3)

a((u,λ,α),(v,μ,β)) = aΩ((u,λ,α),(v,μ,β)) + aΓa((u,λ,α),(v,μ,β)) + aΓb((u,λ,α),(v,μ,β))

ga(x) = VectorValue(0.0,1.0)
gb(x) = VectorValue(0.0,0.0)

la((v,μ,β)) = contribute_vector(intrfA, (v,μ,β), 2, ga)
lb((v,μ,β)) = contribute_vector(intrfB, (v,μ,β), 3, gb)

l((v,μ,β)) = ∫(v⋅f)*dΩ + la((v,μ,β)) + lb((v,μ,β))


#--------------------------------------------------

#A((u,λ),(v,μ)) = ∫( (λ⋅v) + (μ⋅u) )*dΓ
#
##
#UU = get_trial_fe_basis(U)
#VV = get_fe_basis(V)
#A(UU,VV) = ∫( (UU[2]⋅VV[1]) + (VV[2]⋅UU[1]) )*dΓ
#B(UU,VV) = ∫( ∂(VV[1])⊙σ(CTf[1],∂(UU[1])) )*dΩ
#contrA = B(UU,VV)
#elementA = first(contrA.dict).second[1]

op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

sol = solve(op)

uh = sol.single_fe_functions[1]
writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])
