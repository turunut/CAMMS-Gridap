
#Exemple LINEAR ELASTIC generant malla / Dirichelt
#
#Exemple Linear elastic en que generem un model cuadrat de forma
#procedural, el mallem i definim condicions de Dirichlet imposant
#desplacaments als seus extrems esquerra i dret.

push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/Materials")

using Gridap
using GridapGmsh
using GridapGiD
using Gridap.Geometry
using modCT
using modModel
using modSubroutines
using modInterface
using Gridap.TensorValues
using Gridap.Arrays

prblName = "MDC_PS-TB_McCune_UniqueFE"
projFldr = pwd()

model = GmshDiscreteModel( projFldr*"/models/"*prblName*"/"*prblName*".msh" )
#model = GiDDiscreteModel( projFldr*"/models/"*prblName*"/"*prblName )

order  = 2
degree = 2*order

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)

# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


#--------------------------------------------------


modlType = PlaneStress()

CT1 = CT_Isotrop(1000, 0.2)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop( 100, 0.4)
ct2 = modModel.computeCT(modlType, CT2)

CT3 = CT_Isotrop( 300, 0.3)
ct3 = modModel.computeCT(modlType, CT3)


#--------------------------------------------------


dimens  = 2
#matFlag = ["mat_1", "mat_2", "mat_3"]
matFlag = ["low", "mid", "top"]

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1, ct2, ct3) # Posem un sobre els altres [[CT₁],
                          #                            [CT₂],
                          #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


boundary_left  = ["left","left_points"]
#add_tag_from_tags!(labels,"wall_left",[7])
Γa = BoundaryTriangulation(model,tags=boundary_left)
boundary_right = ["right"]
#add_tag_from_tags!(labels,"wall_right",[8])
Γb = BoundaryTriangulation(model,tags=boundary_right)

Ef = get_E_CellField([CT1, CT2, CT3], tags, Ω)
z_coord(x) = x[2]; zf = CellField(z_coord,Ω)

intrfA = Intrf_Timoshenko(Γa, Ω, degree, Ef, zf)
intrfB = Intrf_Timoshenko(Γb, Ω, degree, Ef, zf)

#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},order)

Vu = TestFESpace(model,reffe;
                 conformity=:H1)
Uu = TrialFESpace(Vu)

Vλa, Uλa = get_test_trial_spaces(intrfA, model)
Vλb, Uλb = get_test_trial_spaces(intrfB, model)

V = MultiFieldFESpace([Vu,Vλa,Vλb])
U = MultiFieldFESpace([Uu,Uλa,Vλb])

u,λ,_ = get_trial_fe_basis(U)
v,μ,_ = get_fe_basis(U)


#--------------------------------------------------

# External forces
f(x) = VectorValue(0.0,0.0)

## Valor equivalent a l'integral
#L_fun = sum(∫(1.0)*dΓ)                # Omega o Gamma
#L     = L_fun

aΩ((u,λ,α),(v,μ,β)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

#aΓ((u,λ),(v,μ)) = ∫( get_x∘(λ)*(v⋅VectorValue(1.0,0.0)) + get_x∘(μ)*(u⋅VectorValue(1.0,0.0)) )*dΓ + 
#                  ∫( get_y∘(λ)*(v⋅VectorValue(0.0,1.0)) + get_y∘(μ)*(u⋅VectorValue(0.0,1.0)) )*dΓ
#aΓ((u,λ),(v,μ)) = ∫( (λ⋅v) + (μ⋅u) )*intrf.dΓ

aΓa((u,λ,α),(v,μ,β)) = contribute_matrix(intrfA, (u,λ,α), (v,μ,β), 1, 2)
aΓb((u,λ,α),(v,μ,β)) = contribute_matrix(intrfB, (u,λ,α), (v,μ,β), 1, 3)

a((u,λ,α),(v,μ,β)) = aΩ((u,λ,α),(v,μ,β)) + aΓa((u,λ,α),(v,μ,β)) + aΓb((u,λ,α),(v,μ,β))

ga(x) = VectorValue(0.0,0.0,0.0)
gb(x) = VectorValue(1.0,0.0,0.0)
l((v,μ,β)) = ∫(v⋅f)*dΩ + ∫(μ⋅ga)*intrfA.dΓ + ∫(β⋅gb)*intrfB.dΓ


#--------------------------------------------------

#A((u,λ),(v,μ)) = ∫( (λ⋅v) + (μ⋅u) )*dΓ
#
##
#Uα = get_trial_fe_basis(U)
#Vα = get_fe_basis(V)
#A(UU,VV) = ∫( (UU[2]⋅VV[1]) + (VV[2]⋅UU[1]) )*dΓ
#B(UU,VV) = ∫( ∂(VV[1])⊙σ(CTf[1],∂(UU[1])) )*dΩ
#contr = aΓa((Uα),(Vα))
#element = first(contr.dict).second[1]

op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

sol = solve(op)

uh = sol.single_fe_functions[1]
writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])
