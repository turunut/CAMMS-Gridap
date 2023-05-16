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

prblName = "MDC_PS-TB_Kinematic"
projFldr = pwd()

model = GmshDiscreteModel( projFldr*"/models/"*prblName*"/"*prblName*".msh" )

order  = 1
degree = 2*order

#writevtk(model,"model")

labels = get_face_labeling(model)


#--------------------------------------------------


modlType = PlaneStress()

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
Vu = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["left", "left_points"],
                 dirichlet_masks=[(true,true), (true,true)])
Vλ0 = ConstantFESpace(model)
Vλ1 = ConstantFESpace(model)
Vλ2 = ConstantFESpace(model)
Vλ3 = ConstantFESpace(model)
V  = MultiFieldFESpace([Vu,Vλ0,Vλ1,Vλ2,Vλ3])

g1(x) = VectorValue(0.0,0.0)
g2(x) = VectorValue(0.1,0.0)

Uu  = TrialFESpace(Vu,[g1,g1])
Uλ0 = TrialFESpace(Vλ0)
Uλ1 = TrialFESpace(Vλ1)
Uλ2 = TrialFESpace(Vλ2)
Uλ3 = TrialFESpace(Vλ3)
U   = MultiFieldFESpace([Uu,Uλ0,Uλ1,Uλ2,Uλ3])


#--------------------------------------------------


boundary_tags = ["right", "right_points"]
Γ  = BoundaryTriangulation(model,tags=boundary_tags)
dΓ = Measure(Γ,degree)
n_Γ = get_normal_vector(Γ)

# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


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

# Valor equivalent a l'integral
L__fun = sum(∫(    1.0 )*dΓ)                # Omega o Gamma
L      = L__fun

g0(x) = 00.0/L
g1(x) = 50.0/L
g2(x) = 00.0/L
g3(x) = 00.0/L

## Ordre zero
#a((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ + ∫( λ*(v⋅n_Γ) )*dΓ +
#                 ∫( μ*(u⋅n_Γ) )*dΓ

z_coord(x) = x[2]
z_cf = CellField(z_coord,Ω)


 aΩ((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ
aΓ0((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) = ∫( λ0*(ctˣ(z_cf,0)*v⋅n_Γ) )*dΓ + ∫( μ0*(ctˣ(z_cf,0)*u⋅n_Γ) )*dΓ
aΓ1((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) = ∫( λ1*(ctˣ(z_cf,1)*v⋅n_Γ) )*dΓ + ∫( μ1*(ctˣ(z_cf,1)*u⋅n_Γ) )*dΓ
aΓ2((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) = ∫( λ2*(ctˣ(z_cf,2)*v⋅n_Γ) )*dΓ + ∫( μ2*(ctˣ(z_cf,2)*u⋅n_Γ) )*dΓ
aΓ3((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) = ∫( λ3*(ctˣ(z_cf,3)*v⋅n_Γ) )*dΓ + ∫( μ3*(ctˣ(z_cf,3)*u⋅n_Γ) )*dΓ

  a((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) =  aΩ((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) + 
                                       aΓ0((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) + 
                                       aΓ1((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) + 
                                       aΓ2((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3)) + 
                                       aΓ3((u,λ0,λ1,λ2,λ3),(v,μ0,μ1,μ2,μ3))

l((v,μ0,μ1,μ2,μ3)) = ∫(v⋅f)*dΩ +
                     ∫(μ0*g0)*dΓ +
                     ∫(μ1*g1)*dΓ +
                     ∫(μ2*g2)*dΓ +
                     ∫(μ3*g3)*dΓ


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

sol = solve(op)

uh = sol.single_fe_functions[1]
writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])
