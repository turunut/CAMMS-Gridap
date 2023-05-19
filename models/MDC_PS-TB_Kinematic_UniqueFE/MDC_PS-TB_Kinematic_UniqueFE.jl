
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

boundary_tags = ["right", "right_points"]
Γ  = BoundaryTriangulation(model,tags=boundary_tags)
dΓ = Measure(Γ,degree)
n_Γ = get_normal_vector(Γ)


#--------------------------------------------------


reffe   = ReferenceFE(lagrangian,VectorValue{2,Float64},order)

Vu = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["left", "left_points"],
                 dirichlet_masks=[(true,true), (true,true)])

Vλ = ConstantFESpace(model,field_type=VectorValue{2,Float64})

V  = MultiFieldFESpace([Vu,Vλ])

g0(x) = VectorValue(0.0,0.0)

Uu = TrialFESpace(Vu,[g0,g0])
Uλ = TrialFESpace(Vλ)

U   = MultiFieldFESpace([Uu,Uλ])


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
L_fun = sum(∫(1.0)*dΓ)                # Omega o Gamma
L     = L_fun

g(x) = VectorValue(0.0,1.0)

z_coord(x) = x[2]
z_cf = CellField(z_coord,Ω)

get_x(x) = x[1]
get_y(x) = x[2]

aΩ((u,λ),(v,μ)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

#aΓ((u,λ),(v,μ)) = ∫( get_x∘(λ)*(v⋅VectorValue(1.0,0.0)) + get_x∘(μ)*(u⋅VectorValue(1.0,0.0)) )*dΓ + 
#                  ∫( get_y∘(λ)*(v⋅VectorValue(0.0,1.0)) + get_y∘(μ)*(u⋅VectorValue(0.0,1.0)) )*dΓ
aΓ((u,λ),(v,μ)) = ∫( (λ⋅v) + (μ⋅u) )*dΓ

 a((u,λ),(v,μ)) = aΩ((u,λ),(v,μ)) + aΓ((u,λ),(v,μ))

l((v,μ)) = ∫(v⋅f)*dΩ + ∫(μ⋅g)*dΓ


#--------------------------------------------------

#A((u,λ),(v,μ)) = ∫( (λ⋅v) + (μ⋅u) )*dΓ
#
#UU = get_trial_fe_basis(U)
#VV = get_fe_basis(V)
#contrA = A(UU,VV)
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
