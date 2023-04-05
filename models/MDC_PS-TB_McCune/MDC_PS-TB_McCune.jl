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

prblName = "MDC_PS-TB_McCune"
projFldr = pwd()

model = GmshDiscreteModel( projFldr*"/models/"*prblName*"/"*prblName*".msh" )

order  = 2
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
VU = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["left", "left_points"],
                 dirichlet_masks=[(true,true), (true,true)])
Vu = ConstantFESpace(model)
Vθ = ConstantFESpace(model)
Vω = ConstantFESpace(model)
V = MultiFieldFESpace([VU,Vu,Vθ,Vω])

g1(x) = VectorValue(0.0,0.0)
g2(x) = VectorValue(0.1,0.0)

UU = TrialFESpace(VU,[g1,g1])
Uu = TrialFESpace(Vu)
Uθ = TrialFESpace(Vθ)
Uω = TrialFESpace(Vω)
U   = MultiFieldFESpace([UU,Uu,Uθ,Uω])


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
u_beam(x) = 0.0
θ_beam(x) = 0.0
ω_beam(x) = 1.0

get_x(x) = x[1]
get_y(x) = x[2]

Ef = get_E_CellField([CT1, CT2, CT1], tags, Ω)

z_coord(x) = x[2]
z_cf = CellField(z_coord,Ω)

Da_fun(Ef)      = sum(∫(           Ef )*dΓ) # Omega o Gamma
Db_fun(Ef,z_cf) = sum(∫(      z_cf*Ef )*dΓ) # Omega o Gamma
Dd_fun(Ef,z_cf) = sum(∫( z_cf*z_cf*Ef )*dΓ) # Omega o Gamma
S__fun(Ef)      = sum(∫(           Ef )*dΓ) # Omega o Gamma
L__fun = sum(∫(    1.0 )*dΓ)                # Omega o Gamma
I__fun = sum(∫( z_cf*z_cf )*dΓ)             # Omega o Gamma

function step(z::Float64,z_val::Float64)
  if z <= (z_val)
    return 1.0
  else
    return 0.0
  end
end

step_field(z_cf,z_val) = CellField(step.(z_cf,z_val),Ω)

da_fun(Ef,z_cf,z_val) = sum(∫(      step_field(z_cf,z_val)*Ef )*dΓ)
db_fun(Ef,z_cf,z_val) = sum(∫( z_cf*step_field(z_cf,z_val)*Ef )*dΓ)

Da = Da_fun(Ef)
Db = Db_fun(Ef,z_cf)
Dd = Dd_fun(Ef,z_cf)
da(z_val) = da_fun(Ef,z_cf,z_val)
db(z_val) = db_fun(Ef,z_cf,z_val)
L  = L__fun
I  = I__fun


a((U,u,θ,ω),(V,v,t,w)) =         ∫(          ∂(V) ⊙ σ(CTf[1],∂(U) ) )*dΩ  + (L/Da)*(∫( v*( (Ef)*(get_x∘(U⋅n_Γ)) ) )*dΓ) + (L/Dd)*(∫( t*( (z_cf*Ef)*(get_x∘(U⋅n_Γ)) ) )*dΓ) + (L/Dd)*(∫( w*( (db∘z_cf)*(get_y∘U) ) )*dΓ) +
                         (L/Da)*(∫( u*(      (Ef)*(get_x∘(V⋅n_Γ)) ) )*dΓ) +
                         (L/Dd)*(∫( θ*( (z_cf*Ef)*(get_x∘(V⋅n_Γ)) ) )*dΓ) +
                         (L/Dd)*(∫( ω*(     (db∘z_cf)*(get_y∘(V)) ) )*dΓ)

l((V,v,t,w)) = ∫(f⋅V)*dΩ + ∫(u_beam*v)*dΓ + ∫(θ_beam*t)*dΓ + ∫(ω_beam*w)*dΓ


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

sol = solve(op)

Uh = sol.single_fe_functions[1]
uh = sol.single_fe_functions[2]
θh = sol.single_fe_functions[3]
ωh = sol.single_fe_functions[4]
writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>Uh,
                     "ε"=>∂(Uh),
                     "σ"=>σ(CTf[1],∂(Uh))])
