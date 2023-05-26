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

prblName = "MDC_PS-TB_McCune_Isolated_ML"
projFldr = pwd()

model = GmshDiscreteModel( projFldr*"/models/"*prblName*"/"*prblName*".msh" )

order  = 2
degree = 2*order

#writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)


#--------------------------------------------------


modlType = PlaneStress()

CT1 = CT_Isotrop(1000, 0.2)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop( 100, 0.4)
ct2 = modModel.computeCT(modlType, CT2)

CT3 = CT_Isotrop( 300, 0.3)
ct3 = modModel.computeCT(modlType, CT3)


#--------------------------------------------------


reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
VU = TestFESpace(model,reffe;
                 conformity=:H1,
                 dirichlet_tags=["left", "left_points"],
                 dirichlet_masks=[(true,true), (true,true)])
Vu = ConstantFESpace(model)
Vθ = ConstantFESpace(model)
V = MultiFieldFESpace([VU,Vu,Vθ])

g1(x) = VectorValue(0.0,0.0)
g2(x) = VectorValue(0.1,0.0)

UU = TrialFESpace(VU,[g1,g1])
Uu = TrialFESpace(Vu)
Uθ = TrialFESpace(Vθ)
U   = MultiFieldFESpace([UU,Uu,Uθ])


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
matFlag = ["low", "mid", "top"]

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1, ct2, ct3) # Posem un sobre els altres [[CT₁],
                          #                            [CT₂],
                          #                            [CT₁]]

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#--------------------------------------------------


# External forces
f(x) = VectorValue(0.0,0.0)
u_beam(x) = 1.0
θ_beam(x) = 0.0

get_x(x) = x[1]
get_y(x) = x[2]

Ef = get_E_CellField([CT1, CT2, CT3], tags, Ω)

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
Dinv = 1/((Dd*Da)-(Db^2))

#for i=-5:0.25:5
#  println(i, " ", da(i))
#end
#
#for i=-5:0.25:5
#  println(i, " ", db(i))
#end

a((U,u,θ),(V,v,t)) =           ∫( ∂(V) ⊙ σ(CTf[1],∂(U) ) )*dΩ  + 
                     (L*Dinv)*(∫( v*( (+Dd)*(Ef)*(get_x∘(U⋅n_Γ)) + (-Db)*(z_cf*Ef)*(get_x∘(U⋅n_Γ)) ) )*dΓ) +
                     (L*Dinv)*(∫( u*( (+Dd)*(Ef)*(get_x∘(V⋅n_Γ)) + (-Db)*(z_cf*Ef)*(get_x∘(V⋅n_Γ)) ) )*dΓ) + 
                     (L*Dinv)*(∫( t*( (-Dd)*(Ef)*(get_x∘(U⋅n_Γ)) + (+Da)*(z_cf*Ef)*(get_x∘(U⋅n_Γ)) ) )*dΓ) +
                     (L*Dinv)*(∫( θ*( (-Dd)*(Ef)*(get_x∘(V⋅n_Γ)) + (+Da)*(z_cf*Ef)*(get_x∘(V⋅n_Γ)) ) )*dΓ)

l((V,v,t)) = ∫(f⋅V)*dΩ + ∫(u_beam*v)*dΓ + ∫(θ_beam*t)*dΓ


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

sol = solve(op)

Uh = sol.single_fe_functions[1]
uh = sol.single_fe_functions[2]
writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>Uh,
                     "ε"=>∂(Uh),
                     "σ"=>σ(CTf[1],∂(Uh))])
