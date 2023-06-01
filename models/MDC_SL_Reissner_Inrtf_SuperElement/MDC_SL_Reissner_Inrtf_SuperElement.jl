push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/Materials")

using Gridap
using Gridap.Geometry
using Gridap.Adaptivity
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.CellData
using modCT
using modModel
using modSubroutines
using modInterface

using GridapGmsh

using FillArrays

prblName = "MDC_SL_Reissner_Inrtf_SuperElement"
projFldr = pwd()

order = 2
degree = 2*order

############################################################################################
# Fine model
domain = (0,4,0,4,-0.5,0.5)
partition = (12,12,6)
model = CartesianDiscreteModel(domain,partition)

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"faceY0",[23])
add_tag_from_tags!(labels,"faceY1",[24])
add_tag_from_tags!(labels,"faceX0",[25])
add_tag_from_tags!(labels,"faceX1",[26])


#--------------------------------------------------


modlType = Solid()

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------


Ω  = Triangulation(model)
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1)

CTf = get_CT_CellField(modlType, CTs, tags, Ω)

#---------------------------

modlType_2D = PlaneStress()

ct1_2D = modModel.computeCT(modlType_2D, CT1, true)
ct2_2D = modModel.computeCT(modlType_2D, CT2, true)

CTs_2D = hcat(ct1_2D)

CTf_2D = get_CT_CellField(modlType_2D, CTs_2D, tags, Ω)

z_coord(x) = x[3]; zf = CellField(z_coord,Ω)


#--------------------------------------------------


Γ₀ = BoundaryTriangulation(model,tags=["tag_23"])
Γ₁ = BoundaryTriangulation(model,tags=["tag_26"])
Γ₂ = BoundaryTriangulation(model,tags=["tag_24"])
Γ₃ = BoundaryTriangulation(model,tags=["tag_25"])


############################################################################################
# Coarse model

Dc      = num_cell_dims(model)
topo    = get_grid_topology(model)
f2n_map = Geometry.get_faces(topo,Dc-1,0)
coords  = Geometry.get_vertex_coordinates(topo)
face_coords = map(Reindex(coords),f2n_map)

tol = 1.e-3
com_coords = lazy_map(N->sum(N)/length(N),face_coords) # center of mass
int_coords = map(N->VectorValue(Int(floor(N[1]/tol)),Int(floor(N[2]/tol)),Int(floor(N[3]/tol))),com_coords)

############################################################################################

axis_id = 2; face_0_pos = minimum(lazy_map(c->c[axis_id],int_coords))
intrf₀  = Intrf_Reissner(Ω, Γ₀, int_coords, axis_id, face_0_pos, degree, true)
intrf₀.CTf_2D = CTf_2D[1]; intrf₀.zf = zf

arrayD = zeros(6,6)
z1=-0.5; z2=0.5
arrayDa = (1/1)*(z2^1-z1^1) *get_array(ct1_2D[1])
arrayDb = (1/2)*(z2^2-z1^2) *get_array(ct1_2D[1])
arrayDd = (1/3)*(z2^3-z1^3) *get_array(ct1_2D[1])
arrayD[1:3,1:3] = arrayDa
arrayD[4:6,1:3] = arrayDb; arrayD[1:3,4:6] = arrayDb
arrayD[4:6,4:6] = arrayDd
invarD = inv(arrayD)
tensorD = TensorValue(arrayD)
invtesD = TensorValue(invarD)

intrf₀.invD = invtesD



dempty_fun(Ef,zf,z_val) = sum(∫( 1.0 )*intrf₀.dΓf)
dempty(z_val) = dempty_fun(intrf₀.CTf_2D,intrf₀.zf,z_val)

da_fun(Ef,zf,z_val) = sum(∫(    step_field(zf,z_val,intrf₀.Ω)*Ef )*intrf₀.dΓ)
da(z_val) = da_fun(intrf₀.CTf_2D,intrf₀.zf,z_val)
da(0.5)/4
Aa = sum( ∫( da∘(zf) )*intrf₀.dΓ )

db_fun(Ef,zf,z_val) = sum(∫( zf*step_field(zf,z_val,intrf₀.Ω)*Ef )*intrf₀.dΓ)
db(z_val) = db_fun(intrf₀.CTf_2D,intrf₀.zf,z_val)
db(0.5)/4
Ab = sum( ∫( db∘(zf) )*intrf₀.dΓ )



axis_id = 1; face_1_pos = maximum(lazy_map(c->c[axis_id],int_coords))
intrf₁  = Intrf_Reissner(Ω, Γ₁, int_coords, axis_id, face_1_pos, degree, false)
intrf₁.CTf_2D = CTf_2D[1]; intrf₁.zf = zf; intrf₁.invD = invtesD

axis_id = 2; face_2_pos = maximum(lazy_map(c->c[axis_id],int_coords))
intrf₂  = Intrf_Reissner(Ω, Γ₂, int_coords, axis_id, face_2_pos, degree, true)
intrf₂.CTf_2D = CTf_2D[1]; intrf₂.zf = zf; intrf₂.invD = invtesD

axis_id = 1; face_3_pos = minimum(lazy_map(c->c[axis_id],int_coords))
intrf₃  = Intrf_Reissner(Ω, Γ₃, int_coords, axis_id, face_3_pos, degree, true)
intrf₃.CTf_2D = CTf_2D[1]; intrf₃.zf = zf; intrf₃.invD = invtesD

############################################################################################
# FESpaces 
# Model 3D
reffe_u  = ReferenceFE(lagrangian,VectorValue{3,Float64},order)


Vu = TestFESpace(model,reffe_u;conformity=:H1)
Uu = TrialFESpace(Vu)

#Vλ_X0, Uλ_X0 = get_test_trial_spaces(intrf_X0, model)
Vλ₀, Uλ₀ = get_test_trial_spaces(intrf₀)
Vλ₁, Uλ₁ = get_test_trial_spaces(intrf₁)
Vλ₂, Uλ₂ = get_test_trial_spaces(intrf₂)
Vλ₃, Uλ₃ = get_test_trial_spaces(intrf₃)

V = MultiFieldFESpace([Vu,Vλ₀,Vλ₁,Vλ₂,Vλ₃])
U = MultiFieldFESpace([Uu,Uλ₀,Uλ₁,Uλ₂,Uλ₃])

## Models linea
Ve₀, Ue₀, reffe_e₀ = get_line_test_trial_spaces(intrf₀,order)
Ve₁, Ue₁, reffe_e₁ = get_line_test_trial_spaces(intrf₁,order)
Ve₂, Ue₂, reffe_e₂ = get_line_test_trial_spaces(intrf₂,order)
Ve₃, Ue₃, reffe_e₃ = get_line_test_trial_spaces(intrf₃,order)


#--------------------------------------------------
  

f = VectorValue(0.0,0.0,0.0)

ue_c₀ = get_line_distribution(3, intrf₀, reffe_e₀, Ue₀, 9)
ue_c₁ = get_line_distribution(3, intrf₁, reffe_e₁, Ue₁, 0)
ue_c₂ = get_line_distribution(3, intrf₂, reffe_e₂, Ue₂, 0)
ue_c₃ = get_line_distribution(3, intrf₃, reffe_e₃, Ue₃, 0)

aΩ((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ

aΓ₀((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) = contribute_matrix(intrf₀, (u,λ₀,λ₁,λ₂,λ₃), (v,μ₀,μ₁,μ₂,μ₃), 1, 2)
aΓ₁((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) = contribute_matrix(intrf₁, (u,λ₀,λ₁,λ₂,λ₃), (v,μ₀,μ₁,μ₂,μ₃), 1, 3)
aΓ₂((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) = contribute_matrix(intrf₂, (u,λ₀,λ₁,λ₂,λ₃), (v,μ₀,μ₁,μ₂,μ₃), 1, 4)
aΓ₃((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) = contribute_matrix(intrf₃, (u,λ₀,λ₁,λ₂,λ₃), (v,μ₀,μ₁,μ₂,μ₃), 1, 5)

la₀((v,μ₀,μ₁,μ₂,μ₃)) = contribute_vector(intrf₀, (v,μ₀,μ₁,μ₂,μ₃), 2, ue_c₀)
la₁((v,μ₀,μ₁,μ₂,μ₃)) = contribute_vector(intrf₁, (v,μ₀,μ₁,μ₂,μ₃), 3, ue_c₁)
la₂((v,μ₀,μ₁,μ₂,μ₃)) = contribute_vector(intrf₂, (v,μ₀,μ₁,μ₂,μ₃), 4, ue_c₂)
la₃((v,μ₀,μ₁,μ₂,μ₃)) = contribute_vector(intrf₃, (v,μ₀,μ₁,μ₂,μ₃), 5, ue_c₃)

a((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) =  aΩ((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) + 
                                     aΓ₀((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) + 
                                     aΓ₁((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) + 
                                     aΓ₂((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃)) + 
                                     aΓ₃((u,λ₀,λ₁,λ₂,λ₃),(v,μ₀,μ₁,μ₂,μ₃))

l((v,μ₀,μ₁,μ₂,μ₃)) = ∫(v⋅f)*dΩ + 
                     la₀((v,μ₀,μ₁,μ₂,μ₃)) + 
                     la₁((v,μ₀,μ₁,μ₂,μ₃)) + 
                     la₂((v,μ₀,μ₁,μ₂,μ₃)) + 
                     la₃((v,μ₀,μ₁,μ₂,μ₃))

A = assemble_matrix(a,U,V)
b = assemble_vector(l,V)

##--------------------------------------------------

#op = AffineFEOperator(a,l,U,V)
op = AffineFEOperator(U,V,A,b)

ls = LUSolver()
solver = LinearFESolver(ls)

xh = solve(op);
uh, λh = xh;

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])