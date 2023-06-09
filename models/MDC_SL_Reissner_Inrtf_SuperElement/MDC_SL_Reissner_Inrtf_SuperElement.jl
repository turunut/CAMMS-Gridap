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

order  = 2
degree = 2*order

#--------------------------------------------------
# Definim el model volumetric

domain = (0,4,0,4,-0.5,0.5)
partition = (8,8,4)
model = CartesianDiscreteModel(domain,partition)

writevtk(model,"models/"*prblName*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"faceY0",[23])
add_tag_from_tags!(labels,"faceY1",[24])
add_tag_from_tags!(labels,"faceX0",[25])
add_tag_from_tags!(labels,"faceX1",[26])


#--------------------------------------------------
# Definim les propietats del model solid 

modlType = Solid()

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(modlType, CT1)

CT2 = CT_Isotrop(72000, 0.3)
ct2 = modModel.computeCT(modlType, CT2)


#--------------------------------------------------
# Definim la triangulacio del model solid

Ω  = Triangulation(model)
dΩ = Measure(Ω,degree)


#--------------------------------------------------
# Definim el Cell Field de CT del model solid

dimens  = 3
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1)

CTf = get_CT_CellField(modlType, CTs, tags, Ω)


#---------------------------
# Definim les propietats de cada layer per fer un model Shell equivalent

modlType_2D = PlaneStress()

ct1_2D = modModel.computeCT(modlType_2D, CT1, true)
ct2_2D = modModel.computeCT(modlType_2D, CT2, true)

CTs_2D = hcat(ct1_2D)

CTf_2D = get_CT_CellField(modlType_2D, CTs_2D, tags, Ω)


#--------------------------------------------------
# Acotem les triangulacions de les interficies

# Totes juntes
Γ₉ = BoundaryTriangulation(model,tags=["tag_23","tag_26","tag_24","tag_25"])

# Cada una de les cares
Γ₀ = BoundaryTriangulation(model,tags=["tag_23"])
Γ₁ = BoundaryTriangulation(model,tags=["tag_26"])
Γ₂ = BoundaryTriangulation(model,tags=["tag_24"])
Γ₃ = BoundaryTriangulation(model,tags=["tag_25"])

# Model lineal al llarg de lespessor
Ψ₀ = EdgeTriangulation(model,["tag_17"])
Ψ₁ = EdgeTriangulation(model,["tag_18"])
Ψ₂ = EdgeTriangulation(model,["tag_20"])
Ψ₃ = EdgeTriangulation(model,["tag_19"])


#--------------------------------------------------
# Definim els models 1x1 fine i coarse del conjunt de les interficies
# Creem una llista de coordenades en forma de integer dels CG dels elements de les cares per
# despres poder poder filtrarlos per la seva posicio 

Dc      = num_cell_dims(model)
topo    = get_grid_topology(model)
f2n_map = Geometry.get_faces(topo,Dc-1,0)
coords  = Geometry.get_vertex_coordinates(topo)
face_coords = map(Reindex(coords),f2n_map)

tol = 1.e-3
com_coords = lazy_map(N->sum(N)/length(N),face_coords) # center of mass
int_coords = map(N->VectorValue(Int(floor(N[1]/tol)),Int(floor(N[2]/tol)),Int(floor(N[3]/tol))),com_coords)

# Definim leix i posicio en aquest de cada una de les interficies

z_coord(x) = x[end]; zf = CellField(z_coord,Ω)

axis_id₀ = 2; face_pos₀ = minimum(lazy_map(c->c[axis_id₀],int_coords))

axis_id₁ = 1; face_pos₁ = maximum(lazy_map(c->c[axis_id₁],int_coords))

axis_id₂ = 2; face_pos₂ = maximum(lazy_map(c->c[axis_id₂],int_coords))

axis_id₃ = 1; face_pos₃ = minimum(lazy_map(c->c[axis_id₃],int_coords))

#glue₉ = create_interface_global([Γ₀,Γ₁,Γ₂,Γ₃],
#                                int_coords, 
#                                [axis_id₀, axis_id₁, axis_id₂,axis_id₃],
#                                [face_pos₀,face_pos₁,face_pos₂,face_pos₃],
#                                [false,false,true,true])

# Definim el model global 1x1 i la glue corresponent de totes les interficies

Γ = [Γ₀,Γ₁,Γ₂,Γ₃]
axis_id = [axis_id₀, axis_id₁, axis_id₂,axis_id₃]
axis_int_coord = [face_pos₀,face_pos₁,face_pos₂,face_pos₃]
inv_list = [false,false,true,true]


global global_n2o_faces = Vector{Int}(undef, 0)
global global_child_ids = Vector{Int}(undef, 0)
global partial_glues = AdaptivityGlue[]

global totalColumns = 0
global n_refs = 0 # Divisio vertical
for iintrf in eachindex(Γ)
  n2o_faces, child_ids, o2n_faces = get_comp_glue(Γ[iintrf], int_coords, axis_id[iintrf], axis_int_coord[iintrf], inv_list[iintrf])
  n_coarse = length(o2n_faces)
  global n_refs = length(o2n_faces[1])

  n2o_faces[3] .+= totalColumns
  rrules = Fill(RefinementRule(QUAD,(1,n_refs)),n_coarse)
  glue = AdaptivityGlue(n2o_faces,child_ids,rrules) # From coarse to fine
  push!(partial_glues, glue)

  global global_n2o_faces = vcat(global_n2o_faces, n2o_faces[3])
  global global_child_ids = vcat(global_child_ids, child_ids)
  global totalColumns += n_coarse
end

global_rrules = Fill(RefinementRule(QUAD,(1,n_refs)),totalColumns)
glue = AdaptivityGlue([Int[],Int[],global_n2o_faces],global_child_ids,global_rrules) # From coarse to fine

cface_model = CartesianDiscreteModel((0,1,0,1),(totalColumns,1),isperiodic=(true,false))
Γc  = Triangulation(cface_model)
Γf  = Gridap.Adaptivity.GluedTriangulation(Γ₉,Γc,glue)

# Definim 

Γ₀_glue = Gridap.Adaptivity.GluedTriangulation(Γ₀,Γc,partial_glues[1])
Γ₁_glue = Gridap.Adaptivity.GluedTriangulation(Γ₁,Γc,partial_glues[2])
Γ₂_glue = Gridap.Adaptivity.GluedTriangulation(Γ₂,Γc,partial_glues[3])
Γ₃_glue = Gridap.Adaptivity.GluedTriangulation(Γ₃,Γc,partial_glues[4])

intrf₀  = Intrf_ReissnerV2(Ω, Γf, Ψ₀, Γ₀_glue, CTf_2D[1], degree)
intrf₁  = Intrf_ReissnerV2(Ω, Γf, Ψ₁, Γ₁_glue, CTf_2D[1], degree)
intrf₂  = Intrf_ReissnerV2(Ω, Γf, Ψ₂, Γ₂_glue, CTf_2D[1], degree)
intrf₃  = Intrf_ReissnerV2(Ω, Γf, Ψ₃, Γ₃_glue, CTf_2D[1], degree)


#--------------------------------------------------
# Definim els FE del model volumetric i dels multiplicadors de lagrange (un per cada columna delements) 
# Definim el model 3D
reffe_u  = ReferenceFE(lagrangian,VectorValue{3,Float64},order)

Vu = TestFESpace(model,reffe_u;conformity=:H1)
Uu = TrialFESpace(Vu)

# Definim l'espai de l'acoplament
dofs  = 5
reffe = ReferenceFE(lagrangian,VectorValue{dofs,Float64},0)
Vλ = FESpace(Γc,reffe,conformity=:L2)
Uλ = TrialFESpace(Vλ)

VΨ₀ = ConstantFESpace(model,field_type=VectorValue{dofs,Float64}); UΨ₀ = TrialFESpace(VΨ₀)
VΨ₁ = ConstantFESpace(model,field_type=VectorValue{dofs,Float64}); UΨ₁ = TrialFESpace(VΨ₁)
VΨ₂ = ConstantFESpace(model,field_type=VectorValue{dofs,Float64}); UΨ₂ = TrialFESpace(VΨ₂)
VΨ₃ = ConstantFESpace(model,field_type=VectorValue{dofs,Float64}); UΨ₃ = TrialFESpace(VΨ₃)

# Ajuntem els dos espais anteriors
#V = MultiFieldFESpace([Vu,Vλ])
#U = MultiFieldFESpace([Uu,Uλ])
V = MultiFieldFESpace([Vu,Vλ,VΨ₀,VΨ₁,VΨ₂,VΨ₃])
U = MultiFieldFESpace([Uu,Uλ,UΨ₀,UΨ₁,UΨ₂,UΨ₃])


#--------------------------------------------------
# Model linea fine amb periodicitat als seus extrems

modelΛf = CartesianDiscreteModel((0,1),totalColumns;isperiodic=Tuple(true))
Λf = Triangulation(modelΛf)

dofs = 5
reffeΛ = ReferenceFE(lagrangian,VectorValue{dofs,Float64},order)
VΛf = FESpace(Λf,reffeΛ,conformity=:H1)
UΛf = TrialFESpace(VΛf)


#--------------------------------------------------
interfaces = [intrf₀, intrf₁, intrf₂, intrf₃]
crs_dsc = [1,1,1,1]

# Defineixo el model lineal coarse amb peridicitat als extrems
domainΛc = (0, 1); partitionΛc = (sum(crs_dsc))
modelΛc  = CartesianDiscreteModel(domainΛc, partitionΛc; isperiodic=Tuple(true))
Λc  = Triangulation(modelΛc)
# i el seu FEspace
VΛc = FESpace(modelΛc,reffeΛ,conformity=:H1); UΛc = TrialFESpace(VΛc)


#--------------------------------------------------
# Dic quin DOF vui activar
active_DOF = 3

# Projeccio de la funcio en el model linea al llarg de leix Y del model 2D coarse de les interfaces
_get_y(x) = VectorValue(x[2])
function π_Λe_Γc(f::CellField, Γc::Triangulation)
  _data = CellData.get_data(f)
  _cellmap = Fill(Broadcasting(_get_y),length(_data))
  data = lazy_map(∘,_data,_cellmap) # data = lazy_map(Broadcasting(∘),_data,_cellmap) 
  return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
end

# Defineixo la funcio activant el DOF
xC = zero_free_values(UΛc);
if active_DOF != 0; xC[active_DOF] = 1.0; end
fun_uC = FEFunction(UΛc,xC)

# Interpolo la funcio del model lineal coarse al model lineal fine
uC_intrp = Interpolable(fun_uC)
# Portem la funcio de la linea coarse al model fine
fun_ue = interpolate(uC_intrp,UΛf)
fun_ue(Point(0.0))

ue_c = π_Λe_Γc(fun_ue,Γc)



struct InterfaceProblem
  Ω::Triangulation
  Γf::Triangulation
  Γc::Triangulation
  Λ::Triangulation
  intrf :: Vector{Intrf_ReissnerV2}
end

prob = InterfaceProblem(Ω, Γf, Γc, Λf, [intrf₀, intrf₁, intrf₂, intrf₃])

function func(problem::InterfaceProblem,u,v,λ,μ)

  λ_f = change_domain(λ,problem.Γf,ReferenceDomain())
  μ_f = change_domain(μ,problem.Γf,ReferenceDomain())

  contr = DomainContribution()
  for intrf in problem.intrf
    dΓfi = intrf.dΓi

    invD = intrf.invD
    invA = intrf.invA
  
    function f_da(x); z_val = x[end]; return sum( ∫(          step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end
    function f_db(x); z_val = x[end]; return sum( ∫( intrf.zf*step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end
  
    function _my_tensor(z,CT_2D)
      da_db_arr = zeros(4,2)
  
      α_arr = invD ⋅ ( TensorValue{6,3}( 1.0, 0.0, 0.0,   z, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0,   z, 0.0,
                                         0.0, 0.0, 1.0, 0.0, 0.0,   z) ⋅ CT_2D )
  
      AAA = f_da(z)
      BBB = f_db(z)
      
      da_db_arr[1:2,1:2] = get_array( AAA )[1:2,1:2]
      da_db_arr[3:4,1:2] = get_array( BBB )[1:2,1:2]
      
      β_arr = invA ⋅ ( TensorValue{4,2}( da_db_arr ) )
      
      A_arr = zeros(5,3)
      A_arr[1:4,1:2] .= get_array( α_arr )[[1,3,4,6],[1,3]]
      A_arr[5,3]      = get_array( β_arr )[1,1]
      return TensorValue{5,3}(A_arr)
    end
  
    T = _my_tensor∘(intrf.zf, intrf.CTf_2D)

    contr += ∫(λ_f⋅T⋅v + μ_f⋅T⋅u)*dΓfi
  end
  return contr
end



function func_edge(intrf,u,v,α,β,dΨ)
  invD = intrf.invD
  invA = intrf.invA

  function f_da(x); z_val = x[end]; return sum( ∫(          step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end
  function f_db(x); z_val = x[end]; return sum( ∫( intrf.zf*step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end

  function _my_tensor(z,CT_2D)
    da_db_arr = zeros(4,2)

    α_arr = invD ⋅ ( TensorValue{6,3}( 1.0, 0.0, 0.0,   z, 0.0, 0.0,
                                       0.0, 1.0, 0.0, 0.0,   z, 0.0,
                                       0.0, 0.0, 1.0, 0.0, 0.0,   z) ⋅ CT_2D )

    AAA = f_da(z)
    BBB = f_db(z)
    
    da_db_arr[1:2,1:2] = get_array( AAA )[1:2,1:2]
    da_db_arr[3:4,1:2] = get_array( BBB )[1:2,1:2]
    
    β_arr = invA ⋅ ( TensorValue{4,2}( da_db_arr ) )
    
    A_arr = zeros(5,3)
    A_arr[1:4,1:2] .= get_array( α_arr )[[1,3,4,6],[1,3]]
    A_arr[5,3]      = get_array( β_arr )[1,1]
    return TensorValue{5,3}(A_arr)
  end

  T = _my_tensor∘(intrf.zf, intrf.CTf_2D)

  return ∫(α⋅T⋅v + β⋅T⋅u)*dΨ
end

dΨ₀ = Measure(Ψ₀,degree)
dΨ₁ = Measure(Ψ₁,degree)
dΨ₂ = Measure(Ψ₂,degree)
dΨ₃ = Measure(Ψ₃,degree)

aΩ((u,λ,α₀,α₁,α₂,α₃),(v,μ,β₀,β₁,β₂,β₃)) = ∫( ∂(v)⊙σ(CTf[1],∂(u)) )*dΩ
aΓ((u,λ,α₀,α₁,α₂,α₃),(v,μ,β₀,β₁,β₂,β₃)) = func(prob,u,v,λ,μ)
aΨ((u,λ,α₀,α₁,α₂,α₃),(v,μ,β₀,β₁,β₂,β₃)) = func_edge(intrf₀,u,v,α₀,β₀,dΨ₀) + 
                                          func_edge(intrf₁,u,v,α₁,β₁,dΨ₁) + 
                                          func_edge(intrf₂,u,v,α₂,β₂,dΨ₂) + 
                                          func_edge(intrf₃,u,v,α₃,β₃,dΨ₃)

a((u,λ,α₀,α₁,α₂,α₃),(v,μ,β₀,β₁,β₂,β₃)) =  aΩ((u,λ,α₀,α₁,α₂,α₃),(v,μ,β₀,β₁,β₂,β₃)) + 
                                          aΓ((u,λ,α₀,α₁,α₂,α₃),(v,μ,β₀,β₁,β₂,β₃)) + 
                                          aΨ((u,λ,α₀,α₁,α₂,α₃),(v,μ,β₀,β₁,β₂,β₃))

######fun_mult(u,v,λ,μ) = func(prob,u,v,λ,μ)
######UU = get_trial_fe_basis(U)
######VV = get_fe_basis(V)
######contrA = fun_mult(UU[1],VV[1],UU[2],VV[2])
######elementA = first(contrA.dict).second[1]

intrf = intrf₀
Γi = intrf.Γi

get_cell_dof_ids(Vλ,Γi)

u,λ = get_trial_fe_basis(U);
v,μ = get_fe_basis(U);


pts = get_cell_points(Γi)
λf = change_domain(λ,Γf,ReferenceDomain())
λfi = change_domain(λf,Γi,ReferenceDomain())
λfi(pts)

#ui = change_domain(u,prob.intrf[1].Γi,ReferenceDomain())

A = assemble_matrix(a,U,V)

#for i in 82:101
#  println( maximum(A[i,:]) )
#end

function func_edge_vector(β,g,dΨ)
  return ∫( β⋅g )*dΨ
end

f  = VectorValue(0.0,0.0,0.0)
g₀ = VectorValue(0.0,0.0,1.0,0.0,0.0)
g₁ = VectorValue(0.0,0.0,0.0,0.0,0.0)
g₂ = VectorValue(0.0,0.0,0.0,0.0,0.0)
g₃ = VectorValue(0.0,0.0,0.0,0.0,0.0)

dΓf = Measure(Γf, degree)
l((v,μ,β₀,β₁,β₂,β₃)) = ∫(v⋅f)*dΩ + 
                       ∫(μ⋅ue_c)*dΓf + 
                       #func_edge_vector(intrf₀,β₀,g₀,dΨ₀) + 
                       #func_edge_vector(intrf₁,β₁,g₁,dΨ₁) + 
                       #func_edge_vector(intrf₂,β₂,g₂,dΨ₂) + 
                       #func_edge_vector(intrf₃,β₃,g₃,dΨ₃)
                       func_edge_vector(β₀,g₀,dΨ₀) + 
                       func_edge_vector(β₁,g₁,dΨ₁) + 
                       func_edge_vector(β₂,g₂,dΨ₂) + 
                       func_edge_vector(β₃,g₃,dΨ₃)

b = assemble_vector(l,V)

b[(Uu.nfree+1):end]


#--------------------------------------------------
# Multiplicadors

#op = AffineFEOperator(a,l,U,V)
op = AffineFEOperator(U,V,A,b)

ls = LUSolver()
solver = LinearFESolver(ls)

xh = solve(op);
uh, λh = xh;

for i in 1:5:(5*partition[1])
  println(get_free_dof_values(λh)[i])
end

#x = get_free_dof_values(xh)
#xλ = Gridap.MultiField.restrict_to_field(U,x,2)

#λ₀ = xh.single_fe_functions[2]
#λ₁ = xh.single_fe_functions[3]
#λ₂ = xh.single_fe_functions[4]
#λ₃ = xh.single_fe_functions[5]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u"=>uh,
                     "ε"=>∂(uh),
                     "σ"=>σ(CTf[1],∂(uh))])
