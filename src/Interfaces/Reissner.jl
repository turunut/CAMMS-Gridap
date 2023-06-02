
mutable struct Intrf_Reissner <: inter3D
  Ω::Triangulation

  Γ::Triangulation
  dΓ::Measure

  Γc::Triangulation
  dΓc::Measure

  Γf::Triangulation
  dΓf::Measure
  
  Λe::Triangulation

  rot_arr::CellField

  fix_axis::Int64
  pos_axis::Int64
  glue::Any
  c2f_faces::Any

  
  CTf_2D::CellField
  zf::CellField

  l::Float64
  I::Float64
  h::Float64
  Da::Float64
  Db::Float64
  Dd::Float64
  Aa::Float64
  Ab::Float64
  invD::TensorValue{6, 6, Float64, 36}
  invA::TensorValue{4, 4, Float64, 16}

  Intrf_Reissner() = new()
end
function Intrf_Reissner(Ω, Γ, int_coords::Vector{VectorValue{3, Int64}}, fix_axis::Int64, pos_axis::Int64, degree::Int64, inv::Bool)
  intrf = Intrf_Reissner()

  intrf.Ω = Ω

  intrf.dΓ = Measure(Γ,degree)
  intrf.glue, intrf.c2f_faces = create_interface(Γ, int_coords, fix_axis, pos_axis, inv)

  cface_model = CartesianDiscreteModel((0,1,0,1),(length(intrf.c2f_faces),1))
  intrf.Γc  = Triangulation(cface_model)
  intrf.Γf  = Adaptivity.GluedTriangulation(Γ,intrf.Γc,intrf.glue)
  intrf.dΓc = Measure(intrf.Γc,degree)
  intrf.dΓf = Measure(intrf.Γf,degree)

  intrf.Λe = get_line_model_triangulation(intrf.c2f_faces)
  intrf.rot_arr = get_rot_arr(Γ)
  
  return intrf
end

function contribute_matrix(intrf::Intrf_Reissner, U_basis, V_basis, 
                                                  U_ind::Int64, V_ind::Int64)
  u = U_basis[U_ind]; v = V_basis[U_ind]
  λ = U_basis[V_ind]; μ = V_basis[V_ind]

  invD = intrf.invD
  invA = intrf.invA

  tr_Γf(λ) = change_domain(λ,intrf.Γf,DomainStyle(λ))

  da_fun(CT,zf,z_val) = sum(∫(    step_field(zf,z_val,intrf.Γf)*CT )*intrf.dΓf)
  db_fun(CT,zf,z_val) = sum(∫( zf*step_field(zf,z_val,intrf.Γf)*CT )*intrf.dΓf)
  da(z_val) = da_fun(intrf.CTf_2D,intrf.zf,z_val)
  db(z_val) = db_fun(intrf.CTf_2D,intrf.zf,z_val)

  function _my_tensor(z,CT_2D)
    da_db_arr = zeros(4,2)

    α_arr = invD ⋅ ( TensorValue{6,3}( 1.0, 0.0, 0.0,   z, 0.0, 0.0,
                                       0.0, 1.0, 0.0, 0.0,   z, 0.0,
                                       0.0, 0.0, 1.0, 0.0, 0.0,   z) ⋅ CT_2D )
    
    AAA = da(z)
    BBB = db(z)
    println(AAA)
    da_db_arr[1:2,1:2] = get_array( AAA )[1:2,1:2]
    da_db_arr[3:4,1:2] = get_array( BBB )[1:2,1:2]
    
    β_arr = invA ⋅ ( TensorValue{4,2}( da_db_arr ) )
    
    A_arr = zeros(5,3)
    A_arr[1:4,1:2] .= get_array( α_arr )[[1,3,4,6],[1,3]]
    A_arr[5,3]     = get_array( β_arr )[1,1]
    return TensorValue{5,3}(A_arr)
  end

  T = _my_tensor∘(intrf.zf, intrf.CTf_2D)

  return ∫(tr_Γf(λ)⋅T⋅v)*intrf.dΓf + 
         ∫(tr_Γf(μ)⋅T⋅u)*intrf.dΓf
end

function computeWork(u, λ, invD, CT_2D, z)
  uᵣ =  VectorValue(u[1][1], 0.0, u[2][1], u[3][1])
  λᵣ =  VectorValue(λ[1][1], 0.0, λ[2][1], λ[3][1], 0.0, λ[4][1], λ[5][1])

  A_arr = zeros(7,4)

  α_arr = invD ⋅ ( TensorValue{6,3}( 1.0, 0.0, 0.0, z, 0.0, 0.0,
                                     0.0, 1.0, 0.0, 0.0, z, 0.0,
                                     0.0, 0.0, 1.0, 0.0, 0.0, z ) ⋅ CT_2D )

  A_arr[1:6,1:3] = get_array( α_arr )
  A_arr[7,4] = 1.0
                                   
  return λᵣ⋅(TensorValue{7,4}(A_arr)⋅uᵣ)
end

function _my_tensor(z,CT_2D,invD)
  α_arr = invD ⋅ ( TensorValue{6,3}( 1.0, 0.0, 0.0, z, 0.0, 0.0,
                                     0.0, 1.0, 0.0, 0.0, z, 0.0,
                                     0.0, 0.0, 1.0, 0.0, 0.0, z ) ⋅ CT_2D )
  A_arr = zeros(5,3)
  A_arr[1:4,1:2] .= get_array( α_arr )[[1,3,4,6],[1,3]]
  A_arr[5,3] = 1.0
  T = TensorValue{5,3}(A_arr)
end