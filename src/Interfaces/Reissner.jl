
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

  ext_TrialFESpace::FESpace

  Ψ::Triangulation
  dΨ::Measure
  da_cf::CellField
  db_cf::CellField

  l::Float64
  I::Float64
  h::Float64
  Da::TensorValue{3, 3, Float64, 9}
  Db::TensorValue{3, 3, Float64, 9}
  Dd::TensorValue{3, 3, Float64, 9}
  Aa::TensorValue{3, 3, Float64, 9}
  Ab::TensorValue{3, 3, Float64, 9}
  invD::TensorValue{6, 6, Float64, 36}
  invA::TensorValue{4, 4, Float64, 16}

  Intrf_Reissner() = new()
end
function Intrf_Reissner(Ω, Γ, Ψ, CTf_2D, int_coords::Vector{VectorValue{3, Int64}}, fix_axis::Int64, pos_axis::Int64, degree::Int64, inver::Bool)
  intrf = Intrf_Reissner()

  intrf.Ω = Ω
  z_coord(x) = x[end]; intrf.zf = CellField(z_coord,Ω)

  intrf.dΓ = Measure(Γ,degree)
  intrf.glue, intrf.c2f_faces = create_interface(Γ, int_coords, fix_axis, pos_axis, inver)

  cface_model = CartesianDiscreteModel((0,1,0,1),(length(intrf.c2f_faces),1))
  intrf.Γc  = Triangulation(cface_model)
  intrf.Γf  = Adaptivity.GluedTriangulation(Γ,intrf.Γc,intrf.glue)
  intrf.dΓc = Measure(intrf.Γc,degree)
  intrf.dΓf = Measure(intrf.Γf,degree)

  intrf.Λe = get_line_model_triangulation(intrf.c2f_faces)
  intrf.rot_arr = get_rot_arr(Γ)

  intrf.CTf_2D = CTf_2D

  # edge Ψ
  intrf.Ψ = Ψ
  intrf.dΨ = Measure(Ψ,degree)

  # tensors Da Db Dd i linversa de D
  Da_fun(Ef)    = sum(∫(       Ef )*intrf.dΨ)
  Db_fun(Ef,zf) = sum(∫(    zf*Ef )*intrf.dΨ)
  Dd_fun(Ef,zf) = sum(∫( zf*zf*Ef )*intrf.dΨ)
  S__fun(Ef)    = sum(∫(       Ef )*intrf.dΨ)
  L__fun = sum(∫( 1.0 )*intrf.dΨ);               intrf.h = L__fun
  I__fun = sum(∫( intrf.zf*intrf.zf )*intrf.dΨ); intrf.I = I__fun
  
  intrf.Da = Da_fun(intrf.CTf_2D)
  intrf.Db = Db_fun(intrf.CTf_2D,intrf.zf)
  intrf.Dd = Dd_fun(intrf.CTf_2D,intrf.zf)

  invD = zeros(6,6)
  invD[1:3,1:3] = get_array(intrf.Da)
  invD[4:6,1:3] = get_array(intrf.Db); invD[1:3,4:6] = get_array(intrf.Db)
  invD[4:6,4:6] = get_array(intrf.Dd)
  intrf.invD = inv(invD)

  # da i db cellfields
  function f_da(x); z_val = x[end]; return sum( ∫(          step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end
  function f_db(x); z_val = x[end]; return sum( ∫( intrf.zf*step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end
  
  intrf.da_cf = CellField(f_da,intrf.Ψ)
  intrf.db_cf = CellField(f_db,intrf.Ψ)
  
  intrf.Aa = sum( ∫( intrf.da_cf )*intrf.dΨ )
  intrf.Ab = sum( ∫( intrf.db_cf )*intrf.dΨ )


  invA = zeros(4,4)
  invA[1:2,1:2] = get_array(intrf.Aa)[1:2,1:2]
  invA[3:4,1:2] = get_array(intrf.Da)[1:2,1:2]; invA[1:2,3:4] = get_array(intrf.Ab)[1:2,1:2]
  invA[3:4,3:4] = get_array(intrf.Db)[1:2,1:2]
  intrf.invA = transpose(inv(invA))
  
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



  function f_da(x); z_val = x[end]; return sum( ∫(          step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end
  function f_db(x); z_val = x[end]; return sum( ∫( intrf.zf*step_field(intrf.zf,z_val,intrf.Ψ)*intrf.CTf_2D )intrf.dΨ ); end



  function _my_tensor(z,CT_2D)
    da_db_arr = zeros(4,2)

    α_arr = invD ⋅ ( TensorValue{6,3}( 1.0, 0.0, 0.0,   z, 0.0, 0.0,
                                       0.0, 1.0, 0.0, 0.0,   z, 0.0,
                                       0.0, 0.0, 1.0, 0.0, 0.0,   z) ⋅ CT_2D )
    
    #AAA = da(z)
    #BBB = db(z)

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