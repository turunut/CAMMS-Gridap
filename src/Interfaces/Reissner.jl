
struct Intrf_Reissner <: inter3D
  Γ_3D::Triangulation
  fix_axis::Int64
  pos_axis::Int64
  glue::Any
  c2f_faces::Any
end
function Intrf_Reissner(Γ_3D::Triangulation, int_coords::Vector{VectorValue{3, Int64}}, fix_axis::Int64, pos_axis::Int64, inv::Bool)
  glue, c2f_faces = create_interface(Γ_3D, int_coords, fix_axis, pos_axis, inv)
  return Intrf_Reissner(Γ_3D, fix_axis, pos_axis, glue, c2f_faces)
end

function contribute_matrix(intrf::Intrf_Reissner, U_basis, V_basis, U_ind::Int64, V_ind::Int64)
  function step(z::Float64,z_val::Float64)
    if z <= (z_val)
      return 1.0
    else
      return 0.0
    end
  end
  
  u = U_basis[U_ind]; v = V_basis[U_ind]
  λ = U_basis[V_ind]; μ = V_basis[V_ind]
  
  #step_field(zf,z_val) = CellField(step.(zf,z_val),intrf.Γ)
  step_field(zf,z_val) = CellField(step.(zf,z_val),intrf.Ω)
  
  da_fun(Ef,zf,z_val) = sum(∫(    step_field(zf,z_val)*Ef )*intrf.dΓ)
  db_fun(Ef,zf,z_val) = sum(∫( zf*step_field(zf,z_val)*Ef )*intrf.dΓ)
  da(z_val) = da_fun(intrf.Ef,intrf.zf,z_val)
  db(z_val) = db_fun(intrf.Ef,intrf.zf,z_val)
  
  function comp_c_arr_cf(cₚ,cₘ,cᵥ)
      return TensorValue{3,2}(cₚ,cₘ,0.0,0.0,0.0,cᵥ) # [1,1] [2,1] [3,1] [1,2] ...
  end

  cₚ = (intrf.L/intrf.Da)*intrf.Ef
  cₘ = (intrf.L/intrf.Dd)*intrf.zf*intrf.Ef
  cᵥ = (intrf.L/intrf.Dd)*(db∘intrf.zf)
  c_arr = comp_c_arr_cf∘(cₚ,cₘ,cᵥ)
  return ∫( (λ⋅(c_arr⋅v)) + (μ⋅(c_arr⋅u)) )*intrf.dΓ
end