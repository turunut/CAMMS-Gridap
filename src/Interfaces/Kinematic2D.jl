
struct Intrf_Kinematic2D <: inter2D
  Γ::Triangulation
  dΓ::Measure
  rot_cf::CellField
end
function Intrf_Kinematic2D(Γ::Triangulation, degree::Int64)
  dΓ = Measure(Γ,degree)
  rot_cf = get_rot_arr(Γ)
  return Intrf_Kinematic2D(Γ,dΓ,rot_cf)    
end

function contribute_matrix(intrf::Intrf_Kinematic2D, U_basis, V_basis,
                                                     U_ind::Int64, V_ind::Int64)
  u = U_basis[U_ind]; v = V_basis[U_ind]
  λ = U_basis[V_ind]; μ = V_basis[V_ind]
  return ∫( (λ⋅(intrf.rot_cf⋅v)) + (μ⋅(intrf.rot_cf⋅u)) )*intrf.dΓ
end

function contribute_vector(intrf::Intrf_Kinematic2D, V_basis, V_ind::Int64, f)
  μ = V_basis[V_ind]
  return ∫( μ⋅(intrf.rot_cf⋅f) )*intrf.dΓ
end

function get_arr_g2l_2D(n̂::CellField)
  arr_g2l = TensorValue{2,2}( n̂[1], -n̂[2], n̂[2], n̂[1] ) * sign(n̂[1])
  return arr_g2l
end

function get_rot_arr(Γ::Triangulation{1})
    n̂ = get_normal_vector(Γ)
    return CellField(get_arr_g2l_2D∘(n̂),Γ)
end

