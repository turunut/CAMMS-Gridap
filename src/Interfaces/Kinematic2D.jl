
struct Intrf_Kinematic2D <: inter2D
  Γ::Triangulation
  dΓ::Measure
end
function Intrf_Kinematic2D(Γ::Triangulation, degree::Int64)
  dΓ = Measure(Γ,degree)
  return Intrf_Kinematic2D(Γ,dΓ)    
end

function contribute_matrix(intrf::Intrf_Kinematic2D, U_basis, V_basis,
                                                     U_ind::Int64, V_ind::Int64)
  u = U_basis[U_ind]; v = V_basis[U_ind]
  λ = U_basis[V_ind]; μ = V_basis[V_ind]
  return ∫( (λ⋅v) + (μ⋅u) )*intrf.dΓ
end

function contribute_vector(intrf::Intrf_Kinematic2D, V_basis, V_ind::Int64, f)
  μ = V_basis[V_ind]
  return ∫( μ⋅f )*intrf.dΓ
end