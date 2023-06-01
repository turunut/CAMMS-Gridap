
struct Intrf_Kinematic3D <: inter3D
  Γ::Triangulation
  dΓ::Measure

  Γc::Triangulation
  dΓc::Measure

  Γf::Triangulation
  dΓf::Measure
  
  Λe::Triangulation # line model

  rot_arr::CellField

  fix_axis::Int64
  pos_axis::Int64
  glue::Any
  c2f_faces::Any
end
function Intrf_Kinematic3D(Γ, int_coords::Vector{VectorValue{3, Int64}}, fix_axis::Int64, pos_axis::Int64, degree::Int64, inv::Bool)
  dΓ = Measure(Γ,degree)
  glue, c2f_faces = create_interface(Γ, int_coords, fix_axis, pos_axis, inv)

  cface_model = CartesianDiscreteModel((0,1,0,1),(length(c2f_faces),1))
  Γc  = Triangulation(cface_model)
  Γf  = Adaptivity.GluedTriangulation(Γ,Γc,glue)
  dΓc = Measure(Γc,degree)
  dΓf = Measure(Γf,degree)

  Λe = get_line_model_triangulation(c2f_faces)
  rot_arr = get_rot_arr(Γ)

  return Intrf_Kinematic3D(Γ, dΓ, Γc, dΓc, Γf, dΓf, Λe, rot_arr, fix_axis, pos_axis, glue, c2f_faces)
end

function contribute_matrix(intrf::Intrf_Kinematic3D, U_basis, V_basis,
                                                     U_ind::Int64, V_ind::Int64)
  u = U_basis[U_ind]; v = V_basis[U_ind]
  λ = U_basis[V_ind]; μ = V_basis[V_ind]
  tr_Γf(λ) = change_domain(λ,intrf.Γf,DomainStyle(λ))
  return ∫( (tr_Γf(λ)⋅v) + (tr_Γf(μ)⋅u) )*intrf.dΓf
end

function get_rot_arr(Γ::Triangulation{2})
  n̂ = get_normal_vector(Γ)
  return CellField(get_arr_g2l_3D∘(n̂),Γ)
end
function get_arr_g2l_3D(n̂::VectorValue{3, Float64})
  n̂ *= sign(n̂[1])
  arr_g2l = TensorValue{3,3}( n̂[1], -n̂[2], 0, n̂[2], n̂[1], 0,     0,0,1 )
  return arr_g2l
end
