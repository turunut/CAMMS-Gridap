
struct Intrf_Kinematic3D <: inter3D
  Γ::Triangulation
  dΓ::Measure

  Γc::Triangulation
  dΓc::Measure

  Γf::Triangulation
  dΓf::Measure
  
  Λe::Triangulation # line model

  fix_axis::Int64
  pos_axis::Int64
  glue::Any
  c2f_faces::Table
end
function Intrf_Kinematic3D(Γ, int_coords::Vector{VectorValue{3, Int64}}, fix_axis::Int64, pos_axis::Int64, degree::Int64)
  dΓ = Measure(Γ,degree)
  glue, c2f_faces = create_interface(Γ, int_coords, fix_axis, pos_axis)

  cface_model = CartesianDiscreteModel((0,1,0,1),(length(c2f_faces),1))
  Γc  = Triangulation(cface_model)
  Γf  = Adaptivity.GluedTriangulation(Γ,Γc,glue)
  dΓc = Measure(Γc,degree)
  dΓf = Measure(Γf,degree)

  Λe = get_line_model_triangulation(c2f_faces)

  return Intrf_Kinematic3D(Γ, dΓ, Γc, dΓc, Γf, dΓf, Λe, fix_axis, pos_axis, glue, c2f_faces)
end

function define_corse_fine_triangulation(Γ::Triangulation, int_coords::Vector{VectorValue{3, Int64}}, axis_id::Int64, axis_int_coord::Int64)
  glue, c2f_faces = create_interface(Γ, int_coords, axis_id, axis_int_coord)
  cface_model = CartesianDiscreteModel((0,1,0,1),(length(c2f_faces),1))
  Γc  = Triangulation(cface_model)
  Γf  = Adaptivity.GluedTriangulation(Γ,Γc,glue)
  return c2f_faces, cface_model, Γc, Γf
end

function get_line_model_triangulation(c2f_faces::Table)
  line_model = CartesianDiscreteModel((0,1),(length(c2f_faces)))
  Λe = Triangulation(line_model)
  return Λe
end

function get_test_trial_spaces(intrf::Intrf_Kinematic3D)
  dofs = get_dofs(intrf)
  reffe = ReferenceFE(lagrangian,VectorValue{dofs,Float64},0)
  Vλ = FESpace(intrf.Γc,reffe,conformity=:L2)
  Uλ = TrialFESpace(Vλ)
  return Vλ, Uλ
end

function get_line_test_trial_spaces(intrf::Intrf_Kinematic3D, order)
  dofs = get_dofs(intrf)
  reffe = ReferenceFE(lagrangian,VectorValue{dofs,Float64},order)
  Vλ = FESpace(intrf.Λe,reffe,conformity=:H1)
  Uλ = TrialFESpace(Vλ)
  return Vλ, Uλ
end

function contribute_matrix(intrf::Intrf_Kinematic3D, U_basis, V_basis,
                                                     U_ind::Int64, V_ind::Int64)
  u = U_basis[U_ind]; v = V_basis[U_ind]
  λ = U_basis[V_ind]; μ = V_basis[V_ind]
  tr_Γf(λ) = change_domain(λ,intrf.Γf,DomainStyle(λ))
  return ∫( (tr_Γf(λ)⋅v) + (tr_Γf(μ)⋅u) )*intrf.dΓf
end

function contribute_vector(intrf::Intrf_Kinematic3D, V_basis, V_ind::Int64, f)
  μ = V_basis[V_ind]
  tr_Γf(λ) = change_domain(λ,intrf.Γf,DomainStyle(λ))
  return ∫( tr_Γf(μ)⋅f )*intrf.dΓf
end