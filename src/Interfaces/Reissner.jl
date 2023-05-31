
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
  invD::Float64

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

  db_fun(Ef,zf,z_val) = sum(∫( zf*step_field(zf,z_val,intrf.Ω)*Ef )*intrf.dΓ)
  db(z_val) = db_fun(intrf.CTf_2D,intrf.zf,z_val)

  chos = TensorValue{7,3}(1,0,0, 1,0,0, 0,
                          0,1,0, 0,1,0, 0,
                          0,0,0, 0,0,0, 1)
  
  function comp_c_arr_cf(c_vec)
      return TensorValue{7,1}(cₚ,cₘ,0.0,0.0,0.0,cᵥ) # [1,1] [2,1] [3,1] [1,2] ...
  end

  c_vec = intrf.l*( intrf.invD ⋅ intrf.CTf_2D )

  c_arr = comp_c_arr_cf∘(cₚ,cₘ,cᵥ)
  return ∫( (λ⋅(c_arr ⋅ (chos⋅v))) + (μ⋅(c_arr ⋅ (chos⋅u))) )*intrf.dΓf
end