
 mutable struct Intrf_Kinematic3DV2 <: inter3D
  Ω::Triangulation

  Γ::Triangulation
  dΓ::Measure

  Γi::Triangulation
  dΓi::Measure

  rot_arr::CellField

  c2f_faces::Any

  CTf_2D::CellField
  zf::CellField

  Intrf_Kinematic3DV2() = new()
end
function Intrf_Kinematic3DV2(Ω, Γf, Ψ, Γi, CTf_2D, degree::Int64)
  intrf = Intrf_Kinematic3DV2()

  intrf.Ω = Ω
  z_coord(x) = x[end]; intrf.zf = CellField(z_coord,Ω)

  intrf.Γ = Γf
  intrf.dΓ = Measure(Γf,degree)

  intrf.Γi = Γi
  intrf.dΓi = Measure(Γi,degree)

  intrf.rot_arr = get_rot_arr(Γi)

  intrf.CTf_2D = CTf_2D

  return intrf
end


