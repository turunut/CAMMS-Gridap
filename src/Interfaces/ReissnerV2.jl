
mutable struct Intrf_ReissnerV2 <: inter3D
  Ω::Triangulation

  Γ::Triangulation
  dΓ::Measure

  Γi::Triangulation
  dΓi::Measure
  
  rot_arr::CellField

  c2f_faces::Any

  CTf_2D::CellField
  zf::CellField

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

  Intrf_ReissnerV2() = new()
end

function Intrf_ReissnerV2(Ω, Γf, Ψ, Γi, CTf_2D, degree::Int64)
  intrf = Intrf_ReissnerV2()

  intrf.Ω = Ω
  z_coord(x) = x[end]; intrf.zf = CellField(z_coord,Ω)

  intrf.Γ = Γf
  intrf.dΓ = Measure(Γf,degree)

  intrf.Γi = Γi
  intrf.dΓi = Measure(Γi,degree)

  intrf.rot_arr = get_rot_arr(Γi)

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
  invA[1:2,3:4] = get_array(intrf.Ab)[1:2,1:2]
  invA[3:4,1:2] = get_array(intrf.Da)[1:2,1:2]
  invA[3:4,3:4] = get_array(intrf.Db)[1:2,1:2]
  intrf.invA = transpose(inv(invA))
  
  return intrf
end


