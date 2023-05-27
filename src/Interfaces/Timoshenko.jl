
struct Intrf_Timoshenko <: inter2D
  Γ::Triangulation
  Ω::Triangulation
  dΓ::Measure
  rot_cf::CellField
  Ef::CellField
  zf::CellField
  I::Float64
  L::Float64
  Da::Float64
  Db::Float64
  Dd::Float64
  Dinv::Float64
  Aa::Float64
  Ab::Float64
  Ainv::Float64
end
function Intrf_Timoshenko(Γ::Triangulation, Ω::Triangulation, degree::Int64, Ef::CellField, zf::CellField)
  dΓ = Measure(Γ,degree)
  
  Da_fun(Ef)    = sum(∫(       Ef )*dΓ)
  Db_fun(Ef,zf) = sum(∫(    zf*Ef )*dΓ)
  Dd_fun(Ef,zf) = sum(∫( zf*zf*Ef )*dΓ)
  S__fun(Ef)    = sum(∫(       Ef )*dΓ)
  L__fun = sum(∫(   1.0 )*dΓ)               
  I__fun = sum(∫( zf*zf )*dΓ)
  
  Da = Da_fun(Ef)
  Db = Db_fun(Ef,zf)
  Dd = Dd_fun(Ef,zf)
  Dinv = 1.0/(Dd*Da-Db^2)

  da(z_val) = sum(∫(    step_field(zf,z_val,Γ)*Ef )*dΓ)
  db(z_val) = sum(∫( zf*step_field(zf,z_val,Γ)*Ef )*dΓ)
  Aa = sum( ∫( da∘(zf) )*dΓ )
  Ab = sum( ∫( db∘(zf) )*dΓ )
  Ainv = 1.0/(Ab*Da-Aa*Db)

  L  = L__fun
  I  = I__fun
  
  rot_cf = get_rot_arr(Γ)

  return Intrf_Timoshenko(Γ,Ω,dΓ,rot_cf,Ef,zf,I,L,Da,Db,Dd,Dinv,Aa,Ab,Ainv)    
end

function contribute_matrix(intrf::Intrf_Timoshenko, U_basis, V_basis,
  U_ind::Int64, V_ind::Int64)
  get_i(i,x) = x[i]
  
  #function step(z::Float64,z_val::Float64)
  #  if z <= (z_val)
  #    return 1.0
  #  else
  #    return 0.0
  #  end
  #end
  
  u = U_basis[U_ind]; v = V_basis[U_ind]
  λ = U_basis[V_ind]; μ = V_basis[V_ind]
  
  #step_field(zf,z_val) = CellField(step.(zf,z_val),intrf.Γ)
  #step_field(zf,z_val) = CellField(step.(zf,z_val),intrf.Ω,PhysicalDomain())
  
  da_fun(Ef,zf,z_val) = sum(∫(    step_field(zf,z_val,intrf.Γ)*Ef )*intrf.dΓ)
  db_fun(Ef,zf,z_val) = sum(∫( zf*step_field(zf,z_val,intrf.Γ)*Ef )*intrf.dΓ)
  da(z_val) = da_fun(intrf.Ef,intrf.zf,z_val)
  db(z_val) = db_fun(intrf.Ef,intrf.zf,z_val)
  
  function comp_c_arr_cf(cₚ,cₘ,cᵥ)
    return TensorValue{3,2}(cₚ,cₘ,0.0,0.0,0.0,cᵥ) # [1,1] [2,1] [3,1] [1,2] ...
  end
  
  cₚ = (intrf.L*intrf.Dinv)*(   intrf.Dd*intrf.Ef      - intrf.Db*intrf.zf*intrf.Ef )
  cₘ = (intrf.L*intrf.Dinv)*( - intrf.Db*intrf.Ef      + intrf.Da*intrf.zf*intrf.Ef )
  cᵥ = (intrf.L*intrf.Ainv)*(   intrf.Db*(da∘intrf.zf) - intrf.Da*(db∘intrf.zf)     )
  
  c_arr = comp_c_arr_cf∘(cₚ,cₘ,cᵥ)
  
  return ∫( (λ⋅(c_arr⋅(intrf.rot_cf⋅v))) + (μ⋅(c_arr⋅(intrf.rot_cf⋅u))) )*intrf.dΓ
end

function contribute_vector(intrf::Intrf_Timoshenko, V_basis, V_ind::Int64, f)
  μ = V_basis[V_ind]
  return ∫( μ⋅f )*intrf.dΓ
end

function print_info(intrf::Intrf_Timoshenko)
  println("----------------------------")
  println(intrf.Da)
  println(intrf.Db)
  println(intrf.Dd)
  println(intrf.L )
  println(intrf.I )
  println(intrf.Aa)
  println(intrf.Ab)
  println("----------------------------")
end

