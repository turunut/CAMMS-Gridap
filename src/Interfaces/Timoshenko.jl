
mutable struct Intrf_Timoshenko <: inter2D
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

  da_cf::CellField
  db_cf::CellField

  Aa::Float64
  Ab::Float64
  Ainv::Float64

  Intrf_Timoshenko() = new()
end
function Intrf_Timoshenko(Γ::Triangulation, Ω::Triangulation, degree::Int64, Ef::CellField, zf::CellField)
  intrf = Intrf_Timoshenko()

  intrf.Γ = Γ
  intrf.Ω = Ω
  intrf.Ef = Ef
  intrf.zf = zf

  intrf.dΓ = Measure(Γ,degree)
  
  Da_fun(Ef)    = sum(∫(       Ef )*intrf.dΓ)
  Db_fun(Ef,zf) = sum(∫(    zf*Ef )*intrf.dΓ)
  Dd_fun(Ef,zf) = sum(∫( zf*zf*Ef )*intrf.dΓ)
  S__fun(Ef)    = sum(∫(       Ef )*intrf.dΓ)
  L__fun = sum(∫(   1.0 )*intrf.dΓ)               
  I__fun = sum(∫( zf*zf )*intrf.dΓ)
  
  intrf.Da = Da_fun(Ef)
  intrf.Db = Db_fun(Ef,zf)
  intrf.Dd = Dd_fun(Ef,zf)
  intrf.Dinv = 1.0/(intrf.Dd*intrf.Da-intrf.Db^2)

  function f_da(x); z_val = x[end]; return sum( ∫(    step_field(zf,z_val,Γ)*Ef )intrf.dΓ ); end
  function f_db(x); z_val = x[end]; return sum( ∫( zf*step_field(zf,z_val,Γ)*Ef )intrf.dΓ ); end

  intrf.da_cf = CellField(f_da,Γ)
  intrf.db_cf = CellField(f_db,Γ)

  # Dona el mateix  
  Aa_v2 = sum( ∫( intrf.da_cf )*intrf.dΓ )#; println(Aa_v2)
  Ab_v2 = sum( ∫( intrf.db_cf )*intrf.dΓ )#; println(Ab_v2)

  da(z_val) = sum(∫(    step_field(zf,z_val,Γ)*Ef )*intrf.dΓ)
  db(z_val) = sum(∫( zf*step_field(zf,z_val,Γ)*Ef )*intrf.dΓ)
  intrf.Aa = sum( ∫( da∘(zf) )*intrf.dΓ )
  intrf.Ab = sum( ∫( db∘(zf) )*intrf.dΓ )
  intrf.Ainv = 1.0/(intrf.Ab*intrf.Da-intrf.Aa*intrf.Db)

  intrf.L = L__fun
  intrf.I = I__fun
  
  intrf.rot_cf = get_rot_arr(Γ)

  return intrf
end

function contribute_matrix(intrf::Intrf_Timoshenko, U_basis, V_basis,
                                                    U_ind::Int64, V_ind::Int64)
  get_i(i,x) = x[i]
  
  u = U_basis[U_ind]; v = V_basis[U_ind]
  λ = U_basis[V_ind]; μ = V_basis[V_ind]
  
  da_fun(Ef,zf,z_val) = sum(∫(    step_field(zf,z_val,intrf.Γ)*Ef )*intrf.dΓ)
  db_fun(Ef,zf,z_val) = sum(∫( zf*step_field(zf,z_val,intrf.Γ)*Ef )*intrf.dΓ)
  da(z_val) = da_fun(intrf.Ef,intrf.zf,z_val)
  db(z_val) = db_fun(intrf.Ef,intrf.zf,z_val)
  
  function comp_c_arr_cf(cₚ,cₘ,cᵥ)
    return TensorValue{3,2}(cₚ,cₘ,0.0,0.0,0.0,cᵥ) # [1,1] [2,1] [3,1] [1,2] ...
  end
  
  cₚ = (intrf.L*intrf.Dinv)*(   intrf.Dd*intrf.Ef      - intrf.Db*intrf.zf*intrf.Ef )
  cₘ = (intrf.L*intrf.Dinv)*( - intrf.Db*intrf.Ef      + intrf.Da*intrf.zf*intrf.Ef )
  #cᵥ = (intrf.L*intrf.Ainv)*(   intrf.Db*(da∘intrf.zf) - intrf.Da*(db∘intrf.zf)     )
  cᵥ = (intrf.L*intrf.Ainv)*(   intrf.Db*(intrf.da_cf) - intrf.Da*(intrf.db_cf)     )
  
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

