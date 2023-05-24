module modInterface
  export create_interface, define_corse_fine_Γ, get_line_model
  export get_dofs, get_test_trial_spaces, contribute_matrix, contribute_vector, print_info
  export get_line_model_triangulation
  export Intrf_Timoshenko, Intrf_Kinematic2D
  export Intrf_Reissner,   Intrf_Kinematic3D

  using Gridap
  using Gridap.Arrays
  using Gridap.TensorValues
  using Gridap.Adaptivity
  using Gridap.CellData
  using Gridap.MultiField
  using FillArrays

  abstract type interface end

  abstract type inter2D <: interface end
  abstract type inter3D <: interface end

  get_i(x,i::Int64) = x[i]

  include("Interfaces/Kinematic2D.jl")
  include("Interfaces/Kinematic3D.jl")
  
  struct Intrf_Timoshenko <: inter2D
    Γ::Triangulation
    Ω::Triangulation
    dΓ::Measure
    Ef::CellField
    zf::CellField
    I::Float64
    L::Float64
    Da::Float64
    Db::Float64
    Dd::Float64
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
    L  = L__fun
    I  = I__fun

    return Intrf_Timoshenko(Γ,Ω,dΓ,Ef,zf,I,L,Da,Db,Dd)    
  end

  struct Intrf_Reissner <: inter3D
    Γ_3D::Triangulation
    fix_axis::Int64
    pos_axis::Int64
    glue::Any
    c2f_faces::Table
  end
  function Intrf_Reissner(Γ_3D::Triangulation, int_coords::Vector{VectorValue{3, Int64}}, fix_axis::Int64, pos_axis::Int64)
    glue, c2f_faces = create_interface(Γ_3D, int_coords, fix_axis, pos_axis)
    return Intrf_Reissner(Γ_3D, fix_axis, pos_axis, glue, c2f_faces)
  end

  get_dofs(intrf::Intrf_Kinematic2D) = 2
  get_dofs(intrf::Intrf_Kinematic3D) = 3
  get_dofs(intrf::Intrf_Timoshenko)  = 3
  get_dofs(intrf::Intrf_Reissner)    = 5

  function get_test_trial_spaces(intrf::interface, model::DiscreteModel)
    dofs = get_dofs(intrf)
    Vλ = ConstantFESpace(model,field_type=VectorValue{dofs,Float64})
    Uλ = TrialFESpace(Vλ)
    return Vλ, Uλ
  end

  function contribute_matrix(intrf::Intrf_Timoshenko, U_basis, V_basis,
                                                      U_ind::Int64, V_ind::Int64)
    get_i(i,x) = x[i]
    
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

  function contribute_matrix(intrf::Intrf_Reissner, U_basis, V_basis,
                                                    U_ind::Int64, V_ind::Int64)
    
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

  function print_info(intrf::Intrf_Timoshenko)
    println("----------------------------")
    println(intrf.Da)
    println(intrf.Db)
    println(intrf.Dd)
    println(intrf.L )
    println(intrf.I )
    println("----------------------------")
  end

  function in_matrix(intrf::Intrf_Timoshenko, aΓ)
    aΓ = ∫( (λ⋅v) + (μ⋅u) )*intrf.dΓ
    return aΓ
  end

  function get_line_model(intrf::interface)
    line_model = CartesianDiscreteModel((0,1),(length(intrf.c2f_faces)))
    Λe = Triangulation(line_model)
    return line_model, Λe
  end

  function create_interface(Γ::Triangulation, int_coords::Vector{VectorValue{3, Int64}}, axis_id::Int64, axis_int_coord::Int64)
    axis_p = [2,1][axis_id] # retorna laxis perpendicular al entrat
    println(axis_p)
    face_x = axis_int_coord
    interface_faces = findall(c->c[axis_id]==face_x,int_coords)
    interface_coords = view(int_coords,interface_faces)
    
    y_coords = map(c->c[axis_p],interface_coords)
    #y_unique = sort(unique(map(c->c[axis_p],interface_coords)))
    y_counts = [count(==(y),y_coords) for y in sort(unique(y_coords))]
    y_ptrs = Gridap.Adaptivity.counts_to_ptrs(y_counts)
    
    perm = sortperm(interface_coords,by=x->x[axis_p])
    data = lazy_map(Reindex(interface_faces),perm)
    c2f_faces = Table(data,y_ptrs)
    
    n2o_cells = zeros(Int, length(Γ.glue.face_to_bgface))
    child_ids = zeros(Int, length(Γ.glue.face_to_bgface))
    for (islide, slide_list) in enumerate(c2f_faces)
      for (izpos, num) in enumerate(slide_list)
        ind = findall(x->x==num, Γ.glue.face_to_bgface)[1]
        n2o_cells[ind] = islide
        child_ids[ind] = izpos
      end  
    end
    
    n2o_faces = [Int[],Int[],n2o_cells]
    rrules = Fill(RefinementRule(QUAD,(1,length(c2f_faces[1]))),length(c2f_faces))
    glue = AdaptivityGlue(n2o_faces,child_ids,rrules) # From coarse to fine
    return glue, c2f_faces
  end

end