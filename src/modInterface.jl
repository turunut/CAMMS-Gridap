module modInterface
  export create_interface, define_corse_fine_Γ
  export get_dofs, get_test_trial_spaces, contribute_matrix, contribute_vector, print_info
  export get_line_model_triangulation, get_line_test_trial_spaces, get_line_distribution
  export Intrf_Timoshenko, Intrf_Kinematic2D
  export Intrf_Reissner,   Intrf_Kinematic3D

  export defineExternalLine
  export interpolate_nodal_displ

  using Gridap
  using Gridap.Arrays
  using Gridap.TensorValues
  using Gridap.Adaptivity
  using Gridap.CellData
  using Gridap.MultiField
  using Gridap.FESpaces
  using modSubroutines
  using FillArrays

  abstract type interface end

  abstract type inter2D <: interface end
  abstract type inter3D <: interface end

  get_i(x,i::Int64) = x[i]

  include("Interfaces/Inter3D.jl")
  include("Interfaces/Kinematic2D.jl")
  include("Interfaces/Kinematic3D.jl")
  include("Interfaces/Timoshenko.jl")
  include("Interfaces/Reissner.jl")

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

  #function get_line_model(intrf::interface)
  #  line_model = CartesianDiscreteModel((0,1),(length(intrf.c2f_faces)))
  #  Λe = Triangulation(line_model)
  #  return line_model, Λe
  #end
  
  function get_line_model_triangulation(c2f_faces)
    line_model = CartesianDiscreteModel((0,1),(length(c2f_faces)))
    Λe = Triangulation(line_model)
    return Λe
  end

  function create_interface(Γ::Triangulation, int_coords::Vector{VectorValue{3, Int64}}, axis_id::Int64, axis_int_coord::Int64, inv::Bool)
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

    if inv
      c2f_faces = reverse(c2f_faces)
    end
    
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
