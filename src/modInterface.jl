module modInterface
  export create_interface, define_corse_fine_Γ
  export McCune

  using Gridap
  using Gridap.Arrays
  using Gridap.TensorValues
  using Gridap.Adaptivity
  using FillArrays

  abstract type interface end

  abstract type inter2D <: interface end
  abstract type inter3D <: interface end

  #mutable struct McCune <: interface
  struct McCune <: inter3D
    Γ_3D::Triangulation
    fix_axis::Int64
    pos_axis::Int64
    glue::Any
    c2f_faces::Table
    #McCune() = new()
  end

  function McCune(Γ_3D::Triangulation, int_coords::Vector{VectorValue{3, Int64}}, fix_axis::Int64, pos_axis::Int64)
    glue, c2f_faces = create_interface(Γ_3D, int_coords, fix_axis, pos_axis)
    return McCune(Γ_3D, fix_axis, pos_axis, glue, c2f_faces)
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

  function define_corse_fine_Γ(Γ::Triangulation, int_coords::Vector{VectorValue{3, Int64}}, axis_id::Int64, axis_int_coord::Int64)
    glue, c2f_faces = create_interface(Γ, int_coords, axis_id, axis_int_coord)
  
    cface_model = CartesianDiscreteModel((0,1,0,1),(length(c2f_faces),1))
    Γc = Triangulation(cface_model)
    Γf = Adaptivity.GluedTriangulation(Γ,Γc,glue)
    return c2f_faces, cface_model, Γc, Γf
  end

end