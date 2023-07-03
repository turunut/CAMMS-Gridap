module modSubroutines
  export get_tags, get_E_CellField, ctˣ, step_field, EdgeTriangulation

  using Gridap.Helpers
  using Gridap.Geometry
  using Gridap.Arrays
  using Gridap.Fields
  using Gridap.CellData
  using FillArrays
  using modCT

  function get_tags(matFlag, labels, dimens)
    tags = copy(get_face_tag(labels,dimens))
    if length(matFlag) != 0
      for (iflag, flag) in enumerate(matFlag)
          num_tag = get_tag_from_name(labels,flag)
          replace!(tags, num_tag => iflag)
      end
    else
      tags .= 1
    end
    return tags
  end

  function get_E_CellField(listCTs, tags, Ω::Triangulation)
    listEs = Vector{Float64}()
    for ct in listCTs
      push!(listEs, getE(ct, [1.0, 0.0, 0.0]))
    end
    cts = CompressedArray(listEs, tags)
    return CellField(cts,Ω,PhysicalDomain())
  end

  ctˣ_ope(cf,p) = cf^p
  function ctˣ(cf::CellField,p::Integer)
    return Operation(ctˣ_ope)(cf,p)
  end
  
  function step(z::Float64,z_val::Float64)
    if z <= (z_val)
      return 1.0
    else
      return 0.0
    end
  end
  
  function step_field(zf::CellField,z_val::Float64,Ω::Triangulation)
    return CellField(step.(zf,z_val),Ω,PhysicalDomain())
  end

  function EdgeTriangulation(model::DiscreteModel,tags)
    D = Geometry.num_cell_dims(model)
    labeling = get_face_labeling(model)
    face_to_mask = get_face_mask(labeling,tags,D-2)
    face_to_bgface = findall(face_to_mask)
    return EdgeTriangulation(model,face_to_bgface)
  end

  function EdgeTriangulation(
    model::DiscreteModel,
    face_to_bgface::AbstractVector{<:Integer})

    D = Geometry.num_cell_dims(model)
    topo = get_grid_topology(model)
    bgface_grid = Grid(Geometry.ReferenceFE{D-2},model)
    bgface_to_lcell = Fill(1,Geometry.num_faces(model,D-2))

    face_grid = view(bgface_grid,face_to_bgface)
    cell_grid = get_grid(model)
    glue = Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
    trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)

    BoundaryTriangulation(trian,glue)
  end

  function EdgeTriangulation(Γ::Triangulation{Dc,Dp},tags) where {Dc,Dp}
    model = get_background_model(Γ)
    labeling = get_face_labeling(model)
    face_to_mask = get_face_mask(labeling,tags,Dc-1)
    face_to_bgface = findall(face_to_mask)
    return EdgeTriangulation(Γ,face_to_bgface)
  end

  function EdgeTriangulation(
    Γ::Triangulation{Dc,Dp},
    face_to_bgface::AbstractVector{<:Integer}) where {Dc,Dp}
    @check Dc == Dp-1
    
    Λ  = SkeletonTriangulation(Γ)
    glue = get_glue(Λ,Val(Dc)) # glue from the edges to the attached faces
    
    model = get_background_model(Γ)
    topo = get_grid_topology(model)
    edge_to_face_map = Geometry.get_faces(topo,Dc-1,Dc) # edge to face connectivity
    
    # Find local numbering of the edges we want
    face_pairs = map(p -> [p...],zip(glue.plus.tface_to_mface,glue.minus.tface_to_mface))
    local_edge_ids = map(face_to_bgface) do edge
      edge_faces = edge_to_face_map[edge]
      return findfirst(faces -> faces == edge_faces, face_pairs)
    end
    return view(Λ,local_edge_ids)
  end

end