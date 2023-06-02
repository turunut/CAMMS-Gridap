module modSubroutines
  export get_tags, get_E_CellField, ctˣ, step_field, EdgeTriangulation

  using Gridap.Geometry
  using Gridap.Arrays
  using Gridap.Fields
  using Gridap.CellData
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
    D = num_cell_dims(model)
    labeling = get_face_labeling(model)
    face_to_mask = get_face_mask(labeling,["tag_17"],D-2)
    face_to_bgface = findall(face_to_mask)
    return EdgeTriangulation(model,face_to_bgface)
  end

  function EdgeTriangulation(
    model::DiscreteModel,
    face_to_bgface::AbstractVector{<:Integer})

    D = num_cell_dims(model)
    topo = get_grid_topology(model)
    bgface_grid = Grid(ReferenceFE{D-2},model)
    bgface_to_lcell = Fill(1,num_faces(model,D-2))

    face_grid = view(bgface_grid,face_to_bgface)
    cell_grid = get_grid(model)
    glue = Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
    trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)

    BoundaryTriangulation(trian,glue)
  end

end