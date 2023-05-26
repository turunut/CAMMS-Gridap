module modSubroutines
  export get_tags, get_E_CellField, ctˣ, step_field

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

end