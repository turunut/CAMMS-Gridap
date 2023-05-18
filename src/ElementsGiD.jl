module ElementsGiD
  export factoryElemGiD, getGridapType, sortConn, getElementsTypeMap,
         ElementGiD, linear, quadrilateral, hexahedra,
         linear2N,
         quadrilateral4N,
         hexahedra8N

  using Gridap.Geometry
  using Gridap.ReferenceFEs
  #using DataStructures

  abstract type ElementGiD end

  abstract type surface <: ElementGiD end
  abstract type volume  <: ElementGiD end
  
  struct point                <: ElementGiD end
  abstract type linear        <: ElementGiD end
  abstract type quadrilateral <: surface    end
  abstract type hexahedra     <: volume    end

  struct linear2N        <: linear        end
  struct quadrilateral4N <: quadrilateral end
  struct hexahedra8N     <: hexahedra end

  function get_set_nodes(elem::ElementGiD); return Set(elem.nodes); end

  function get_boun_type(elem::point);           return nothing;    end
  function get_boun_type(elem::linear2N);        return point();    end
  function get_boun_type(elem::quadrilateral4N); return linear2N(); end
  function get_boun_type(elem::hexahedra8N);     return quadrilateral4N(); end

  function get_type_index(elem::Any);     return 0; end
  function get_type_index(elem::point);   return 1; end
  function get_type_index(elem::linear);  return 2; end
  function get_type_index(elem::surface); return 3; end
  function get_type_index(elem::volume);  return 4; end

  function get_bouns(elem::point, nodes::Vector{Integer})
    return Vector{Vector{Integer}}()
  end
  function get_bouns(elem::linear2N, nodes::Vector{Integer})
    bounds = Vector{Vector{Integer}}([[nodes[1]],
                                      [nodes[2]]])
    return bounds
  end
  function get_bouns(elem::quadrilateral4N, nodes::Vector{Integer})
    bounds = Vector{Vector{Integer}}([[nodes[1], nodes[2]],
                                      [nodes[3], nodes[4]],
                                      [nodes[1], nodes[3]],
                                      [nodes[2], nodes[4]]])
    return bounds
  end
  function get_bouns(elem::hexahedra8N, nodes::Vector{Integer})
    bounds = Vector{Vector{Integer}}([[nodes[1], nodes[2], nodes[3], nodes[4]],
                                      [nodes[5], nodes[6], nodes[7], nodes[8]],
                                      [nodes[1], nodes[2], nodes[5], nodes[6]],
                                      [nodes[3], nodes[4], nodes[7], nodes[8]],
                                      [nodes[1], nodes[3], nodes[5], nodes[7]],
                                      [nodes[2], nodes[4], nodes[6], nodes[8]]])
    return bounds
  end

  function compute_bouns(elemType::ElementGiD, listsOfbounds::Vector{Vector{Set{Integer}}}, nodes::Vector{Integer})
    boutyp = get_boun_type(elemType)
    index = get_type_index(boutyp)
    for boun in get_bouns(elemType, nodes)
      if Set(boun) ∉ listsOfbounds[index]
        push!(listsOfbounds[index], Set(boun) )
      end
      compute_bouns(boutyp, listsOfbounds, boun)
    end
    return listsOfbounds
  end

  function factoryElemGiD(elemType::String, numNodes::Integer)
    if     elemType == "linear"
      if numNodes == 2
        return linear2N()
      end
    elseif elemType == "quadrilateral"
      if numNodes == 4
        return quadrilateral4N()
      end
    elseif elemType == "hexahedra"
      if numNodes == 8
        return hexahedra8N()
      end
    end
  end

  function getGridapType(type::linear);        return SEGMENT; end
  function getGridapType(type::quadrilateral); return QUAD;    end
  function getGridapType(type::hexahedra);     return HEX;     end

  function sortConn(elemgid::linear2N,        con::Vector{Integer})
    return con
  end
  function sortConn(elemgid::quadrilateral4N, con::Vector{Integer})
    con[1], con[2], con[4], con[3] = con[1], con[2], con[3], con[4]
    return con
  end
  function sortConn(elemgid::hexahedra8N, con::Vector{Integer})
    con[1], con[2], con[4], con[3], con[5], con[6], con[8], con[7] = con[1], con[2], con[3], con[4], con[5], con[6], con[7], con[8]
    return con
  end

  function getElementsTypeMap(Elements::Vector{Integer},
                              ElementsMap::Vector{Integer},
                              ElementsTypeGiD::Vector{ElementGiD},
                              ElementsTypeMap::Vector{Integer})

    entities_lists = Vector{Vector{Set{Integer}}}([[],[],[],[]])

    for (ielem, elem) in enumerate(ElementsTypeMap)
      elemgid = ElementsTypeGiD[elem]
      nodes   = Elements[ElementsMap[ielem]:ElementsMap[ielem+1]-1]
      index = get_type_index(elemgid)

      entities_lists = compute_bouns(elemgid, entities_lists, nodes)

      push!(entities_lists[index], Set(nodes))
    end

    # Els nodes han destar ordenats...
    for i in 1:length(entities_lists[1])
      entities_lists[1][i] = Set(i)
    end

    for i in reverse(1:4)
      if length(entities_lists[i]) == 0
        pop!(entities_lists)
      end
    end

    return entities_lists
  end
    
end