module GridapGiD
  export GiDDiscreteModel

  using Gridap.TensorValues
  using Gridap.ReferenceFEs
  using Gridap.Arrays
  using Gridap.Geometry
  using ElementsGiD

  function GiDDiscreteModel(fileName::String)

    nodes    = readFileNodes(fileName*".msh")
    Elements, ElementsMap, ElementsType, ElementsTypeMap, entities_lists = readFileModel(fileName*".msh")
    d_to_dface_to_entity = create_face_to_entity(entities_lists)
    d_to_dface_to_entity, tag_to_entities, tag_to_name = readSets(fileName*".set", entities_lists, d_to_dface_to_entity)

    c2n_map = Table(Elements, ElementsMap)

    cell_type = ElementsTypeMap
    polys = ElementsType
    reffes = map(p->LagrangianRefFE(Float64,p,1),polys)
    orientation = NonOriented()

    topo   = UnstructuredGridTopology(nodes,c2n_map,cell_type,polys,orientation)
    grid   = UnstructuredGrid(nodes,c2n_map,reffes,cell_type,orientation)
    
    labels = FaceLabeling(d_to_dface_to_entity, tag_to_entities, tag_to_name)

    model = UnstructuredDiscreteModel(grid,topo,labels)

    return model
  end

  function readFileNodes(problemName::String)
    ndime = 0

    open(problemName) do f
      while ! eof(f)
        splittedLine = readAndSplit(f)
        if splittedLine[1] == "mesh"
          ndime = parse(Int, splittedLine[3])
        end
      end
    end

    Nodes = Vector{VectorValue{ndime, Float64}}()

    open(problemName) do f
      while ! eof(f) 
        
        splittedLine = readAndSplit(f)

        if splittedLine[1] == "coordinates"
          println("Reading coordinates.")

          splittedLine = readAndSplit(f)

          while splittedLine[1] != "end"

            idNode = parse(Int, splittedLine[1])

            coords  = parse.(Float64, splittedLine[2:1+ndime])

            push!(Nodes, VectorValue(coords))

            splittedLine = readAndSplit(f)
          end
        end

      end
    end

    return Nodes
  end

  function readFileModel(problemName::String)
    Elements = Vector{Integer}()
    ElementsMap = Vector{Integer}()
    ElementsType = Vector{Polytope}() # Vector{ExtrusionPolytope{1}}() # 
    ElementsTypeMap = Vector{Integer}()
    
    ElementsTypeGiD = Vector{ElementGiD}()
    
    push!(ElementsMap, 1)

    open(problemName) do f
      numNodes = 0
      elemType = ""
      elemTypeString = ""
      elemGiD = ""
      while ! eof(f)
        splittedLine = readAndSplit(f)

        if splittedLine[1] == "mesh"
          elemTypeString = string(splittedLine[5])
          numNodes = parse(Int, splittedLine[7])

          elemGiD = factoryElemGiD(elemTypeString, numNodes)

          elemType = getGridapType(elemGiD)
          push!(ElementsTypeGiD, elemGiD)
          push!(ElementsType, elemType)
        end

        if splittedLine[1] == "elements"
          println("  Reading elements.")

          splittedLine = readAndSplit(f)

          nodeList = Vector{Integer}(undef, numNodes)
          while splittedLine[1] != "end"

            for (inum, num) in enumerate(splittedLine[2:(1+numNodes)])
              nodeList[inum] = parse(Int, num)
            end

            nodeList = sortConn(elemGiD, nodeList)

            append!(Elements,        nodeList)
            append!(ElementsMap,     ElementsMap[end] + numNodes)
            append!(ElementsTypeMap, length(ElementsType))

            splittedLine = readAndSplit(f)

          end

        end

      end
    end

    entities_lists = getElementsTypeMap(Elements, ElementsMap, ElementsTypeGiD, ElementsTypeMap)

    return Elements, ElementsMap, ElementsType, ElementsTypeMap, entities_lists
  end

  function readSets(problemName::String, entities_lists::Vector{Vector{Set{Integer}}}, d_to_dface_to_entity::Vector{Vector{Integer}})
    tag_to_name     = [] # = ["interior"]
    tag_to_entities = [] # = [[1]       ]
    ibody = 1
  
    open(problemName) do f
      while ! eof(f)
        
        splittedLine = readAndSplit(f)

        max_imat = 0

        if splittedLine[1] == "set_definition"
          println("  Reading set.")
          splittedLine = readAndSplit(f)
          while splittedLine[1] != "set_end"
            ielem = parse(Int, splittedLine[1])
            imat  = parse(Int, splittedLine[2])
            max_imat = max(imat, max_imat)
            d_to_dface_to_entity[end][ielem] = imat
            splittedLine = readAndSplit(f)
          end
          println("  Finish reading set.")

          for i=1:max_imat
            push!(tag_to_name, "mat_" * string(i))
            push!(tag_to_entities, [ibody])
            ibody += 1
          end
        end

        if splittedLine[1] == "tag_definition"
          println("  Reading tag.")
          index = parse(Int, splittedLine[2])

          if string(splittedLine[3]) in tag_to_name
            istri = findfirst(item -> item == string(splittedLine[3]), tag_to_name)
            itag = tag_to_entities[istri][:]
          else
            push!(tag_to_name, string(splittedLine[3]))
            push!(tag_to_entities, [ibody])
            itag = ibody
            ibody += 1
          end
          
          splittedLine = readAndSplit(f)
          while splittedLine[1] != "tag_end"
            nodes = Vector{Integer}()
            for num in splittedLine; push!(nodes, parse(Int, num)); end
            setNodes = Set(nodes)

            icon = findfirst(item -> item == setNodes, entities_lists[index])

            d_to_dface_to_entity[index][icon] = itag
            splittedLine = readAndSplit(f)
          end
          println("  Finish reading tag.")
        end

      end
    end

    return d_to_dface_to_entity, tag_to_entities, tag_to_name
  end

  function readAndSplit(file)
    line = readline(file)
    splitGivenLine(line, true)
  end

  function splitGivenLine(line::String, lowerFlag::Bool)
    if (lowerFlag); line = lowercase(line); end
    if (length(line) != 0)
      if line[1:1] == "#"
        return repeat(["emptyLine"], 5)
      else
        line = replace(line,":"=>" ")
        splittedLine = split(line)
        return splittedLine
      end
    else
      return repeat(["emptyLine"], 5)
    end
  end

  function create_face_to_entity(entities_lists::Vector{Vector{Set{Integer}}})
    d_to_dface_to_entity = Vector{Vector{Integer}}()

    for list in entities_lists
      push!(d_to_dface_to_entity, ones(Integer, length(list)))
    end

    return d_to_dface_to_entity
  end

end
