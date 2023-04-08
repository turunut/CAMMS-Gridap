module GridapGiD
  export GiDDiscreteModel

  using Gridap.TensorValues
  using Gridap.ReferenceFEs
  using Gridap.Arrays
  using Gridap.Geometry

  function GiDDiscreteModel(fileName::String)

    nodes    = readFileNodes(fileName*".msh")
    Elements, ElementsMap, ElementsType, ElementsTypeMap = readFileModel(fileName*".msh")
    d_to_dface_to_entity = Vector{Vector{Int8}}()
    push!(d_to_dface_to_entity, ones(Int8, length(nodes)))
    push!(d_to_dface_to_entity, ones(Int8, length(ElementsMap) - 1))
    d_to_dface_to_entity, tag_to_entities, tag_to_name = readSets(fileName*".set", d_to_dface_to_entity)

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
    Nodes = Vector{VectorValue{2, Float64}}()
    ndime = 0

    open(problemName) do f
      while ! eof(f) 
        
        splittedLine = readAndSplit(f)

        if splittedLine[1] == "mesh"
          ndime = parse(Int64, splittedLine[3])
        end

        if splittedLine[1] == "coordinates"
          println("Reading coordinates.")

          splittedLine = readAndSplit(f)

          while splittedLine[1] != "end"

            idNode = parse(Int64, splittedLine[1])

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
    
    push!(ElementsMap, 1)

    open(problemName) do f
      numNodes = 0
      elemType = ""
      while ! eof(f)
        splittedLine = readAndSplit(f)

        if splittedLine[1] == "mesh"
          elemType = getElemType(string(splittedLine[5]))
          push!(ElementsType, elemType)
          numNodes = parse(Int8, splittedLine[7])
        end

        if splittedLine[1] == "elements"
          println("  Reading elements.")

          splittedLine = readAndSplit(f)

          nodeList = Vector{Int8}(undef, numNodes)
          while splittedLine[1] != "end"

            for (inum, num) in enumerate(splittedLine[2:(1+numNodes)])
              nodeList[inum] = parse(Int8, num)
            end

            append!(Elements,        nodeList)
            append!(ElementsMap,     ElementsMap[end] + numNodes)
            append!(ElementsTypeMap, length(ElementsType))

            splittedLine = readAndSplit(f)

          end

        end

      end
    end

    return Elements, ElementsMap, ElementsType, ElementsTypeMap
  end

  function readSets(problemName::String, d_to_dface_to_entity::Vector{Vector{Int8}})
    tag_to_name     = ["interior"]
    tag_to_entities = [[1]       ]
    ibody = 2
  
    open(problemName) do f
      while ! eof(f)
        
        splittedLine = readAndSplit(f)

        if splittedLine[1] == "set_definition"
          println("  Reading set.")
          splittedLine = readAndSplit(f)
          while splittedLine[1] != "set_end"
            ielem = parse(Int8, splittedLine[1])
            imat  = parse(Int8, splittedLine[2])
            d_to_dface_to_entity[2][ielem] = imat
            splittedLine = readAndSplit(f)
          end
          println("  Finish reading set.")
        end

        if splittedLine[1] == "tag_definition"
          println("  Reading tag.")
          push!(tag_to_name, string(splittedLine[2]))
          push!(tag_to_entities, [ibody])
          splittedLine = readAndSplit(f)
          while splittedLine[1] != "tag_end"
            inode = parse(Int64, splittedLine[1])
            d_to_dface_to_entity[1][inode] = ibody
            splittedLine = readAndSplit(f)
          end
          ibody += 1
          println("  Finish reading tag.")
        end

      end
    end

    return d_to_dface_to_entity, tag_to_entities, tag_to_name
  end

  function getElemType(geometry::String)
    if geometry == "linear"
      return SEGMENT        
    end
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

end
