module GridapGiD
  export GiDDiscreteModel

  using Gridap.TensorValues
  using Gridap.ReferenceFEs

  function GiDDiscreteModel(fileName::String)

    Nodes    = readFileNodes(fileName)
    Elements, ElementsMap, ElementsType, ElementsTypeMap = readFileModel(fileName)
    d_to_dface_to_entity, tag_to_entities, tag_to_name = readSets(fileName)

    return
  end

  function readFileNodes(problemName::String)
    Nodes = Vector{VectorValue{}}()
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
    ElementsType = Vector{Polytope}()
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

            println( splittedLine[2:(1+numNodes)] )

            for (inum, num) in enumerate(splittedLine[2:(1+numNodes)])
              #push!(nodeList, parse(Int8, num))
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

  function readSets(problemName::String)
    d_to_dface_to_entity = [[1,3,3,3,3,3,3,3,2],[3,3,3,3,3,3,3,3]] # nodes edges
    tag_to_entities = [[3],[1,2],[1],[2]]
    tag_to_name     = ["interior","boundary","bc1","bc2"]
  
    open(problemName) do f
      while ! eof(f)
        
        splittedLine = readAndSplit(f)

        if splittedLine[1] == "set_definition"
          splittedLine = readAndSplit(f)
          while splittedLine[1] != "set_end"
            ielem = parse(Int64, splittedLine[1])
            imat  = parse(Int64, splittedLine[2])
            splittedLine = readAndSplit(f)
          end
        end

        if splittedLine[1] == "body_definition"
          
          splittedLine = readAndSplit(f)
          while splittedLine[1] != "body_end"
            ielem = parse(Int64, splittedLine[1])
            imat  = parse(Int64, splittedLine[2])
            splittedLine = readAndSplit(f)
          end
        end

      end
    end

    return Entities
  end

  function getElemType(geometry::String)
    println(geometry)
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
