#using Gridap
#using Gridap.FESpaces
#using Gridap.CellData
#
#using FillArrays


order = 2
interfaces = [1,1,1,1]
ext_discrt = [3,3,3,3]
intrf_map = []
for (iintrf,intrf) in enumerate(interfaces)
  if order == 1
    mapping = 1:ext_discrt[iintrf] + 1
  elseif order == 2
    mapping1 = 1:ext_discrt[iintrf] + 1
    mapping2 = (length(mapping1)+1):(length(mapping1)*2-1)
    mapping = zeros(Int, length(mapping1) + length(mapping2))
    mapping[1] = mapping1[1]
    for i in 1:1:length(mapping2)
      mapping[ (i)*2 + 0 ] = mapping2[i]
      mapping[ (i)*2 + 1 ] = mapping1[i]
    end
  else
    println("NOT IMPLEMENTED")
  end
  push!(intrf_map, mapping)
end





















order = 1

# Create a Cartesian grid
domainF    = (0, 1)
partitionF = (5)
modelF = CartesianDiscreteModel(domainF, partitionF)
ΩF = Triangulation(modelF)

# Create a Cartesian grid
domainC    = (0, 1)
partitionC = (3)
modelC = CartesianDiscreteModel(domainC, partitionC)
ΩC = Triangulation(modelC)



reffe = ReferenceFE(lagrangian,Float64,order)

VF = FESpace(ΩF,reffe,conformity=:H1)
UF = TrialFESpace(VF)

#reffe = ReferenceFE(lagrangian,Float64,order)
VC = FESpace(ΩC,reffe,conformity=:H1)
UC = TrialFESpace(VC)



#xF = zero_free_values(UF); xF[3] = 1.0
#uF = FEFunction(UF,xF)

xC = zero_free_values(UC); xC[3] = 1.0
uC = FEFunction(UC,xC)

uCi = Interpolable(uC)

uF = interpolate(uCi,UF)


#_get_y(x) = VectorValue(x[1])
#function π_Λe_Γc(f::CellField, Γc::Triangulation)
#    _data = CellData.get_data(f)
#    _cellmap = Fill(Broadcasting(_get_y),length(_data))
#    data = lazy_map(∘,_data,_cellmap)
#    return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
#end



xF = zero_free_values(UF)
uF = π_Λe_Γc(fC,ΩF)

for i in 0:0.02:1
  println( i, "  ", uF( Point(i) ) )
end
