using Gridap
using Gridap.FESpaces

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


fC = interpolate_everywhere(uC, FESpace(ΩC,reffe))
fF = interpolate_everywhere(uC, FESpace(ΩF,reffe))



for i in 0:0.02:1
  println( i, "  ", uC( Point(i) ) )
end

uF = interpolate(uC, VF)

for i in 0:0.02:1
  println( i, "  ",uF( Point(i) ) )
end
