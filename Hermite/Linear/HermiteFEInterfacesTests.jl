module HermiteFEInterfacesTests

using Test
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs

#struct PointValue{P} <: Dof
#  point::P
#end

D = 2
T = Float64
order = 1
prebasis = MonomialBasis{D}(T,order)

polytope = SEGMENT
x = get_vertex_coordinates(polytope)

dofsLag = LagrangianDofBasis(T,x)

##################################################################

struct HermiteDofBasis{P,V} <: AbstractVector{ReferenceFEs.PointValue{P}}
  nodes::Vector{P}
  dof_to_node::Vector{Int}
  dof_to_comp::Vector{Int}
  node_and_comp_to_dof::Vector{V}
end

function _generate_dof_layout_node_major(::Type{<:Real},nnodes::Integer)
  ndofs = nnodes
  dof_to_comp = [1, 1, 1, 1]
  dof_to_node = [1, 1, 2, 2]
  node_and_comp_to_dof = [[1,3],[2,4]]
  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
end

function HermiteDofBasis(::Type{T},nodes::Vector{<:Point}) where T
  r = _generate_dof_layout_node_major(T,length(nodes))
  HermiteDofBasis(nodes,r...)
end

function ReferenceFEs.evaluate!(cache,b::HermiteDofBasis,field)
  c, cf = cache
  vals = evaluate!(cf,field,b.nodes)
  ndofs = length(b.dof_to_node)
  T = eltype(vals)
  ncomps = num_components(T)
  #@check ncomps == num_components(eltype(b.node_and_comp_to_dof)) """\n
  #Unable to evaluate LagrangianDofBasis. The number of components of the
  #given Field does not match with the LagrangianDofBasis.
  #
  #If you are trying to interpolate a function on a FESpace make sure that
  #both objects have the same value type.
  #
  #For instance, trying to interpolate a vector-valued funciton on a scalar-valued FE space
  #would raise this error.
  #"""
  _evaluate_herm_dof!(c,vals,b.node_and_comp_to_dof,ndofs,ncomps)
end

function _evaluate_herm_dof!(c::AbstractVector,node_comp_to_val,node_and_comp_to_dof,ndofs,ncomps)
  setsize!(c,(ndofs,))
  r = c.array
  for node in LinearIndices(node_and_comp_to_dof)
    comp_to_dof = node_and_comp_to_dof[node]
    comp_to_val = node_comp_to_val[node]
    for comp in 1:ncomps
      dof = comp_to_dof[comp]
      val = comp_to_val[comp]
      r[dof] = val
    end
  end
  r
end



##################################################################

dofs = HermiteDofBasis(T,x)


v = 3.0
f = GenericField(x->v*x[1])
fx = evaluate(f,x)
test_dof_array(dofs,f,fx)

ndof = 4
b = fill(f,ndof)
bx = evaluate(b,x)
test_dof_array(dofs,b,bx)

ndofs = 2

#face_own_dofs = [[1],[2],[3],[4],Int[],Int[],Int[],Int[],Int[]]
#a = [1]
#b = Int[]
#face_own_dofs_permutations = [[a],[a],[a],[a],[b,b],[b,b],[b,b],[b,b],fill(b,8)]
face_dofs = [[1],[2],[1,2]]
metadata = nothing

struct MockConformity <: Conformity end

conf = MockConformity()
reffe = GenericRefFE{typeof(conf)}(
  ndofs, polytope, prebasis, dofs, conf, metadata ,face_dofs)

function ReferenceFEs.get_face_own_dofs(reffe::GenericRefFE{MockConformity},::MockConformity)
  face_own_dofs
end

function ReferenceFEs.get_face_own_dofs_permutations(reffe::GenericRefFE{MockConformity},::MockConformity)
  face_own_dofs_permutations
end

test_reference_fe(reffe)

@test get_face_own_dofs(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3, 4]]
@test get_face_own_dofs_permutations(reffe,L2Conformity()) == [[[]], [[]], [[]], [[]], [[]], [[]], [[]], [[]], [[1, 2, 3, 4]]]

shapefuns = get_shapefuns(reffe)

@test evaluate(shapefuns,x) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

@test get_own_dofs_permutations(reffe) == face_own_dofs_permutations[end]

end # module
