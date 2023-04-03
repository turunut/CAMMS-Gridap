module modModel
  export Model, PlaneStress, Timoshenko, Reissner
  export computeCT, ε, σ, γ
  export get_CT_CellField

  using modCT
  using Gridap.CellData
  using Gridap.Fields
  using Gridap.TensorValues
  using Gridap.Geometry
  using Gridap.Arrays

  abstract type Model end

  mutable struct PlaneStress <: Model end

  mutable struct Timoshenko <: Model
    b::Float64
    z2::Float64
    z1::Float64
  end

  mutable struct Reissner <: Model
    z2::Float64
    z1::Float64
  end

  function computeCT(mod::PlaneStress, ct::CT)
    CT_2D = modCT.compute2D(ct)

    ct = SymFourthOrderTensorValue( (CT_2D[1,1],        0.0, CT_2D[1,2],
                                            0.0, CT_2D[3,3],        0.0, 
                                     CT_2D[2,1],        0.0, CT_2D[2,2],) )

    return ct
  end

  function computeCT(mod::Timoshenko, ct::CT)
    CT_1D = modCT.compute1D(ct)

    A = SymFourthOrderTensorValue( ( CT_1D[1,1]*(1/1)*mod.b*(mod.z2^1-mod.z1^1) ) )
    B = SymFourthOrderTensorValue( ( CT_1D[1,1]*(1/2)*mod.b*(mod.z2^2-mod.z1^2) ) )
    D = SymFourthOrderTensorValue( ( CT_1D[1,1]*(1/3)*mod.b*(mod.z2^3-mod.z1^3) ) )
    S = SymFourthOrderTensorValue( ( CT_1D[2,2]*(1/1)*mod.b*(mod.z2^1-mod.z1^1)*(5/6) ) )

    return [A, B, D, S]
  end

  function computeCT(mod::Reissner, ct::CT)
    #CT_3D = modCT.compute(ct)

    IParray = zeros((3,3))
    IParray[1:3,1:3] = CT_3D[1:3,1:3]

    OParray = zeros((2,2))
    OParray[1:2,1:2] = CT_3D[4:5,4:5]

    A = SymTensorValue( ( IParray[1,1]*(1/1)*b*(z2^1+z1^1) ) )
    B = SymTensorValue( ( IParray[1,1]*(1/2)*b*(z2^2+z1^2) ) )
    D = SymTensorValue( ( IParray[1,1]*(1/3)*b*(z2^3+z1^3) ) )
    S = SymTensorValue( ( OParray[1,1]*(1/1)*b*(z2^1+z1^1) ) )

    return [A, B, D, S]
  end

  function σ(ct::CellField, ϵ::CellField)
    return ct ⊙ ϵ
  end

  function ε(u::CellField)
    return ε(u)
  end

  tosym_op(θ) = SymTensorValue(θ[1])
  tosym(θ::CellField) = Operation(tosym_op)(θ)
  function γₘ(mod::Timoshenko, ω::CellField, θ::CellField)
    return ε(ω) + tosym(θ)
  end

  function γ(mod::Reissner,   ω::CellField, θ::CellField)
    return ε(ω)⋅VectorValue(1.0) + θ
  end

  function get_CT_CellField(mod::Model, listCTs, tags, Ω::Triangulation)
    cts = CompressedArray(listCTs[1,:], tags)
    return CellField(cts,Ω)
  end

  function get_CT_CellField(mod::Timoshenko, listCTs, tags, Ω::Triangulation)
    cts1 = CompressedArray(listCTs[1,:], tags)
    cts2 = CompressedArray(listCTs[2,:], tags)
    cts3 = CompressedArray(listCTs[3,:], tags)
    cts4 = CompressedArray(listCTs[4,:], tags)
    return [CellField(cts1,Ω), CellField(cts2,Ω), CellField(cts3,Ω), CellField(cts4,Ω)]
  end

end