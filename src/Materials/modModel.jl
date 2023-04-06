module modModel
  export Model, Beam, Shell, PlaneStress, Solid, Timoshenko, TimoshenkoLayout, Reissner
  export computeCT, ∂, σ, σₑ, γ, γγ, toten, tovec
  export get_CT_CellField

  using modCT
  using Gridap.CellData
  using Gridap.Fields
  using Gridap.TensorValues
  using Gridap.Geometry
  using Gridap.Arrays

  abstract type Model end

  abstract type Beam  <: Model end
  abstract type Shell <: Model end

  mutable struct PlaneStress <: Model end

  mutable struct Solid <: Model end

  mutable struct Timoshenko <: Beam
    b::Float64
    z2::Float64
    z1::Float64
  end

  mutable struct TimoshenkoLayout <: Beam
    layers::Vector{Timoshenko}
    cts::Vector{CT}
  end

  mutable struct Reissner <: Shell
    z2::Float64
    z1::Float64
  end

  function computeCT(mod::PlaneStress, ct::CT)
    CT_2D = modCT.compute2D(ct)

    return [CT_2D]
  end

  function computeCT(mod::Solid, ct::CT)
    CT_3D = modCT.compute3D(ct)

    return [CT_3D]
  end

  function computeCT(mod::Timoshenko, ct::CT)
    CT_1Dip = modCT.compute1D(ct)
    CT_1Dop = modCT.compute1DoutPlane(ct)

    A = CT_1Dip*(1/1)*mod.b*(mod.z2^1-mod.z1^1)
    B = CT_1Dip*(1/2)*mod.b*(mod.z2^2-mod.z1^2)
    D = CT_1Dip*(1/3)*mod.b*(mod.z2^3-mod.z1^3)
    S = CT_1Dop*(1/1)*mod.b*(mod.z2^1-mod.z1^1)*(5/6)

    return [A, B, D, S]
  end

  function computeCT(mod::TimoshenkoLayout, ct::CT)
    for (layer, ct) in zip(mod.layers)
      CT_1Dip = layer.compute1D(ct)
      CT_1Dop = layer.compute1DoutPlane(ct)
  
      A += CT_1Dip*(1/1)*mod.b*(mod.z2^1-mod.z1^1)
      B += CT_1Dip*(1/2)*mod.b*(mod.z2^2-mod.z1^2)
      D += CT_1Dip*(1/3)*mod.b*(mod.z2^3-mod.z1^3)
      S += CT_1Dop*(1/1)*mod.b*(mod.z2^1-mod.z1^1)*(5/6)
    end

    return [A, B, D, S]
  end

  function computeCT(mod::Reissner, ct::CT)
    CT_2Dip = modCT.compute2D(ct)
    CT_2Dop = modCT.compute2DoutPlane(ct)

    A = CT_2Dip*(1/1)*(mod.z2^1-mod.z1^1)
    B = CT_2Dip*(1/2)*(mod.z2^2-mod.z1^2)
    D = CT_2Dip*(1/3)*(mod.z2^3-mod.z1^3)
    S = CT_2Dop*(1/1)*(mod.z2^1-mod.z1^1)#*(5/6)

    return [A, B, D, S]
  end

  function σ(ct::CellField, ϵ::CellField)
    return ct ⊙ ϵ
  end
  
  function σₑ(ct::CellField, ϵ::CellField)
    return ∘(ct ⋅ ϵ)
  end

  function ∂(u::CellField)
    return ε(u)
  end

  toten_op(a) = TensorValue(a[1])
  toten(a::CellField) = Operation(toten_op)(a)
  tovec_op(a) = VectorValue(a[1])
  tovec(a::CellField) = Operation(tovec_op)(a)


  γ_op_TB(∂ω, θ) = VectorValue(∂ω[1] - θ[1])
  function γ(mod::Timoshenko, ω::CellField, θ::CellField)
    return Operation(γ_op_TB)(gradient(ω), θ)
  end
  γ_op_TB(∂ω, θ) = VectorValue(∂ω[1] - θ[1])
  function γγ(mod::Timoshenko, ω::CellField, θ::CellField)
    return Operation(γ_op_TB)(ω, θ)
  end

  γ_op_RS(∂ω, θ) = VectorValue(∂ω[1] - θ[1], ∂ω[2] - θ[2])
  function γ(mod::Reissner,   ω::CellField, θ::CellField)
    return Operation(γ_op_RS)(gradient(ω), θ)             # ∂(ω) + tosym(θ)
  end

  function get_CT_CellField(mod::Model, listCTs, tags, Ω::Triangulation)
    cts = CompressedArray(listCTs[1,:], tags)
    return [CellField(cts,Ω)]
  end

  function get_CT_CellField(mod::Timoshenko, listCTs, tags, Ω::Triangulation)
    cts1 = CompressedArray(listCTs[1,:], tags)
    cts2 = CompressedArray(listCTs[2,:], tags)
    cts3 = CompressedArray(listCTs[3,:], tags)
    cts4 = CompressedArray(listCTs[4,:], tags)
    return [CellField(cts1,Ω), CellField(cts2,Ω), CellField(cts3,Ω), CellField(cts4,Ω)]
  end

  function get_CT_CellField(mod::Reissner, listCTs, tags, Ω::Triangulation)
    cts1 = CompressedArray(listCTs[1,:], tags)
    cts2 = CompressedArray(listCTs[2,:], tags)
    cts3 = CompressedArray(listCTs[3,:], tags)
    cts4 = CompressedArray(listCTs[4,:], tags)
    return [CellField(cts1,Ω), CellField(cts2,Ω), CellField(cts3,Ω), CellField(cts4,Ω)]
  end

end