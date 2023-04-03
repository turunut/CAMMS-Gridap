module modBehaviour
  using modCT
  using modVoiOpe

  export Behaviour,
         Elastic,
         Laminate

  abstract type Behaviour end

  mutable struct Elastic  <: Behaviour
    ct::CT
    vo::VoiOpe
    Elastic() = new()
  end

  mutable struct Damage <: Behaviour
    ct::CT
    layers::Vector{Material}
    vo::VoiOpe
    Laminate() = new()
  end

  function computeStress(bhv::Elastic, mat, ε)
    ε.data

  end

end
