module modVoiOpe
  export VoiOpe, VoiOpePS, VoiOpeTS
  
  abstract type VoiOpe end

  struct VoiOpeTS <: VoiOpe
    h::Float64
    b::Float64
  end

  struct VoiOpePS  <: VoiOpe end

  function ten2vec(vo::VoiOpePS, u)
    gu = gradient(u)
    return VectorValue(1.0)
  end

end