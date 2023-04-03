module modCT
  export CT, CT_Isotrop, CT_Orthotrop
  export compute1D, compute2D, compute3D
  
  abstract type CT end

  struct CT_Isotrop  <: CT
    E::Float64
    n::Float64
  end

  struct CT_Orthotrop <: CT
    E11::Float64
    E22::Float64
    E33::Float64
    G12::Float64
    G13::Float64
    G23::Float64
    nu12::Float64
    nu13::Float64
    nu23::Float64
  end

  function compute1D(ct::CT_Isotrop)
    Sarray = zeros((2,2))
    G = ct.E/(2*(1+ct.n))
    Sarray[1,1] = 1/ct.E; Sarray[2,2] = 1/G
    return inv(Sarray)
  end

  function compute2D(ct::CT_Isotrop)
    Sarray = zeros((3,3))
    G = ct.E/(2*(1+ct.n))
    Sarray[1,1] =     1/ct.E; Sarray[1,2] = -ct.n/ct.E
    Sarray[2,1] = -ct.n/ct.E; Sarray[2,2] =     1/ct.E
    Sarray[3,3] =     1/G
    return inv(Sarray)
  end

  function compute3D(ct::CT_Isotrop)
    Sarray = zeros((6,6))
    G = ct.E/(2*(1+ct.n))
    Sarray[1,1] =     1/ct.E; Sarray[1,2] = -ct.n/ct.E; Sarray[1,3] = -ct.n/ct.E
    Sarray[2,1] = -ct.n/ct.E; Sarray[2,2] =     1/ct.E; Sarray[2,3] = -ct.n/ct.E
    Sarray[3,1] = -ct.n/ct.E; Sarray[3,2] = -ct.n/ct.E; Sarray[3,3] =     1/ct.E
    Sarray[4,4] =     1/G
    Sarray[5,5] =     1/G
    Sarray[6,6] =     1/G
    return inv(Sarray)
  end

  function compute3D(ct::CT_Orthotrop)
    Sarray = zeros((6,6))
    n21 = 1.0
    n31 = 1.0
    n32 = 1.0
    Sarray[1,1] =       1/ct.E11; Sarray[1,2] = -ct.n21/ct.E22; Sarray[1,3] = -ct.n31/ct.E33
    Sarray[2,1] = -ct.n12/ct.E11; Sarray[2,2] =       1/ct.E22; Sarray[2,3] = -ct.n32/ct.E33
    Sarray[3,1] = -ct.n13/ct.E11; Sarray[3,2] = -ct.n23/ct.E22; Sarray[3,3] =       1/ct.E33
    Sarray[4,4] =       1/ct.G12
    Sarray[5,5] =       1/ct.G13
    Sarray[6,6] =       1/ct.G23
    return inv(Sarray)
  end

end