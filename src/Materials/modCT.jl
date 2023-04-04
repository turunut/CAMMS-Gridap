module modCT
  export CT, CT_Isotrop, CT_Orthotrop
  export compute1D, compute2D, compute3D,
         compute1DoutPlane,
         compute2DoutPlane
  
  using Gridap.TensorValues
  
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
    return SymFourthOrderTensorValue( ( ct.E ) )
  end

  function compute1DoutPlane(ct::CT_Isotrop)
    G = ct.E/(2*(1+ct.n))
    return TensorValue(G)
  end

  function compute2D(ct::CT_Isotrop)
    Sarray = zeros((3,3))
    G = ct.E/(2*(1+ct.n))
    Sarray[1,1] =     1/ct.E; Sarray[1,2] = -ct.n/ct.E
    Sarray[2,1] = -ct.n/ct.E; Sarray[2,2] =     1/ct.E
    Sarray[3,3] =     1/G
    
    Carray = inv(Sarray)

    ct = SymFourthOrderTensorValue( (Carray[1,1],         0.0, Carray[1,2],
                                             0.0, Carray[3,3],         0.0, 
                                     Carray[2,1],         0.0, Carray[2,2],) )

    return ct
  end

  function compute2DoutPlane(ct::CT_Isotrop)
    Sarray = zeros((2,2))
    G = ct.E/(2*(1+ct.n))
    Sarray[1,1] = 1/G
    Sarray[2,2] = 1/G
    
    Carray = inv(Sarray)

    #ct = SymFourthOrderTensorValue( (Carray[1,1], 0.0, Carray[2,2]) )
    
    # 11 = xz; 22 = yz
    # | Sxz 0.0 | = xz
    # | 0.0 Syz | = yz
    # SymFourthOrderTensorValue(1111,1121,1122, 2111,2121,2122, 2211,2221,2222)
    # SymFourthOrderTensorValue(1111, 0.0,1122,  0.0, 0.0, 0.0, 2211, 0.0,2222)
    ct = TensorValue(Carray)

    return ct
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
    
    Carray = inv(Sarray)

    ct = SymFourthOrderTensorValue( (Carray[1,1],         0.0, Carray[1,2],
                                             0.0, Carray[3,3],         0.0, 
                                     Carray[2,1],         0.0, Carray[2,2],) )

    return ct
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