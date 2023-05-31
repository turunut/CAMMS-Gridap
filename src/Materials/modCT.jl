module modCT
  export CT, CT_Isotrop, CT_Orthotrop
  export compute1D, compute2D, compute3D,
         compute1DoutPlane,
         compute2DoutPlane,
         getE
  
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
    n12::Float64
    n13::Float64
    n23::Float64
  end

  function compute1D(ct::CT_Isotrop)
    return TensorValue(ct.E)
  end

  function compute1DoutPlane(ct::CT_Isotrop)
    G = ct.E/(2*(1+ct.n))
    return TensorValue(G)
  end

  function compute2D(ct::CT_Isotrop, array_form::Bool)
    Sarray = zeros((3,3))
    G = ct.E/(2*(1+ct.n))
    Sarray[1,1] =     1/ct.E; Sarray[1,2] = -ct.n/ct.E
    Sarray[2,1] = -ct.n/ct.E; Sarray[2,2] =     1/ct.E
    Sarray[3,3] =     1/G
    
    CTarray = inv(Sarray)
    
    if array_form ; CTtensor = CT_array2array_2D(CTarray)
    else;           CTtensor = CT_array2tensor_2D(CTarray)
    end

    return CTtensor
  end

  function compute2D(ct::CT_Orthotrop)
    Sarray = zeros((3,3))
    n21 = ct.E22*(ct.n12/ct.E11)
    Sarray[1,1] =       1/ct.E11; Sarray[1,2] =    -n21/ct.E22
    Sarray[2,1] = -ct.n12/ct.E11; Sarray[2,2] =       1/ct.E22
    Sarray[3,3] =     1/G12
    
    CTarray = inv(Sarray)

    CTtensor = CT_array2tensor_2D(CTarray)

    return CTtensor
  end

  function compute2DoutPlane(ct::CT_Isotrop)
    Sarray = zeros((2,2))
    G = ct.E/(2*(1+ct.n))
    Sarray[1,1] = 1/G
    Sarray[2,2] = 1/G
    
    Carray = inv(Sarray)

    # 11 = xz; 22 = yz
    # | Sxz 0.0 | = xz
    # | 0.0 Syz | = yz
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
    
    CTarray = inv(Sarray)

    CTtensor = CT_array2tensor_3D(CTarray)

    return CTtensor
  end

  function compute3D(ct::CT_Orthotrop)
    Sarray = zeros((6,6))
    n21 = ct.E22*(ct.n12/ct.E11)
    n31 = ct.E33*(ct.n13/ct.E11)
    n32 = ct.E33*(ct.n23/ct.E22)
    Sarray[1,1] =       1/ct.E11; Sarray[1,2] =    -n21/ct.E22; Sarray[1,3] =    -n31/ct.E33
    Sarray[2,1] = -ct.n12/ct.E11; Sarray[2,2] =       1/ct.E22; Sarray[2,3] =    -n32/ct.E33
    Sarray[3,1] = -ct.n13/ct.E11; Sarray[3,2] = -ct.n23/ct.E22; Sarray[3,3] =       1/ct.E33
    Sarray[4,4] =       1/ct.G12
    Sarray[5,5] =       1/ct.G13
    Sarray[6,6] =       1/ct.G23
    return inv(Sarray)
  end

  function getE(ct::CT_Isotrop, flag::Vector{Float64})
    return ct.E
  end

  function CT_array2tensor_2D(array::Array{Float64})
    # 11 12 22
    #                             1      3      2
    # SymFourthOrderTensorValue(11_11, 11_12, 11_22,   1   
    #                           12_11, 12_12, 12_22,   3
    #                           22_11, 22_12, 22_22,   2
    #
    tensor = SymFourthOrderTensorValue( (array[1,1], array[1,3], array[1,2],
                                         array[3,1], array[3,3], array[3,2], 
                                         array[2,1], array[2,3], array[2,2],) )
    return tensor
  end

  function CT_array2array_2D(array::Array{Float64})
    tensor = TensorValue{3,3}( array[1,1], array[2,1], array[3,1],
                               array[1,2], array[2,2], array[3,2], 
                               array[1,3], array[2,3], array[3,3] )
    return tensor
  end

  function CT_array2tensor_3D(a::Array{Float64})
    # 11 12 13 22 23 33
    #                             1      4      5      2      6      3
    # SymFourthOrderTensorValue(11_11, 11_12, 11_13, 11_22, 11_23, 11_33,   1   
    #                           12_11, 12_12, 12_13, 12_22, 12_23, 12_33,   4
    #                           13_11, 13_12, 13_13, 13_22, 13_23, 13_33,   5
    #                           22_11, 22_12, 22_13, 22_22, 22_23, 22_33,   2
    #                           23_11, 23_12, 23_13, 23_22, 23_23, 23_33,   6
    #                           33_11, 33_12, 33_13, 33_22, 33_23, 33_33)   3
    #
    t = SymFourthOrderTensorValue( (a[1,1], a[1,4], a[1,5], a[1,2], a[1,6], a[1,3],
                                    a[4,1], a[4,4], a[4,5], a[4,2], a[4,6], a[4,3],
                                    a[5,1], a[5,4], a[5,5], a[5,2], a[5,6], a[5,3],
                                    a[2,1], a[2,4], a[2,5], a[2,2], a[2,6], a[2,3],
                                    a[6,1], a[6,4], a[6,5], a[6,2], a[6,6], a[6,3],
                                    a[3,1], a[3,4], a[3,5], a[3,2], a[3,6], a[3,3]) )
    return t
  end

end