module modMaterial_computeMatCT
  export computeMatCT

  using modCT
  using modMaterial

  function computeMatCT(mat::MaterialTSLaminate)
    arrayABD = zeros((3,3))
    for layer in mat.layers 
      computeMatCT(layer)
      arrayABD += layer.ct.CTarray
    end

    #Calculem matrix D respecte la linea neutre
    DdNL = 0.0
    for layer in mat.layers 
      zk   = layer.ct.z1
      zk1  = layer.ct.z2
      E    = layer.ct.E11  
      DdNL   += E*(zk1^3-zk^3)*(1/3)
    end

    factorK = 0.0
    SUM = 0.0
    for (indexR, LayerActual) in enumerate(mat.layers)
      R = 0.0
      for LayerAnterior in mat.layers[begin:indexR-1]
        z1 = LayerAnterior.ct.z1
        z2 = LayerAnterior.ct.z2
        E  = LayerAnterior.ct.E11
        G  = LayerAnterior.ct.G13
        R += (z2^2-z1^2)*E
      end
  
      z1 = LayerActual.ct.z1
      z2 = LayerActual.ct.z2
      E  = LayerActual.ct.E11
      G  = LayerActual.ct.G13
  
      R = - E*z1^2 + R
  
      Si = ((1/5)*E^2*z2^5+(2/3)*R*E*z2^3+R^2*z2) - ((1/5)*E^2*z1^5+(2/3)*R*E*z1^3+R^2*z1)
      Si = Si*1/(4*G)
      SUM += Si
    end
    factorK = ((DdNL^2)/(arrayABD[3,3]))*(1/SUM)

    arrayABD[3,3] = factorK*arrayABD[3,3]
    mat.ct = arrayABD

    printSparse(arrayABD, 3, 3)
    return nothing
  end

  function computeMatCT(mat::MaterialZZBeamLaminate)
    #computeBeta
    h = mat.layers[end].ct.z2 - mat.layers[1].ct.z1
    sum = 0.0
    for layer in mat.layers
      sum += (layer.ct.z2-layer.ct.z1)/layer.ct.G13
    end
    
    Ghat = h/sum
    
    for layer in mat.layers
      layer.ct.beta = (Ghat/layer.ct.G13) - 1
    end

    arrayABD = zeros((5,5))

    Da      = 0.0
    Db      = 0.0
    Dd      = 0.0
    Daphi   = 0.0
    Dbphi   = 0.0
    Ddphi   = 0.0
    Ds      = 0.0
    Dsbeta  = 0.0
    Dsbeta2 = 0.0

    for (ilayer, layer) in enumerate(mat.layers)
      zk   = layer.ct.z1
      zk1  = layer.ct.z2
      E    = layer.ct.E11
      G    = layer.ct.G13
      beta = layer.ct.beta
  
      Da    += E*(zk1  -zk  )
      Db    += E*(zk1^2-zk^2)*(1/2)
      Dd    += E*(zk1^3-zk^3)*(1/3)
  
      Ds      += (zk1 - zk)*G
      Dsbeta  += (zk1 - zk)*G*beta
      Dsbeta2 += (zk1 - zk)*G*beta^2
  
      auxsum = 0.0
      for layer in mat.layers[begin:(ilayer-1)]
          izk   = layer.ct.z1
          izk1  = layer.ct.z2
          ibeta = layer.ct.beta
          auxsum += (izk1-izk)*ibeta
      end
  
      Daphi += E*(zk1^2-zk^2)*(1/2)*beta - E*(zk1-zk)*zk*beta + E*(zk1-zk)*auxsum
  
      Dbphi += E*(zk1^3-zk^3)*(1/3)*beta - E*(zk1^2-zk^2)*(1/2)*zk*beta + E*(zk1^2-zk^2)*(1/2)*auxsum
  
      Ddphi += E*(zk1^3-zk^3)*(1/3)*beta^2 - 2*E*(zk1^2-zk^2)*(1/2)*zk*beta^2 + E*(zk1-zk)*zk^2*beta^2 +
               2*E*(zk1^2-zk^2)*(1/2)*beta*auxsum - (zk1-zk)*2*E*zk*beta*auxsum +
               (zk1-zk)*E*auxsum^2
    end

    arrayABD[1,1] = Da
    arrayABD[1,2] = Db
    arrayABD[1,3] = Daphi

    arrayABD[2,1] = Db
    arrayABD[2,2] = Dd
    arrayABD[2,3] = Dbphi

    arrayABD[3,1] = Daphi
    arrayABD[3,2] = Dbphi
    arrayABD[3,3] = Ddphi

    arrayABD[4,4] = Ds#*(0.932833144)^-1
    arrayABD[4,5] = Dsbeta#*(0.188931942)^-1
    arrayABD[5,4] = Dsbeta#*(0.028141402)^-1
    arrayABD[5,5] = Dsbeta2#*(0.848168271)^-1

    printSparse(arrayABD, 5, 5)

    mat.ct = arrayABD
  end

  function computeMatCT(mat::MaterialPSElastic)
    CTarray = zeros((3,3))

    CT_3D = compute(mat.ct)

    CTarray[1,1] = CT_3D[1,1]
    CTarray[1,2] = CT_3D[1,2]
    CTarray[2,1] = CT_3D[2,1]
    CTarray[2,2] = CT_3D[2,2]
    CTarray[3,3] = CT_3D[3,3]

    return inv(CTarray)
  end

end