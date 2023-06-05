


function get_line_distribution(num_divisions::Int64, intrf::inter3D, reffe_line, U_FE, active_DOF)
  _get_y(x) = VectorValue(x[2])
  function π_Λe_Γc(f::CellField, Γc::Triangulation)
    _data = CellData.get_data(f)
    _cellmap = Fill(Broadcasting(_get_y),length(_data))
    data = lazy_map(∘,_data,_cellmap)
    return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
  end
  
  domainC = (0, 1); partitionC = (num_divisions)
  modelC  = CartesianDiscreteModel(domainC, partitionC)
  ΩC = Triangulation(modelC)
  VC = FESpace(ΩC,reffe_line,conformity=:H1); UC = TrialFESpace(VC)
  # Definim la funcio a partir del DOFs
  xC = zero_free_values(UC);
  if active_DOF != 0; xC[active_DOF] = 1.0; end
  uC = FEFunction(UC,xC)
  # -----------------
  uC_intrp = Interpolable(uC)
  
  ue = interpolate(uC_intrp,U_FE)
  ue_c = π_Λe_Γc(ue,intrf.Γc)
  return ue_c
end



function define(num_divisions::Int64, intrf::inter3D, reffe_line, U_FE, active_DOF)
  _get_y(x) = VectorValue(x[2])
  function π_Λe_Γc(f::CellField, Γc::Triangulation)
    _data = CellData.get_data(f)
    _cellmap = Fill(Broadcasting(_get_y),length(_data))
    data = lazy_map(∘,_data,_cellmap)
    return CellData.similar_cell_field(f,data,Γc,CellData.DomainStyle(f))
  end
  
  domainC = (0, 1); partitionC = (num_divisions)
  modelC  = CartesianDiscreteModel(domainC, partitionC)
  ΩC = Triangulation(modelC)
  VC = FESpace(ΩC,reffe_line,conformity=:H1); UC = TrialFESpace(VC)
  # Definim la funcio a partir del DOFs
  xC = zero_free_values(UC);
  if active_DOF != 0; xC[active_DOF] = 1.0; end
  uC = FEFunction(UC,xC)
  # -----------------
  uC_intrp = Interpolable(uC)
  
  ue = interpolate(uC_intrp,U_FE)
  ue_c = π_Λe_Γc(ue,intrf.Γc)
  return ue_c
end








function get_test_trial_spaces(intrf::inter3D)
  dofs = get_dofs(intrf)
  reffe = ReferenceFE(lagrangian,VectorValue{dofs,Float64},0)
  Vλ = FESpace(intrf.Γc,reffe,conformity=:L2)
  Uλ = TrialFESpace(Vλ)
  return Vλ, Uλ
end

function get_line_test_trial_spaces(intrf::inter3D, order)
  dofs = get_dofs(intrf)
  reffeλ = ReferenceFE(lagrangian,VectorValue{dofs,Float64},order)
  Vλ = FESpace(intrf.Λe,reffeλ,conformity=:H1)
  Uλ = TrialFESpace(Vλ)
  return Vλ, Uλ, reffeλ
end

function contribute_vector(intrf::inter3D, V_basis, V_ind::Int64, f)
  μ = V_basis[V_ind]
  tr_Γf(λ) = change_domain(λ,intrf.Γf,DomainStyle(λ))
  return ∫( tr_Γf(μ)⋅f )*intrf.dΓf
end
