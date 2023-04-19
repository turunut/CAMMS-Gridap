#=
Exemple LINEAR ELASTIC generant malla / Dirichelt

Exemple Linear elastic en que generem un model cuadrat de forma
procedural, el mallem i definim condicions de Dirichlet imposant
desplacaments als seus extrems esquerra i dret.
=#
push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/Materials")

using Gridap
#using GridapGmsh
using GridapGiD
using Gridap.Geometry
using modCT
using modModel
using modSubroutines
using Gridap.TensorValues
using Gridap.Arrays

prblName = "Reissner3D_2x2"
projFldr = pwd()

#n = 2
#domain = (0,10,0,10)
#partition = (n,n)
#model = CartesianDiscreteModel(domain,partition)

model = GiDDiscreteModel( projFldr*"/models/"*prblName*"/"*prblName )

writevtk(model,"model")

order = 1
degree = 2*order

#writevtk(model,"model")

labels = get_face_labeling(model)

#add_tag_from_tags!(labels,"diri_l",[7, 1, 3]) # left  edge
#add_tag_from_tags!(labels,"diri_r",[4]) # right edge


#--------------------------------------------------


MT = Reissner(5.0, -5.0)

CT1 = CT_Isotrop(72000, 0.3)
ct1 = modModel.computeCT(MT, CT1)


#--------------------------------------------------


reffe3 = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
reffe2 = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
Vv = TestFESpace(model,reffe3;
                 conformity=:H1,
                 dirichlet_tags=["left","leftp","right"],
                 dirichlet_masks=[(true,true,true),(true,true,true),(false,false,true)])
Vt = TestFESpace(model,reffe2;
                 conformity=:H1,
                 dirichlet_tags=["left","left"],
                 dirichlet_masks=[(true,true),(true,true)])
V = MultiFieldFESpace([Vv,Vt])

g0_3(x) = VectorValue( 0.0, 0.0, 0.0)
g1_3(x) = VectorValue( 0.0, 0.0,-1.0)
g0_2(x) = VectorValue( 0.0, 0.0)
g1_2(x) = VectorValue(-1.0, 0.0)

Uv = TrialFESpace(Vv, [g0_3,g0_3,g1_3])
Ut = TrialFESpace(Vt, [g0_2,g0_2])
U = MultiFieldFESpace([Uv,Ut])


# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 2 # surfaces
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1) # Posem un sobre els altres [[A₁,B₁,D₁,S₁],
                #                            [A₂,B₂,D₂,S₂]]

CTf = get_CT_CellField(MT, CTs, tags, Ω)


#--------------------------------------------------

function symmetrize(t)
  return SymTensorValue(get_array(0.5*(transpose(t) + t)))
end

function my_tangent(m)
  ê₁ = m(VectorValue(1.0, 0.0)) - m(VectorValue(0.0, 0.0)); ê₁ = ê₁/norm(ê₁)
  ê₂ = m(VectorValue(0.0, 1.0)) - m(VectorValue(0.0, 0.0)); ê₂ = ê₂/norm(ê₂)
  inplane  = TensorValue([ ê₁[1] ê₂[1];   #ê₃[1];
                           ê₁[2] ê₂[2];   #ê₃[2];
                           ê₁[3] ê₂[3] ]) #ê₃[3] ])
  return inplane
end
function my_tangentθ(m)
ê₁ = m(VectorValue(1.0, 0.0)) - m(VectorValue(0.0, 0.0)); ê₁ = ê₁/norm(ê₁)
ê₂ = m(VectorValue(0.0, 1.0)) - m(VectorValue(0.0, 0.0)); ê₂ = ê₂/norm(ê₂)
inplaneθ  = TensorValue([ ê₁[1] ê₂[1];
                          ê₁[2] ê₂[2] ])
return inplaneθ
end
function my_normal(m)
  ê₁ = m(VectorValue(1.0, 0.0)) - m(VectorValue(0.0, 0.0)); ê₁ = ê₁/norm(ê₁)
  ê₂ = m(VectorValue(0.0, 1.0)) - m(VectorValue(0.0, 0.0)); ê₂ = ê₂/norm(ê₂)
  outplane = cross(ê₁,ê₂)
  return outplane
end

function get_tan_ref_sys(Ω::Triangulation{2})
    cmaps = get_cell_map(Ω)
    return CellField(lazy_map(my_tangent,cmaps),Ω)
end
function get_tan_ref_sysθ(Ω::Triangulation{2})
    cmaps = get_cell_map(Ω)
    return CellField(lazy_map(my_tangentθ,cmaps),Ω)
end
function get_nor_ref_sys(Ω::Triangulation{2})
    cmaps = get_cell_map(Ω)
    return CellField(lazy_map(my_normal,cmaps),Ω)
end

tf  = get_tan_ref_sys(Ω)
tfθ = get_tan_ref_sysθ(Ω)
nf  = get_nor_ref_sys(Ω)

function getₑ(x,ê)
  return symmetrize( ê'⋅x )
end
#getₑ(x,ê) = symmetrize( ê'⋅x )

getₙ(x) = VectorValue( sign(sum(x))*norm(x) )
getᵥ(x) = VectorValue( x )

function ∂ₙ(u,ê,d)
  return getₑ ∘ ( ∇(u)'⋅ê, d )
end

#∂ₙ(u,ê,d) = getₑ ∘ ( ∇(u)⋅ê, d )
∂ᵥ(θ,ê)   = getᵥ ∘ ( ∇(θ)⋅ê )

q(x) = VectorValue(-1.0)


function getₑ₂(x,ê)
  return ê'⋅x
end
function ∂ₙ₂(u,ê,d)
  return getₑ₂ ∘ ( ∇(u)'⋅ê, d )
end



A((u,θ),(v,t)) = ∫( ∂ₙ(v,tf,tf)⊙σ(CTf[1],∂ₙ(u,tf,tf)) + ∂ₙ(t,tf,tfθ)⊙σ(CTf[2],∂ₙ(θ,tf,tfθ)) )*dΩ
S((u,θ),(v,t)) = ∫( γ(MT,∂ₙ₂(v,nf,tf),t) ⊙ σₑ(CTf[4], γ(MT,∂ₙ₂(u,nf,tf),θ)) )*dΩ 

UU = get_trial_fe_basis(U)
VV = get_fe_basis(V)
contrA = A(UU,VV)
elementA = first(contrA.dict).second[1]
contrS = S(UU,VV)
elementS = first(contrS.dict).second[1]





a((u,θ),(v,t)) = ∫( ∂ₙ(v,tf,tf)⊙σ(CTf[1],∂ₙ(u,tf,tf)) + ∂ₙ(t,tf,tfθ)⊙σ(CTf[2],∂ₙ(θ,tf,tfθ)) )*dΩ + # Axial         + Axial/Bending
                 ∫( ∂ₙ(v,tf,tf)⊙σ(CTf[2],∂ₙ(u,tf,tf)) + ∂ₙ(t,tf,tfθ)⊙σ(CTf[3],∂ₙ(θ,tf,tfθ)) )*dΩ + # Bending/Axial + Bending
                 ∫( γ(MT,∂ₙ₂(v,nf,tf),t) ⊙ σₑ(CTf[4], γ(MT,∂ₙ₂(u,nf,tf),θ)) )*dΩ 

l((v,t)) = 0.0
#l((v,t)) = ∫(w⋅q)*dΩ


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)
vh = uh.single_fe_functions[1]
th = uh.single_fe_functions[2]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u" =>vh,
                     "θ" =>th])
                     #"∂ω"=>∂(wh)]) # "∂ω"=>ε(wh)⋅VectorValue(1.0),
                     #"γ" =>γ(MT,wh,th)])


#--------------------------------------------------


A((u,θ),(v,t)) = ∫( ∂ₙ(v,tf,tf)⊙σ(CTf[1],∂ₙ(u,tf,tf)) + ∂ₙ(t,tf,tfθ)⊙σ(CTf[2],∂ₙ(θ,tf,tfθ)) )*dΩ
D((u,θ),(v,t)) = ∫( ∂ₙ(v,tf,tf)⊙σ(CTf[2],∂ₙ(u,tf,tf)) + ∂ₙ(t,tf,tfθ)⊙σ(CTf[3],∂ₙ(θ,tf,tfθ)) )*dΩ
S((u,θ),(v,t)) = ∫( γ(MT,∂ₙ₂(v,nf,tf),t) ⊙ σₑ(CTf[4], γ(MT,∂ₙ₂(u,nf,tf),θ)) )*dΩ 

UU = get_trial_fe_basis(U)
VV = get_fe_basis(V)
contrA = A(UU,VV)
elementA = first(contrA.dict).second[1]
contrD = D(UU,VV)
elementD = first(contrD.dict).second[1]
contrS = S(UU,VV)
elementS = first(contrS.dict).second[1]

contr = a(UU,VV)
element = first(contr.dict).second[1]

