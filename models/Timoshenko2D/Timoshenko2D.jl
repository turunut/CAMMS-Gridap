#=
Exemple LINEAR ELASTIC generant malla / Dirichelt

Exemple Linear elastic en que generem un model cuadrat de forma
procedural, el mallem i definim condicions de Dirichlet imposant
desplacaments als seus extrems esquerra i dret.
=#
push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/Materials")

using Gridap
using GridapGmsh
using Gridap.Geometry
using modCT
using modModel
using modSubroutines
using Gridap.TensorValues
using Gridap.Arrays

prblName = "Timoshenko"
projFldr = pwd()

n = 10
domain = (0,2000)
partition = (n)
model = CartesianDiscreteModel(domain,partition)

order = 1

#writevtk(model,"model")

labels = get_face_labeling(model)


#--------------------------------------------------


MT = Timoshenko(1.0, 50.0, -50.0)

CT1 = CT_Isotrop(72400, 0.3)
ct1 = modModel.computeCT(MT, CT1)

CT2 = CT_Isotrop(7240, 0.3)
ct2 = modModel.computeCT(MT, CT2)


#--------------------------------------------------


reffe2 = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
reffe1 = ReferenceFE(lagrangian,VectorValue{1,Float64},order)
Vv = TestFESpace(model,reffe2;
                 conformity=:H1,
                 dirichlet_tags=["tag_1", "tag_2"],
                 dirichlet_masks=[(true,true), (true,true)])
Vt = TestFESpace(model,reffe1;
                 conformity=:H1,
                 dirichlet_tags=["tag_1", "tag_2"],
                 dirichlet_masks=[(true), (false)])
V = MultiFieldFESpace([Vv,Vt])

g0_1(x) = VectorValue(  0.0)
g1_1(x) = VectorValue(-10.0)
g0_2(x) = VectorValue(0.0, 0.0)
g1_2(x) = VectorValue(0.0,-10.0)

Uv = TrialFESpace(Vv, [g0_2,g0_2])
Ut = TrialFESpace(Vt, [g0_1,g1_1])
U = MultiFieldFESpace([Uv,Ut])


#--------------------------------------------------


degree = 1*order
# Definim l'integration mesh
Ω = Triangulation(model)
# Contruim el l'espai de mesura de Lebesgues de ordre "degree"
dΩ = Measure(Ω,degree)


#--------------------------------------------------


dimens  = 1 # linies
matFlag = []

tags = get_tags(matFlag, labels, dimens)

CTs = hcat(ct1, ct2) # Posem un sobre els altres [[A₁,B₁,D₁,S₁],
                     #                            [A₂,B₂,D₂,S₂]]

CTf = get_CT_CellField(MT, CTs, tags, Ω)


#--------------------------------------------------


function my_tangent(m)
    t = m(VectorValue(1.0)) - m(VectorValue(0.0))
    return t/norm(t)
end

function get_tangent_vector(Ω::Triangulation{1})
    cmaps = get_cell_map(Ω)
    CellField(lazy_map(my_tangent,cmaps),Ω)
end

tangent_cf = get_tangent_vector(Ω)

function g2l(u::CellField)
    return
end

a((u,θ),(v,t)) = ∫( ∂(g2l(v))⊙σ(CTf[1],∂(g2l(u))) + ∂(t)⊙σ(CTf[2],∂(θ)) )*dΩ + # Axial         + Axial/Bending
                 ∫( ∂(g2l(v))⊙σ(CTf[2],∂(g2l(u))) + ∂(t)⊙σ(CTf[3],∂(θ)) )*dΩ + # Bending/Axial + Bending
                 ∫( γ(MT,g2l(v),t) ⊙ σ(CTf[4], γ(MT,g2l(u),θ)) )*dΩ        # Shear

l((v,t)) = 0


#--------------------------------------------------


op = AffineFEOperator(a,l,U,V)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(op)
vh = uh.single_fe_functions[1]
wh = uh.single_fe_functions[2]
th = uh.single_fe_functions[3]

writevtk(Ω,"models/"*prblName*"/"*prblName,
         cellfields=["u" =>vh,
                     "ω" =>wh,
                     "θ" =>th,
                     "∂ω"=>∂(wh), # "∂ω"=>ε(wh)⋅VectorValue(1.0),
                     "γ" =>γ(MT,wh,th)])

