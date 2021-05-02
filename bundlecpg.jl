using Test
using Catlab
using Catlab.Graphs
using Catlab.WiringDiagrams
import Catlab.WiringDiagrams: oapply
using Catlab.WiringDiagrams.CPortGraphs
import Catlab.WiringDiagrams.CPortGraphs: CPortGraph
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.CSets
import Catlab.CategoricalAlgebra.CSets: migrate!
using Catlab.Theories

# function CPortGraph(g::Graph)
#     cg = CPortGraph()
#     copy_parts!(cg, g, 
#       W=parts(g, :E), P = parts(g,:V), B=parts(g, :B))
#     return cg
# end

@present ThBundledCPG <: ThOpenCPortGraph begin
    Bundle::Ob
    bun::Hom(OP, Bundle)
end

const BundledCPG = ACSetType(ThBundledCPG)

migrate!(b::BundledCPG, g::OpenCPortGraph) = migrate!(b,g,
    Dict(:B=>:B, :P=>:P, :W=>:W, :OP=>:OP, :Bundle=>:OP),
    Dict(:src=>:src, :tgt=>:tgt, :box=>:box,
         :con=>:con, :bun=>id(ThOpenCPortGraph[:OP]))
)

function BundledCPG(g::OpenCPortGraph)
    bg = BundledCPG()
    copy_parts!(bg, g, 
      W=parts(g,:W), P = parts(g,:P), B=parts(g, :B), OP=parts(g,:OP))
    return bg
end

function oapply(g::BundledCPG, xs::Vector)
    u = coproduct(xs)
    xsum=apex(u)
    for e in parts(g, :W)
        s = g[e,:src]
        t = g[e,:tgt]
        sbox = g[s, :box]
        tbox = g[t, :box]

        localport_src = findfirst(s .== incident(g, sbox, :box))
        localport_tgt = findfirst(t .== incident(g, tbox, :box))

        println("$s:$sbox→$t:$tbox")
        sbun = incident(xs[sbox], localport_src, :bun)
        println("sbun: $sbun")
        tbun = incident(xs[tbox], localport_tgt, :bun)
        println("tbun: $tbun")
        for thread in zip(sbun, tbun)
            ι_srcport = legs(u.cocone)[sbox][:P](thread[1])
            ι_tgtport = legs(u.cocone)[tbox][:P](thread[2])
            add_part!(xsum, :W; src=ι_srcport, tgt=ι_tgtport) 
        end
    end
    rem_parts!(xsum, :OP, parts(xsum, :OP))
    rem_parts!(xsum, :Bundle, parts(xsum, :Bundle))
    add_parts!(xsum, :Bundle, nparts(g, :Bundle))
    for op in parts(g, :OP)
        i = g[op, [:con, :box]]
        localport = findfirst(op .== incident(g, i, :box))
        newop = legs(u.cocone)[i][:P](incident(xs[i], localport, :bun))
        add_parts!(xsum, :OP, length(newop), con=newop, bun=op)
    end
    return xsum
end


EdgeGraph() = begin
g = Graph()
add_parts!(g, :V, 2)
add_parts!(g, :E, 1, src=[1], tgt=[2])
return g
end

e₁ = EdgeGraph()
oe₁ = migrate!(OpenCPortGraph(), (migrate!(CPortGraph(), e₁)))
@test subpart(oe₁, :con) == 1:2
mboe₁ = migrate!(BundledCPG(), oe₁)
@test subpart(mboe₁, :bun) == 1:2

boe₁ = BundledCPG(oe₁)
add_parts!(boe₁, :Bundle, 1)
set_subpart!(boe₁, :bun, [1,1])
boe₁

f = OpenCPortGraph(migrate!(CPortGraph(), e₁))
f = BundledCPG(f)
fee = oapply(f, [boe₁, boe₁])
@test fee[:, :src] == [1,3,1,2]
@test fee[:, :tgt] == [2,4,3,4]

add_parts!(f, :Bundle, 2)
add_parts!(f, :OP, 2, con=[1,2], bun=[1,2])
fee′ = oapply(f, [boe₁, boe₁])
@test fee′[:, :con] == 1:4
@test fee′[:, :bun] == [1,1,2,2]