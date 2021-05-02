module OGraphs
using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
import Catlab.CategoricalAlgebra.CSets: migrate!
using Catlab.Present
using Catlab.Graphs
using Catlab.Graphs.BasicGraphs
import Catlab.Graphs.BasicGraphs: TheoryGraph


@present TheoryOGraph <: TheoryGraph begin
    IP::Ob
    OP::Ob
    input::Hom(IP, V)
    output::Hom(OP, V)
end

const OGraph = ACSetType(TheoryOGraph)

# construct an open graph using pullback migration from a graph
function OGraph(g::Graph)
    G = OGraph()
    migrate!(G, g, Dict(:V=>:V, :E=>:E, :IP=>:V, :OP=>:V),
       Dict(:src=>:src, :tgt=>:tgt,
           :input=>id(TheoryOGraph[:V]),
           :output=>id(TheoryOGraph[:V])))
    return G
end

Graph(G::OGraph) = begin
    g = Graph()
    migrate!(g, G, Dict(:V=>:V, :E=>:E), Dict(:src=>:src, :tgt=>:tgt))
    return g
end

function ocompose(g::OGraph, args::Vector)
    u = coproduct(args)
    @show u
    # @show apex(u)
    # @show cocone(u)
    @show legs(u)[1]
    @show legs(u)[2]
end

g = Graph()
add_vertices!(g, 4)
add_edges!(g, [1,2,3,4], [2,3,4,1])
G = @show OGraph(g)

bar = Graph()
add_vertices!(bar, 2)
add_edge!(bar, 1,2)
Bar = OGraph(bar)
ocompose(G, [Bar,Bar, Bar, Bar])

end