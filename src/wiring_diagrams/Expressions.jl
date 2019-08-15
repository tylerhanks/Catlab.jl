""" Convert wiring diagrams to and from syntactic expressions.

Wiring diagrams are a graphical syntax for morphisms in a monoidal category. As
mathematical objects, they are intermediate between the morphisms themselves and
expressions in the textual syntax: a single morphism may correspond to many
wiring diagrams, and a single wiring diagram may correspond to many syntactic
expressions.

It is trivial to convert a morphism expression to a wiring diagram, but not to
go the other way around.
"""
module WiringDiagramExpressions
export to_hom_expr, to_wiring_diagram

using LightGraphs

using ...Syntax
using ...Doctrines: id, compose, otimes, munit
using ..WiringDiagramCore, ..WiringLayers
using ..WiringDiagramAlgorithms: crossing_minimization_by_sort

# Expression -> Diagram
#######################

""" Convert a morphism expression into a wiring diagram.

The morphism expression should belong to the doctrine of symmetric monoidal
categories, possibly with diagonals and codiagonals. Thus, the doctrines of
cartesian, cocartesian, and biproduct categories are supported.
"""
function to_wiring_diagram(expr::GATExpr)
  functor((Ports, WiringDiagram), expr;
    terms = Dict(
      :Ob => expr -> Ports([first(expr)]),
      :Hom => expr -> WiringDiagram(expr),
    )
  )
end

# Diagram -> Expression
#######################

""" Convert a wiring diagram into a morphism expression.
"""
function to_hom_expr(Ob::Type, Hom::Type, d::WiringDiagram)
  box_to_expr(box::Box) = to_hom_expr(Ob, Type, box)
  
  d = copy(d)
  n = nboxes(d)
  while true
    parallel_reduction!(Ob, Hom, d)
    series_reduction!(Ob, Hom, d)
    
    # Repeat until the box count stops decreasing.
    nboxes(d) >= n && break
    n = nboxes(d)
  end
  
  @assert nboxes(d) == 1
  v = first(box_ids(d))
  f = to_hom_expr(Ob, Hom, box(d,v))
  in_expr = hom_expr_between(Ob, d, input_id(d), v)
  out_expr = hom_expr_between(Ob, d, v, output_id(d))
  foldl(compose_simplify_id, [in_expr, f, out_expr])
end
function to_hom_expr(Syntax::Module, d::WiringDiagram)
  to_hom_expr(Syntax.Ob, Syntax.Hom, d)
end

""" Perform all possible parallel reductions in wiring diagram.
"""
function parallel_reduction!(Ob::Type, Hom::Type, d::WiringDiagram)
  parallel = collect(find_parallel(graph(d)))
  groups = map(parallel) do ((source, target), vs)
    exprs = Hom[ to_hom_expr(Ob, Hom, box(d,v)) for v in vs ]
    # FIXME: Do two-sided crossing minimization, not one-sided.
    σ = crossing_minimization_by_sort(d, [source], vs)
    vs[σ] => foldl(otimes_simplify_id, exprs[σ])
  end
  if !isempty(groups)
    encapsulate!(d, map(first, groups),
                 discard_boxes=true, values=map(last, groups))
  end
end

""" Perform all possible series reductions in wiring diagram.
"""
function series_reduction!(Ob::Type, Hom::Type, d::WiringDiagram)
  box_to_expr(v::Int) = to_hom_expr(Ob, Hom, box(d,v))
  
  series = find_series(graph(d), source=input_id(d), sink=output_id(d))
  composites = map(series) do vs
    exprs = Hom[ box_to_expr(vs[1]) ]
    for i in 2:length(vs)
      layer_expr = hom_expr_between(Ob, d, vs[i-1], vs[i])
      append!(exprs, [layer_expr, box_to_expr(vs[i])])
    end
    foldl(compose_simplify_id, exprs)
  end
  
  if !isempty(series)
    encapsulate!(d, series, discard_boxes=true, values=composites)
  end
end

function hom_expr_between(Ob::Type, diagram::WiringDiagram, v1::Int, v2::Int)
  layer = wiring_layer_between(diagram, v1, v2)
  inputs = ports_to_obs(Ob, output_ports(diagram, v1))
  outputs = ports_to_obs(Ob, input_ports(diagram, v2))
  to_hom_expr(layer, inputs, outputs)
end

function to_hom_expr(Ob::Type, Hom::Type, box::Box)
  if box.value isa Hom
    box.value
  else
    dom = otimes(ports_to_obs(Ob, input_ports(box)))
    codom = otimes(ports_to_obs(Ob, output_ports(box)))
    Hom(box.value, dom, codom)
  end
end

coerce_ob(Ob::Type, x) = x isa Ob ? x : Ob(Ob,x)
ports_to_obs(Ob::Type, ports) = Ob[ coerce_ob(Ob, port) for port in ports ]

function compose_simplify_id(f::GATExpr, g::GATExpr)
  if head(f) == :id; g
  elseif head(g) == :id; f
  else compose(f,g) end
end
function otimes_simplify_id(f::GATExpr, g::GATExpr)
  if head(f) == :id && head(g) == :id; id(otimes(dom(f),dom(g)))
  else otimes(f,g) end
end

# Layer -> Expression
#####################

""" Convert a wiring layer into a morphism expression.

The `inputs` and `outputs` are corresponding vectors of object expressions.
"""
function to_hom_expr(layer::WiringLayer, inputs::Vector, outputs::Vector)
  σ = to_permutation(layer)
  @assert !isnothing(σ) "Conversion of non-permutation not implemented"
  @assert inputs == outputs
  permutation_to_hom_expr(σ, inputs)
end

""" Convert a typed permutation into a morphism expression.
"""
function permutation_to_hom_expr(σ::Vector{Int}, objects::Vector)
  @assert all(σ[i] == i for i in eachindex(σ)) "Conversion of non-identity not implemented"
  id(otimes(objects))
end

""" Convert a wiring layer into a permutation, if it is one.

Otherwise, return nothing.
"""
function to_permutation(layer::WiringLayer)::Union{Nothing,Vector{Int}}
  nin, nout = layer.ninputs, layer.noutputs
  if nin != nout; return nothing end
  σ = zeros(Int, nin)
  for i in 1:nin
    wires = out_wires(layer, i)
    if length(wires) == 1
      σ[i] = last(first(wires))
    else
      return nothing
    end
  end
  return σ
end

# Graph operations
##################

""" Find parallel compositions in a graph.
"""
function find_parallel(g::DiGraph)::Dict{Pair{Int,Int},Vector{Int}}
  parallel = Dict{Pair{Int,Int},Vector{Int}}()
  for v in 1:nv(g)
    if length(inneighbors(g,v)) == 1 && length(outneighbors(g,v)) == 1
      src, tgt = first(inneighbors(g,v)), first(outneighbors(g,v))
      push!(get!(parallel, src => tgt, Int[]), v)
    end
  end
  filter(pair -> length(last(pair)) > 1, parallel)
end

""" Find series compositions in a graph.
"""
function find_series(g::DiGraph; source=nothing, sink=nothing)::Vector{Vector{Int}}
  reduced = DiGraph(nv(g))
  for edge in edges(g)
    if (length(outneighbors(g,src(edge))) == 1 &&
        length(inneighbors(g,dst(edge))) == 1 &&
        src(edge) != source && dst(edge) != sink)
      add_edge!(reduced, edge)
    end
  end
  series = Vector{Int}[]
  for component in weakly_connected_components(reduced)
    if length(component) > 1
      sub, vmap = induced_subgraph(reduced, component)
      push!(series, Int[ vmap[v] for v in topological_sort_by_dfs(sub) ])
    end
  end
  series
end

""" Transitive reduction of a directed acyclic graph.

Algorithm described at CS.SE: https://cs.stackexchange.com/a/29133
"""
function transitive_reduction!(g::DiGraph)::DiGraph
  for u in 1:nv(g)
    for v in outneighbors(g, u)
      for w in findall(!iszero, LightGraphs.dfs_parents(g, v))
        if w != v
          rem_edge!(g, u, w)
        end
      end
    end
  end
  return g
end

end
