using Catlab
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.CPortGraphs
using  Catlab.CategoricalAlgebra
import Catlab.CategoricalAlgebra: coproduct
import Catlab.WiringDiagrams: oapply

coproduct(xs::Vector{T}) where T<:AbstractVector = (collect ∘ Iterators.flatten)(xs)
#     if length(xs) == 0
#         return eltype(T)[]
#     end
#     y = zeros(eltype(T), sum(length.(xs)))
#     i = 1
#     for x in xs
#         y[i:i+length(x)-1] .= x
#         i += length(x)
#     end
#     return y
# end


struct ParaEucFunc{T}
    update::Function
    readout::Function
    nparams::Int
    state::Vector{T}
end

# (f::ParaEucFunc)(x::Vector, p::Vector) = begin
update!(f::ParaEucFunc, p::Vector, h) = begin
    @assert length(p) == f.nparams
    f.state .+= h*f.update(f.state, p)
end

readout(f::ParaEucFunc) = f.readout(f.state)
readout(f::ParaEucFunc, i) = readout(f)[i]


function oapply(d, xs::Vector)
    function u(h::Real)
        readouts = coproduct(map(f->f.readout(f.state), xs))
        map(parts(d, :B)) do b
            @show ports = incident(d, b, :box)
            @show wires = incident(d, ports, :tgt) |> coproduct
            @show feeders = d[wires, :src]
            @show params = readouts[feeders]
            update!(xs[b], params, h)
        end
    end
    function r()
        map(parts(d, :OP)) do op
            p = d[op, :con]
            b = d[p, :box]
            readout(xs[b])
        end |> Iterators.flatten |> collect
    end
    return u,r
end

barbell(k::Int) = begin
  g = OpenCPortGraph()
  add_parts!(g, :B, 2)
  add_parts!(g, :P, 2k; box=[fill(1,k); fill(2,k)])
  add_parts!(g, :W, k; src=1:k, tgt=k+1:2k)
  add_parts!(g, :W, k; tgt=1:k, src=k+1:2k)
  return g
end

d = barbell(2)
xs = [
    ParaEucFunc((x,p)->[p[1]*x[1], p[2]*x[2]], x->x, 2, [1,1.0]),
    ParaEucFunc((x,p)->[p[1]*x[1], -p[2]*x[2]], x->x, 2, [1,1.0]),
]

f, r = oapply(d, xs)
@show r()
@show f(0.1)
@show f(0.1)