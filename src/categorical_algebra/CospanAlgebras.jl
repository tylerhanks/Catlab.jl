module CospanAlgebras

export CospanAlgebra, Open, DecoratedCospanAlgebra, hom_map, laxator, data, portmap

# Interface for Abstract Cospan Algebras
########################################

"""     AbstractCospanAlgebra{T}

A Cospan(FinSet)-algebra on a type T.
"""
abstract type AbstractCospanAlgebra{T} end

# Cospan algebras need to implement the following methods.

"""    ob_map(a::AbstractCospanAlgebra{T}, N::FinSet, X::T)::Bool where T  

Predicate to check whether X ∈ a(N).
"""
function ob_map(a::AbstractCospanAlgebra{T}, N::FinSet, X::T)::Bool where T
    error("Object map not implemented for $(typeof(a))")
end

"""     hom_map(a::AbstractCospanAlgebra{T}, ϕ::Cospan, X::T)::T where T

Apply the morphism map a(ϕ)(X). It is up to the user to verify that ob_map(a, dom(ϕ), X) == true.
"""
function hom_map(a::AbstractCospanAlgebra{T}, ϕ::Cospan, X::T)::T where T
    error("Morphism map not implemented for $(typeof(a))")
end

"""     laxator(a::AbstractCospanAlgebra{T}, Xs::Vector{T})::T where T

Apply the laxator φ(Xs)
"""
# TODO: should probably add the ability to type check laxator applications.
function laxator(a::AbstractCospanAlgebra{T}, Xs::Vector{T})::T where T
    error("Laxator not implemented for $(typeof(a))")
end

"""     oapply(a::AbstractCospanAlgebra{T}, ϕ::Cospan, Xs::Vector{T})::T where T

Apply the whole operad algebra.
"""
function oapply(a::AbstractCospanAlgebra{T}, ϕ::Cospan, Xs::Vector{T})::T where T
    X = laxator(a, Xs)
    ob_map(a, dom(left(ϕ)), X) || error("Failed to laxate to a correctly typed object")
    return hom_map(a, ϕ, X)
end

# Cospan Algebras from Decorating Functors
##########################################

abstract type FinSetAlgebra{T} end


function hom_map(::FinSetAlgebra{T}, ϕ::FinFunction, X::T)::T where T
    error("Morphism map not implemented.")
end

function laxator(::FinSetAlgebra{T}, Xs::Vector{T})::T where T
    error("Laxator not implemented.")
end

function oapply(A::FinSetAlgebra{T}, ϕ::FinFunction, Xs::Vector{T})::T where T
    return hom_map(A, ϕ, laxator(A, Xs))
end



struct Open{T}
    S::FinSet
    o::T
    m::FinFunction    
    Open{T}(S, o, m) where T = 
        S != codom(m) || dom(o) != S ? error("Invalid portmap.") : new(S, o, m)
end

data(obj::Open{T}) where T = obj.o
portmap(obj::Open{T}) where T = obj.m
dom(obj::Open{T}) where T = dom(obj.m)


#abstract type DecoratedCospanAlgebra{T, a<:FinSetAlgebra{T}} <: AbstractCospanAlgebra{T} end
struct DecoratedCospanAlgebra{Open{T}} <: AbstractCospanAlgebra{Open{T}}
    decorating_algebra::FinSetAlgebra{T}
end

function hom_map(d::DecoratedCospanAlgebra{Open{T}}, ϕ::Cospan, X::Open{T})::Open{T} where T
    l = left(ϕ)
    r = right(ϕ)
    p = pushout(X.m, l)
    pL = legs(p)[1]
    pR = legs(p)[2]
    return Open{T}(apex(p), hom_map(d.decorating_algebra, pL, X.o), compose(r,pR))
end

function laxator(::CospanAlgebra{Open{T}}, A::FinSetAlgebra{T}, Xs::Vector{Open{T}})::Open{T} where T
    S = coproduct([X.S for X in Xs])
    inclusions(i::Int) = legs(S)[i]
    m = copair([compose(Xs[i].m, inclusions(i)) for i in 1:length(Xs)])
    o = laxator(A, [X.o for X in Xs])
    return Open{T}(apex(S), o, m)
end

function oapply(CA::CospanAlgebra{Open{T}}, FA::FinSetAlgebra{T}, ϕ::Cospan, Xs::Vector{Open{T}})::Open{T} where T
    return hom_map(CA, FA, ϕ, laxator(CA, FA, Xs))
end

function uwd_to_cospan(d::AbstractUWD)
    # Build the left leg
    left_dom = vcat([length(ports(d, i)) for i in boxes(d)])
    left_codom = njunctions(d)

    #println(cp_dom)
    ports_to_junctions = FinFunction[]
    total_portmap = subpart(d, :junction)

    for box in ports.([d], boxes(d))
        push!(ports_to_junctions, FinFunction([total_portmap[p] for p in box], length(box), left_codom))
    end
    #println(ports_to_junctions)
    #cp = CompositionPattern(cp_dom, cp_codom, ports_to_junctions)

    left = copair(ports_to_junctions)
    right = FinFunction(subpart(d, :outer_junction), left_codom)
    
    return Cospan(left, right)  
end

function oapply(CA::CospanAlgebra{Open{T}}, FA::FinSetAlgebra{T}, d::AbstractUWD, Xs::Vector{Open{T}})::Open{T} where T
    return oapply(CA, FA, uwd_to_cospan(d), Xs)
end


end