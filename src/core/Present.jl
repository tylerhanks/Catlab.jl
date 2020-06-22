""" Finite presentations of a model of a generalized algebraic theory (GAT).

We support two methods for defining models of a GAT: as Julia objects using the
`@instance` macro and as syntactic objects using the `@present` macro.
Instances are useful for casting generic data structures, such as matrices,
abstract tensor systems, and wiring diagrams, in categorical language.
Presentations define small categories by generators and relations and are useful
in applications like knowledge representation.
"""
module Present
export @present, Presentation, Equation, generator, generators, has_generator,
  equations, add_generator!, add_generators!, add_definition!, add_equation!

using Base.Meta: ParseError
using Compat
using Match

using ..Meta, ..Syntax
import ..Syntax: parse_json_sexpr, to_json_sexpr

# Data types
############

const GenExpr = GATExpr{:generator}
const Equation = Pair{<:GATExpr}{<:GATExpr}

mutable struct Presentation{T}
  generators::Vector{GenExpr}
  generator_name_index::Dict{T,Int}
  equations::Vector{Equation}
end
Presentation(T::Type) = Presentation{T}(GenExpr[], Dict{T,Int}(), Equation[])
Presentation() = Presentation(Symbol)

# Presentation
##############

""" Get all generators of a presentation.
"""
function generators(pres::Presentation)::AbstractVector
  pres.generators
end
function generators(pres::Presentation, type::Type)::AbstractVector
  filter(x -> isa(x, type), pres.generators)
end

""" Retrieve generators by name.
"""
function generator(pres::Presentation{T}, name) where T
  pres.generators[pres.generator_name_index[convert(T, name)]]
end
function generators(pres::Presentation, names)::AbstractVector
  [ generator(pres, name) for name in names ]
end

""" Does the presentation contain a generator with the given name?
"""
function has_generator(pres::Presentation{T}, name) where T
  haskey(pres.generator_name_index, convert(T, name))
end

""" Add a generator to a presentation.
"""
function add_generator!(pres::Presentation{T}, expr::GenExpr) where T
  name = first(expr)
  if !isnothing(name)
    name = convert(T, name)
    if haskey(pres.generator_name_index, name)
      error("Name $name already defined in presentation")
    end
    pres.generator_name_index[name] = length(pres.generators) + 1
  end
  push!(pres.generators, expr)
  expr
end

""" Add multiple generators to a presentation.
"""
function add_generators!(pres::Presentation, exprs)
  for expr in exprs
    add_generator!(pres, expr)
  end
end

""" Get all equations of a presentation.
"""
function equations(pres::Presentation)::AbstractVector
  pres.equations
end

""" Add an equation between terms to a presentation.
"""
function add_equation!(pres::Presentation, lhs::GATExpr, rhs::GATExpr)
  push!(pres.equations, lhs => rhs)
end

""" Add a generator defined by an equation.
"""
function add_definition!(pres::Presentation, name::Symbol, rhs::GATExpr)
  generator = Syntax.generator_like(rhs, name)
  add_generator!(pres, generator)
  add_equation!(pres, generator, rhs)
  generator
end

# Presentation macro
####################

""" Define a presentation using a convenient syntax.
"""
macro present(head, body)
  name, syntax_name = @match head begin
    Expr(:call, [name::Symbol, arg::Symbol]) => (name, arg)
    _ => throw(ParseError("Ill-formed presentation header $head"))
  end
  let_expr = Expr(:let, Expr(:block), translate_presentation(syntax_name, body))
  esc(Expr(:(=), name, let_expr))
end

""" Translate a presentation in the DSL to a block of Julia code.
"""
function translate_presentation(syntax_name::Symbol, body::Expr)::Expr
  @assert body.head == :block
  code = Expr(:block)
  append_expr!(code, :(_presentation = $(module_ref(:Presentation))()))
  for expr in strip_lines(body).args
    append_expr!(code, translate_expr(syntax_name, expr))
  end
  append_expr!(code, :(_presentation))
  code
end

""" Translate a single statement in the presentation DSL to Julia code.
"""
function translate_expr(syntax_name::Symbol, expr::Expr)::Expr
  @match expr begin
    Expr(:(::), [name::Symbol, type_expr]) =>
      translate_generator(syntax_name, name, type_expr)
    Expr(:(::), [type_expr]) =>
      translate_generator(syntax_name, nothing, type_expr)
    Expr(:(:=), [name::Symbol, def_expr]) =>
      translate_definition(name, def_expr)
    Expr(:(=), _) => expr
    Expr(:call, [:(==), lhs, rhs]) => translate_equation(lhs, rhs)
    _ => throw(ParseError("Ill-formed presentation statement $expr"))
  end
end

""" Translate declaration of a generator.
"""
function translate_generator(syntax_name::Symbol, name::Union{Symbol,Nothing},
                             type_expr)::Expr
  function rewrite(expr)
    @match expr begin
      Expr(:call, [name::Symbol, args...]) =>
        Expr(:call, GlobalRef(Syntax, :invoke_term), syntax_name,
             QuoteNode(name), map(rewrite, args)...)
      _ => expr
    end
  end
  type_name, args = @match type_expr begin
    name::Symbol => (name, [])
    Expr(:call, [name::Symbol, args...]) => (name, args)
    _ => throw(ParseError("Ill-formed type expression $type_expr"))
  end
  call_expr = Expr(
    :call, module_ref(:add_generator!), :_presentation,
    rewrite(Expr(:call, type_name,
                 isnothing(name) ? :nothing : QuoteNode(name), args...)))
  isnothing(name) ? call_expr : :($(name) = $call_expr)
end

""" Translate definition of generator in terms of other generators.
"""
function translate_definition(name::Symbol, def_expr)::Expr
  call_expr = Expr(:call, module_ref(:add_definition!), :_presentation,
                   QuoteNode(name), def_expr)
  :($name = $call_expr)
end

""" Translate declaration of equality between terms.
"""
function translate_equation(lhs, rhs)::Expr
  Expr(:call, module_ref(:add_equation!), :_presentation, lhs, rhs)
end

module_ref(sym::Symbol) = GlobalRef(Present, sym)

# Serialization
###############

function to_json_sexpr(pres::Presentation, expr::GATExpr)
  to_json_sexpr(expr;
    by_reference = name -> has_generator(pres, name))
end

function parse_json_sexpr(pres::Presentation, syntax_module::Module, sexpr)
  parse_json_sexpr(syntax_module, sexpr;
    parse_reference = name -> generator(pres, name))
end
function parse_json_sexpr(pres::Presentation{Symbol}, syntax_module::Module, sexpr)
  parse_json_sexpr(syntax_module, sexpr;
    symbols = true,
    parse_reference = name -> generator(pres, name))
end

end
