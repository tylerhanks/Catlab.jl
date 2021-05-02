# cpg_uwd_migration.jl

using Catlab.WiringDiagrams
const UWD = UndirectedWiringDiagram

function migrate!(uwd::UWD, pg::OpenCPortGraph) 
    migrate!(uwd, pg,
        Dict(:Box=>:B, :Port=>:W, :Junction=>:P, :OuterPort=>:OP),
        Dict(:box=>[:src, :box],
             :outer_junction=>[:con],
             :junction=>[:src],
        )
    )
end

migrate!(UntypedUWD(), grid(2,2))
to_graphviz(migrate!(UntypedUWD(), grid(3,3)), port_labels=true, junction_labels=false)
