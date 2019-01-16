""" Serialize abstract wiring diagrams as JSON.

JSON data formats are convenient when programming for the web. Unfortunately, no
standard for JSON graph formats has gained any kind of widespread adoption. We
adopt a format compatible with that used by the KEILER project and its successor
ELK (Eclipse Layout Kernel). This format is roughly feature compatible with
GraphML, supporting nested graphs and ports. It also supports layout information
like node position and size.

References:

- KEILER's JSON graph format:
  https://rtsys.informatik.uni-kiel.de/confluence/display/KIELER/JSON+Graph+Format
- ELK's JSON graph format:
  https://www.eclipse.org/elk/documentation/tooldevelopers/graphdatastructure/jsonformat.html
"""
module JSONWiringDiagrams
export read_json_graph, write_json_graph,
  convert_from_json_graph_data, convert_to_json_graph_data

using DataStructures: OrderedDict

using ..WiringDiagramCore
import ..WiringDiagramCore: PortEdgeData

const JSONObject = OrderedDict{String,Any}

# Serialization
###############

""" Serialize a wiring diagram in JSON graph format.
"""
function write_json_graph(diagram::WiringDiagram)::AbstractDict
  write_json_box(diagram, "n")
end

function write_json_box(diagram::WiringDiagram, id::String)
  # Define node, port, and edge IDs.
  in_id, out_id = input_id(diagram), output_id(diagram)
  node_id(v::Int) = v in (in_id, out_id) ? id : "$id:n$v"
  port_id(port::Port) = begin
    is_input = port.box in (in_id, out_id) ?
      port.box == in_id : port.kind == InputPort
    string(is_input ? "in" : "out", ":", port.port)
  end
  #edge_id(i::Int) = "$id:e$i"

  JSONObject(
    "id" => id,
    "ports" => write_json_ports(diagram),
    "children" => [
      write_json_box(box(diagram, v), node_id(v)) for v in box_ids(diagram)
    ],
    "edges" => [
      json_object_with_value(wire.value,
        #"id" => edge_id(i),
        "source" => node_id(wire.source.box),
        "sourcePort" => port_id(wire.source),
        "target" => node_id(wire.target.box),
        "targetPort" => port_id(wire.target),
      )
      for (i, wire) in enumerate(wires(diagram))
    ],
  )
end

function write_json_box(box::Box, id::String)
  json_object_with_value(box.value,
    "id" => id,
    "ports" => write_json_ports(box),
  )
end

function write_json_ports(box::AbstractBox)::AbstractArray
  [
    [
      json_object_with_value(port,
        "id" => "in:$i",
        "portkind" => "input",
      )
      for (i, port) in enumerate(input_ports(box))
    ];
    [
      json_object_with_value(port,
        "id" => "out:$i",
        "portkind" => "output",
      )
      for (i, port) in enumerate(output_ports(box))
    ];
  ]
end

function json_object_with_value(value, args...)
  obj = JSONObject(args...)
  props = convert_to_json_graph_data(value)
  if !isempty(props)
    obj["properties"] = props
  end
  obj
end

convert_to_json_graph_data(value::AbstractDict) = value
convert_to_json_graph_data(value) = Dict("value" => value)
convert_to_json_graph_data(::Nothing) = Dict()

# Deserialization
#################

""" Deserialize a wiring diagram from JSON graph format.
"""
function read_json_graph(
    ::Type{BoxValue}, ::Type{PortValue}, ::Type{WireValue},
    node::AbstractDict)::WiringDiagram where {BoxValue, PortValue, WireValue}
  diagram, ports = read_json_box(BoxValue, PortValue, WireValue, node)
  diagram
end

function read_json_box(
    BoxValue::Type, PortValue::Type, WireValue::Type, node::AbstractDict)
  # Read the ports of the box.
  ports, input_ports, output_ports = read_json_ports(PortValue, node)

  # Dispense with case when node is an atomic box.
  if !haskey(node, "children")
    value = read_json_graph_data(BoxValue, node)
    return (Box(value, input_ports, output_ports), ports)
  end

  # If we get here, we're reading a wiring diagram.
  diagram = WiringDiagram(input_ports, output_ports)
  all_ports = Dict{Tuple{String,String},Port}()
  for (key, port_data) in ports
    all_ports[key] = port_data.kind == InputPort ?
      Port(input_id(diagram), OutputPort, port_data.port) : 
      Port(output_id(diagram), InputPort, port_data.port)
  end

  # Read the nodes recursively.
  for subnode in node["children"]
    box, subports = read_json_box(BoxValue, PortValue, WireValue, subnode)
    v = add_box!(diagram, box)
    for (key, port_data) in subports
      all_ports[key] = Port(v, port_data.kind, port_data.port)
    end
  end

  # Read the edges.
  for edge in node["edges"]
    value = read_json_graph_data(WireValue, edge)
    source = all_ports[(edge["source"], edge["sourcePort"])]
    target = all_ports[(edge["target"], edge["targetPort"])]
    add_wire!(diagram, Wire(value, source, target))
  end

  (diagram, ports)
end

function read_json_ports(PortValue::Type, node::AbstractDict)
  ports = Dict{Tuple{String,String},PortEdgeData}()
  input_ports, output_ports = PortValue[], PortValue[]
  node_id = node["id"]
  for port in node["ports"]
    port_id = port["id"]
    port_kind = port["portkind"]
    value = read_json_graph_data(PortValue, port)
    if port_kind == "input"
      push!(input_ports, value)
      ports[(node_id, port_id)] = PortEdgeData(InputPort, length(input_ports))
    elseif port_kind == "output"
      push!(output_ports, value)
      ports[(node_id, port_id)] = PortEdgeData(OutputPort, length(output_ports))
    else
      error("Invalid port kind: $portkind")
    end
  end
  (ports, input_ports, output_ports)
end

function read_json_graph_data(Value::Type, obj::AbstractDict)
  convert_from_json_graph_data(Value, get(obj, "properties", Dict()))
end

convert_from_json_graph_data(::Type{T}, data::AbstractDict) where T <: AbstractDict = convert(T, data)
convert_from_json_graph_data(::Type{Nothing}, data::AbstractDict) = nothing

function convert_from_json_graph_data(Value::Type, data::AbstractDict)
  @assert length(data) == 1
  first(values(data))::Value
end
function convert_from_json_graph_data(::Type{Symbol}, data::AbstractDict)
  @assert length(data) == 1
  Symbol(first(values(data)))
end

end
