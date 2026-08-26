module GWCSASDFExt

using ASDF
using GWCS

GWCS._tag(node::ASDF.TaggedMapping) = node.tag
GWCS._tag(node::ASDF.TaggedSequence) = node.tag
GWCS._tag(node::ASDF.TaggedScalar) = node.tag
GWCS._unwrap(node::ASDF.TaggedScalar) = String(node)

function GWCS._materialize_array(array::ASDF.NDArray)
    value = array[]
    return value isa AbstractMatrix ? permutedims(value) : value
end

end
