using GWCS
using GWCS: ArrayLabelMapper, ConstantTransform, DictionaryLabelMapper, RangeLabelMapper, RegionSelector, Scale, Shift, map_label
using StaticArrays
using YAML

include("python.jl")

function main()
    python = validation_python()
    reference_script = joinpath(@__DIR__, "python_selector_reference.py")
    reference = YAML.load(read(`$python $reference_script`, String))

    array_mapper = ArrayLabelMapper([0 2; 1 0]; no_label=0)
    regions = RegionSelector(array_mapper, (1 => (Shift(1) & Scale(1)), 2 => (Shift(2) & Scale(2))); undefined_value=-100)
    dictionary_mapper = DictionaryLabelMapper((-2.0 => ConstantTransform{2}(1.0), 3.0 => ConstantTransform{2}(2.0)), 2; inputs_mapping=(1,), atol=0.001)
    range_mapper = RangeLabelMapper(((0.0, 2.0) => ConstantTransform{2}(1.0), (2.0, 4.0) => ConstantTransform{2}(2.0)), 2; inputs_mapping=(1,))

    @assert [map_label(array_mapper, 0.0, 1.0), map_label(array_mapper, 1.0, 0.0), something(map_label(array_mapper, 0.0, 0.0), 0)] == reference["array_labels"]
    @assert [collect(regions(SVector(0.0, 1.0))), collect(regions(SVector(1.0, 0.0)))] == reference["region_outputs"][1:2]
    # Python returns undefined_value where no region applies; GWCS.jl deliberately returns NaN and stores the value only for ASDF round trips.
    @assert all(isnan, regions(SVector(0.0, 0.0))) && reference["region_outputs"][3] == fill(regions.undefined_value, 2)
    @assert [map_label(dictionary_mapper, -2.0005, 8.0), map_label(dictionary_mapper, 3.0, 8.0), something(map_label(dictionary_mapper, 0.0, 8.0), 0.0)] == reference["dictionary_labels"]
    @assert [map_label(range_mapper, 1.0, 8.0), map_label(range_mapper, 3.0, 8.0), something(map_label(range_mapper, 2.0, 8.0), 0.0)] == reference["range_labels"]

    println("selector scalar cases compared: 12")
    println("GWCS.jl uses nothing internally for no label and does not reserve numeric or string label values")
end

main()
