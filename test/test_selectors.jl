using GWCS: AffineTransform, ArrayLabelMapper, BoundingBox, CallableLabelMapper, CallableTransform, ConstantTransform, DictionaryLabelMapper,
            Frame2D, NoInverse, RangeLabelMapper, RegionSelector, Scale, Shift, Step, TransformLabelMapper, WithInverseLabelMapper,
            annotate_transform, compile_label_mapper, compile_transform, input_dimension, input_names, inverse_label_mapper, inverse_transform,
            label_dependencies, label_mapper_name, label_type, map_label, output_dimension, output_names, selected_transform, transform,
            transform_dependencies, transform_jacobian, transform_name

allocated_selector(transform, coords) = @allocated transform(coords)

@testset "Label mappers" begin
    labels = [0 2; 1 0]
    mapper = ArrayLabelMapper(labels; no_label=0, inputs=("x", "y"), name="detector labels")
    @test input_dimension(mapper) == 2
    @test label_type(mapper) == Int
    @test input_names(mapper) == ("x", "y")
    @test label_mapper_name(mapper) == "detector labels"
    @test map_label(mapper, 0, 1) == 1
    @test map_label(mapper, 1, 0) == 2
    @test map_label(mapper, -0.5, 0) === nothing
    @test map_label(mapper, -0.5001, 0) === nothing
    @test map_label(mapper, 2, 0) === nothing
    @test map_label(mapper, NaN, 0) === nothing
    @test map_label(mapper, floatmax(Float64), 0) === nothing
    @test label_dependencies(mapper) == SVector(true, true)

    zero_is_label = ArrayLabelMapper(labels)
    @test map_label(zero_is_label, 0, 0) == 0

    string_mapper = ArrayLabelMapper(["" "right"; "left" ""]; no_label="")
    @test map_label(string_mapper, 0, 1) == "left"
    @test map_label(string_mapper, 1, 0) == "right"
    @test map_label(string_mapper, 0, 0) === nothing

    projected = ArrayLabelMapper{3}(labels; inputs_mapping=(1, 3))
    @test map_label(projected, 1, 99, 0) == 2
    @test label_dependencies(projected) == SVector(true, false, true)
    @test_throws DimensionMismatch ArrayLabelMapper{3}(labels; inputs_mapping=(1,))
    @test_throws ArgumentError ArrayLabelMapper{3}(labels; inputs_mapping=(1, 1))

    callable = CallableLabelMapper{2,Symbol}((x, y) -> x < y ? :lower : :upper; dependencies=(true, true), inputs=("x", "y"))
    @test map_label(callable, 1, 2) == :lower
    @test map_label(callable, 3, 2) == :upper
    @test label_dependencies(callable) == SVector(true, true)
    @test_throws DimensionMismatch CallableLabelMapper{2,Symbol}((x, y) -> :label; dependencies=(true,))

    transform_mapper = TransformLabelMapper{Int}(ConstantTransform{2}(3.0))
    @test map_label(transform_mapper, 1, 2) == 3
    @test label_dependencies(transform_mapper) == SVector(false, false)
    @test compile_label_mapper(transform_mapper) isa TransformLabelMapper

    dictionary = DictionaryLabelMapper((-2.0 => CallableTransform{2,1}((x, y) -> x + y), 3.0 => ConstantTransform{2}(7.0)), 2)
    @test map_label(dictionary, -2, 1) == -1
    @test map_label(dictionary, 3, 1) == 7
    @test map_label(dictionary, 0, 1) === nothing
    @test map_label(DictionaryLabelMapper((1.0 => ConstantTransform{2}(2.0),), 2; inputs_mapping=(2,), atol=0.1, rtol=0), 8, 1.05) == 2
    numpy_tolerance = DictionaryLabelMapper((1000.0 => ConstantTransform{1}(1.0),), 1; atol=0.1)
    @test map_label(numpy_tolerance, 1000.109) == 1
    @test map_label(numpy_tolerance, 1000.111) === nothing
    @test label_dependencies(dictionary) == SVector(true, true)
    @test compile_label_mapper(dictionary) isa DictionaryLabelMapper
    ambiguous = DictionaryLabelMapper((1.0 => ConstantTransform{1}(1.0), 1.05 => ConstantTransform{1}(2.0)), 1; atol=0.1, rtol=0)
    @test map_label(ambiguous, 1.02) == 2
    @test_throws DimensionMismatch DictionaryLabelMapper((1.0 => Shift(1),), 2)

    ranges = RangeLabelMapper(((0.0, 2.0) => ConstantTransform{2}(1.0), (2.0, 4.0) => CallableTransform{2,1}((x, y) -> x + y)), 2)
    @test map_label(ranges, 1, 5) == 1
    @test map_label(ranges, 3, 5) == 8
    @test map_label(ranges, 0, 5) === nothing
    @test map_label(ranges, 2, 5) === nothing
    @test map_label(ranges, 4, 5) === nothing
    @test map_label(ranges, Inf, 5) === nothing
    @test label_dependencies(ranges) == SVector(true, true)
    @test compile_label_mapper(ranges) isa RangeLabelMapper
    @test_throws ArgumentError RangeLabelMapper(((0.0, 2.0) => ConstantTransform{1}(1.0), (1.0, 3.0) => ConstantTransform{1}(2.0)), 1)

    reverse_callable = CallableLabelMapper{1,Symbol}(x -> x < 0 ? :negative : :positive)
    paired = WithInverseLabelMapper(CallableLabelMapper{2,Symbol}((x, y) -> x < 0 ? :negative : :positive), reverse_callable)
    @test map_label(paired, -1, 2) == :negative
    @test input_dimension(inverse_label_mapper(paired)) == 1
    @test inverse_label_mapper(inverse_label_mapper(paired)) isa WithInverseLabelMapper
    @test inverse_label_mapper(mapper) isa NoInverse

    @test_throws DimensionMismatch map_label(mapper, [1.0])
    @test_throws MethodError map_label(mapper, "x", "y")

    zero_dimensional = CallableLabelMapper{0,Int}(() -> 1)
    @test map_label(zero_dimensional) == 1
end

@testset "Region selectors" begin
    mapper = CallableLabelMapper{2,Symbol}((x, y) -> x < 0 ? :left : x > 10 ? :missing : :right; dependencies=(true, false))
    left = CallableTransform{2,2}((x, y) -> (x + 1, y - 1); jacobian=(x, y) -> [1 0; 0 1], dependencies=[true false; false true])
    right = AffineTransform([2.0 0; 0 3.0], [0.0, 0.0])
    selector = RegionSelector(mapper, (:left => left, :right => right))

    @test input_dimension(selector) == 2
    @test output_dimension(selector) == 2
    @test selector(SVector(-2.0, 4.0)) == SVector(-1.0, 3.0)
    @test selector(SVector(2.0, 4.0)) == SVector(4.0, 12.0)
    @test all(isnan, selector(SVector(11.0, 4.0)))
    @test selected_transform(selector, :left) === left
    @test_throws KeyError selected_transform(selector, :missing)
    @test transform_jacobian(selector, SVector(-2.0, 4.0)) == SMatrix{2,2}(1, 0, 0, 1)
    @test transform_jacobian(selector, SVector(2.0, 4.0)) == SMatrix{2,2}(2, 0, 0, 3)
    @test all(isnan, transform_jacobian(selector, SVector(11.0, 4.0)))
    @test transform_dependencies(selector) == SMatrix{2,2,Bool}(true, true, false, true)
    @test compile_transform(selector) isa RegionSelector
    @test inverse_transform(selector) isa NoInverse

    # A custom undefined value is stored for ASDF round trips; evaluation still marks the miss with NaN.
    filled = RegionSelector(mapper, (:left => left, :right => right); undefined_value=-100)
    @test filled.undefined_value == -100
    @test all(isnan, filled(SVector(11.0, 4.0)))

    forward_mapper = CallableLabelMapper{2,Int}((x, y) -> x < 5 ? 1 : 2)
    backward_mapper = CallableLabelMapper{2,Int}((x, y) -> x < 6 ? 1 : 2)
    paired_mapper = WithInverseLabelMapper(forward_mapper, backward_mapper)
    reversible = RegionSelector(paired_mapper, (1 => (Shift(1) & Shift(1)), 2 => (Shift(2) & Shift(2))))
    backward = inverse_transform(reversible)
    @test backward isa RegionSelector
    # As in Python GWCS, the inverse selector does not inherit a custom undefined value.
    @test isnan(inverse_transform(RegionSelector(paired_mapper, (1 => (Shift(1) & Shift(1)), 2 => (Shift(2) & Shift(2))); undefined_value=-100)).undefined_value)
    @test reversible(SVector(1.0, 2.0)) == SVector(2.0, 3.0)
    @test backward(SVector(2.0, 3.0)) == SVector(1.0, 2.0)

    no_branch_inverse = RegionSelector(paired_mapper, (1 => (Shift(1) & Shift(1)), 2 => CallableTransform{2,2}((x, y) -> (x + 2, y + 2))))
    @test inverse_transform(no_branch_inverse) isa NoInverse
    # Branch inverses are stored only to seed the numerical inverse of a selector whose label mapper has no inverse.
    @test reversible.inverses === nothing
    @test selector.inverses isa Tuple

    @test_throws ArgumentError RegionSelector(mapper, ())
    @test_throws ArgumentError RegionSelector(mapper, (:left => left, :left => right))
    @test_throws ArgumentError RegionSelector(mapper, (1 => left, 2 => right))
    @test_throws DimensionMismatch RegionSelector(mapper, (:left => Shift(1), :right => right))

    annotated = annotate_transform(selector; name="regions", inputs=("x", "y"), outputs=("u", "v"))
    @test transform_name(annotated) == "regions"
    @test input_names(annotated) == ("x", "y")
    @test output_names(annotated) == ("u", "v")

    coords = SVector(-2.0, 4.0)
    selector(coords)
    @test @inferred(selector(coords)) == SVector(-1.0, 3.0)
    @test allocated_selector(selector, coords) == 0

    zero_label = RegionSelector(ArrayLabelMapper(reshape([0], 1, 1)), (0 => (Shift(1) & Shift(2)),))
    @test zero_label(SVector(0.0, 0.0)) == SVector(1.0, 2.0)
    empty_label = RegionSelector(ArrayLabelMapper(reshape([""], 1, 1)), ("" => (Shift(1) & Shift(2)),))
    @test empty_label(SVector(0.0, 0.0)) == SVector(1.0, 2.0)
    zero_input = RegionSelector(CallableLabelMapper{0,Int}(() -> 1), (1 => ConstantTransform{0}(4.0),))
    @test zero_input(SVector{0,Float64}()) == SVector(4.0)
end

@testset "Region selector WCS" begin
    detector = Frame2D(; name="detector")
    world = Frame2D(; name="world")
    mapper = ArrayLabelMapper([0 2; 1 0]; no_label=0)
    selector = RegionSelector(mapper, (1 => (Shift(1) & Scale(1)), 2 => (Shift(2) & Scale(2))))
    box = BoundingBox(((0.0, 1.0), (0.0, 1.0)))
    wcs = WCS(Step(detector, selector), Step(world); bounding_box=box)

    @test pixel_to_world(wcs, 0, 1) == SVector(1.0, 1.0)
    @test pixel_to_world(wcs, 1, 0) == SVector(3.0, 0.0)
    @test all(isnan, pixel_to_world(wcs, 0, 0))
    @test all(isnan, pixel_to_world(wcs, 2, 0))
    # A pixel in no region is undefined; strict mode reports it, and the bounding-box miss separately.
    @test_throws "RegionSelector in the step from frame \"detector\" to \"world\" is undefined" pixel_to_world(wcs, 0, 0; strict=true)
    @test_throws "outside the WCS bounding box" pixel_to_world(wcs, 2, 0; strict=true)
    @test in_domain(wcs, 0, 0)
    @test !in_domain(wcs, 2, 0)

    pixels = [0.0 1 0; 1 0 0]
    result = pixel_to_world(wcs, pixels)
    @test result[:, 1] == [1.0, 1.0]
    @test result[:, 2] == [3.0, 0.0]
    @test all(isnan, result[:, 3])

    # Without an inverse label mapper the WCS inverts numerically. The default start (1, 1) is in no region, so
    # Newton starts from the branch inverse that lands back in its own region instead.
    unboxed = WCS(Step(detector, selector), Step(world))
    @test unboxed.backward isa NoInverse
    for pixel in (SVector(0.0, 1.0), SVector(1.0, 0.0), SVector(0.2, 0.9), SVector(0.9, -0.3))
        @test world_to_pixel(unboxed, pixel_to_world(unboxed, pixel)) ≈ pixel atol=1e-12
        @test in_image(unboxed, pixel_to_world(unboxed, pixel))
    end
    target = pixel_to_world(unboxed, 0.2, 0.9)
    @test @inferred(world_to_pixel(unboxed, target)) ≈ SVector(0.2, 0.9) atol=1e-12
    allocated_inverse(wcs, value) = @allocated world_to_pixel(wcs, value)
    allocated_inverse(unboxed, target)
    @test allocated_inverse(unboxed, target) == 0
end
