using GWCS: AbstractTransform, AffineTransform, AnnotatedTransform, BoundingBox, CartesianToSpherical, CelestialFrame, CompositeFrame,
            CoordinateFrame, Frame2D, IdentityTransform, Rotation3D, Scale, Shift, SpectralFrame, SphericalToCartesian, Step, StokesFrame,
            TemporalFrame, append_frame, append_transform, axis_names, axis_units, compose, fix_inputs, frame_dimension, get_frame, get_transform,
            has_inverse, insert_frame, insert_transform, parallel, prepend_frame, prepend_transform, remove_frame, rename_frame, set_transform,
            transform

struct PipelineCompileFailure <: AbstractTransform{2,2} end

(::PipelineCompileFailure)(coords::StaticVector{2}) = coords
GWCS.compile_transform(::PipelineCompileFailure) = error("pipeline compilation failed")

@testset "WCS pipeline" begin
    detector = Frame2D(name="detector", axes_names=("x", "y"), unit=("pixel", "pixel"))
    focal = Frame2D(name="focal", axes_names=("u", "v"), unit=("deg", "deg"))
    world = CelestialFrame(name="world")

    detector_to_focal = parallel(Shift(-10.0), Shift(-20.0))
    focal_to_world = compose(parallel(Scale(0.1), Scale(0.1)), SphericalToCartesian(), Rotation3D([10.0], "z"), CartesianToSpherical())
    wcs = WCS(Step(detector, detector_to_focal), Step(focal, focal_to_world), Step(world); name="test", pixel_shape=(100, 200), bounding_box=BoundingBox(((-0.5, 99.5), (-0.5, 199.5))))

    @test available_frames(wcs) == ("detector", "focal", "world")
    @test input_frame(wcs) === detector
    @test output_frame(wcs) === world
    @test wcs.name == "test"
    @test wcs.pixel_shape == (100, 200)
    mutable_shape = Union{Int,Nothing}[100, nothing]
    normalized_shape = with_pixel_shape(wcs, mutable_shape)
    mutable_shape[1] = 50
    @test normalized_shape.pixel_shape === (100, nothing)
    @test SVector(10.0, 20.0) in wcs.bounding_box
    @test !(SVector(101.0, 20.0) in wcs.bounding_box)
    c_order_box = BoundingBox(((10.0, 20.0), (1.0, 2.0)); order=:C)
    @test c_order_box.order === :F
    @test SVector(1.5, 15.0) in c_order_box
    @test !(SVector(15.0, 1.5) in c_order_box)

    focal_position = transform(wcs, "detector", "focal", 11.0, 22.0)
    @test focal_position == SVector(1.0, 2.0)
    world_position = transform(wcs, "detector", "world", 11.0, 22.0)
    @test transform(wcs, "world", "detector", world_position) ≈ SVector(11.0, 22.0) atol=1e-10
    @test get_transform(wcs, "focal", "focal") isa IdentityTransform{2}

    compact = sprint(show, wcs)
    @test compact == "WCS(\"test\", 2 → 2, 3 steps)"
    @test !occursin("Tuple{", sprint(show, wcs.steps))
    detailed = sprint(show, MIME("text/plain"), wcs)
    @test occursin("WCS \"test\" (2 → 2), 3 steps", detailed)
    @test occursin("pixel_shape: (100, 200)", detailed)
    @test occursin("bounding_box:", detailed)
    @test occursin("1. frame: CoordinateFrame", detailed)
    @test occursin("AffineTransform", detailed)
    @test occursin("offset=[-10.0, -20.0]", detailed)
    annotated = Step(detector, AnnotatedTransform(AffineTransform([1.0 2.0; 3.0 4.0], [5.0, 6.0]); name="velocity_aberration", inputs=("x", "y"), outputs=("u", "v")))
    annotated_text = sprint(show, MIME("text/plain"), annotated)
    @test occursin("name=\"velocity_aberration\"", annotated_text)
    @test occursin("matrix=[1.0 2.0; 3.0 4.0]", annotated_text)
    @test occursin("offset=[5.0, 6.0]", annotated_text)

    @test_throws "frame names must be unique" WCS(Step(detector, detector_to_focal), Step(detector))
    @test_throws DimensionMismatch WCS(Step(detector, Shift(1.0)), Step(world))

    sliced = fix_inputs(wcs, "y" => 20.0)
    @test pixel_n_dim(sliced) == 1
    @test world_n_dim(sliced) == 2
    @test input_frame(sliced).axes_names == ("x",)
    @test input_frame(sliced).unit == ("pixel",)
    @test sliced.pixel_shape == (100,)
    @test sliced.bounding_box.intervals == ((-0.5, 99.5),)
    @test pixel_to_world(sliced, 11.0) ≈ pixel_to_world(wcs, 11.0, 20.0)
    @test world_to_pixel(sliced, pixel_to_world(sliced, 11.0)) ≈ SVector(11.0) atol=1e-10
    @test all(isnan, world_to_pixel(sliced, pixel_to_world(wcs, 11.0, 21.0)))
    @test_throws DomainError world_to_pixel(sliced, pixel_to_world(wcs, 11.0, 21.0); strict=true)
    @test_throws "lies outside the WCS bounding box" fix_inputs(wcs, "y" => 300.0)

    constant_wcs = fix_inputs(wcs, ("x" => 11.0, "y" => 20.0))
    @test pixel_n_dim(constant_wcs) == 0
    @test input_frame(constant_wcs).axes_names == ()
    @test constant_wcs.pixel_shape == ()
    @test constant_wcs.bounding_box === nothing
    @test pixel_to_world(constant_wcs) ≈ pixel_to_world(wcs, 11.0, 20.0)
    @test pixel_to_world(constant_wcs, ()) ≈ pixel_to_world(wcs, 11.0, 20.0)
    @test world_to_pixel(constant_wcs, pixel_to_world(constant_wcs)) == SVector{0,Float64}()
    @test eltype(world_to_pixel(constant_wcs, pixel_to_world(constant_wcs))) === Float64
end

@testset "Pipeline editing" begin
    input = Frame2D(name="input", axes_names=("x", "y"))
    middle = Frame2D(name="middle", axes_names=("u", "v"))
    output = Frame2D(name="output", axes_names=("a", "b"))
    first_transform = Shift(1.0) & Shift(2.0)
    second_transform = Scale(2.0) & Scale(3.0)
    box = BoundingBox(((0.0, 9.0), (0.0, 19.0)))
    wcs = WCS(Step(input, first_transform), Step(middle, second_transform), Step(output); name="editable", pixel_shape=(10, 20), bounding_box=box)
    pixel = SVector(1.0, 2.0)
    world = SVector(4.0, 12.0)

    @test get_frame(wcs, "middle") === middle
    @test get_frame(wcs, :middle) === middle
    @test get_frame(wcs, Frame2D(name="middle")) === middle
    @test_throws "frame missing is not in the pipeline" get_frame(wcs, :missing)

    replaced = set_transform(wcs, "input", "middle", IdentityTransform{2}())
    @test pixel_to_world(replaced, pixel) == SVector(2.0, 6.0)
    @test pixel_to_world(wcs, pixel) == world
    @test_throws "not adjacent in forward order" set_transform(wcs, "input", "output", IdentityTransform{2}())
    @test_throws DimensionMismatch set_transform(wcs, "input", "middle", Shift(1.0))
    singular = set_transform(wcs, "middle", "output", AffineTransform(@SMatrix [1.0 0.0; 0.0 0.0]))
    @test !has_inverse(singular.forward)

    before_middle = insert_transform(wcs, "middle", Scale(10.0) & Scale(10.0))
    after_middle = insert_transform(wcs, "middle", Scale(10.0) & Scale(10.0); after=true)
    @test pixel_to_world(before_middle, pixel) == SVector(40.0, 120.0)
    @test pixel_to_world(after_middle, pixel) == SVector(40.0, 120.0)
    @test transform(before_middle, "input", "middle", pixel) == SVector(20.0, 40.0)
    @test transform(after_middle, "input", "middle", pixel) == SVector(2.0, 4.0)
    @test transform(after_middle, "middle", "output", 2.0, 4.0) == SVector(40.0, 120.0)
    @test_throws "before the first frame" insert_transform(wcs, "input", IdentityTransform{2}())
    @test_throws "after the last frame" insert_transform(wcs, "output", IdentityTransform{2}(); after=true)
    @test pixel_to_world(prepend_transform(wcs, Scale(10.0) & Scale(10.0)), pixel) == SVector(22.0, 66.0)
    @test pixel_to_world(append_transform(wcs, Shift(1.0) & Shift(1.0)), pixel) == SVector(5.0, 13.0)

    inserted_before = insert_frame(wcs, Frame2D(name="inserted_before"), Scale(10.0) & Scale(10.0), middle)
    @test available_frames(inserted_before) == ("input", "inserted_before", "middle", "output")
    @test transform(inserted_before, "input", "inserted_before", pixel) == SVector(2.0, 4.0)
    @test transform(inserted_before, "inserted_before", "middle", 2.0, 4.0) == SVector(20.0, 40.0)
    @test pixel_to_world(inserted_before, pixel) == SVector(40.0, 120.0)

    inserted_after = insert_frame(wcs, middle, Scale(10.0) & Scale(10.0), Frame2D(name="inserted_after"))
    @test available_frames(inserted_after) == ("input", "middle", "inserted_after", "output")
    @test transform(inserted_after, "middle", "inserted_after", 2.0, 4.0) == SVector(20.0, 40.0)
    @test transform(inserted_after, "inserted_after", "output", 20.0, 40.0) == SVector(40.0, 120.0)
    @test pixel_to_world(inserted_after, pixel) == SVector(40.0, 120.0)
    @test_throws "exactly one inserted-frame endpoint" insert_frame(wcs, input, IdentityTransform{2}(), output)
    @test_throws "exactly one inserted-frame endpoint" insert_frame(wcs, Frame2D(name="new1"), IdentityTransform{2}(), Frame2D(name="new2"))
    @test_throws "new input endpoint must be a coordinate frame" insert_frame(wcs, "new", IdentityTransform{2}(), middle)

    removable = insert_frame(wcs, middle, IdentityTransform{2}(), Frame2D(name="temporary"))
    removed = remove_frame(removable, "temporary")
    @test available_frames(removed) == available_frames(wcs)
    @test pixel_to_world(removed, pixel) == world
    @test world_to_pixel(removed, world) ≈ pixel
    @test_throws "only an intermediate frame" remove_frame(wcs, "input")
    @test_throws "only an intermediate frame" remove_frame(wcs, "output")

    renamed = rename_frame(wcs, "middle", "focal")
    @test available_frames(renamed) == ("input", "focal", "output")
    @test pixel_to_world(renamed, pixel) == world
    @test_throws "already in the pipeline" rename_frame(wcs, "middle", "input")
    @test_throws "must not be empty" rename_frame(wcs, "middle", "")

    frame_families = (
        CoordinateFrame{1}(name="coordinate"),
        CelestialFrame(name="celestial"),
        SpectralFrame(name="spectral"),
        TemporalFrame(name="temporal"),
        StokesFrame(name="stokes"),
        CompositeFrame((CelestialFrame(name="sky"), SpectralFrame(name="wave", axes_order=3)); name="composite")
    )
    for frame in frame_families
        dimension = frame_dimension(frame)
        family_wcs = WCS(Step(frame, IdentityTransform{dimension}()), Step(CoordinateFrame{dimension}(name="terminal")))
        family_renamed = rename_frame(family_wcs, frame, "renamed")
        @test typeof(input_frame(family_renamed)) == typeof(frame)
        @test input_frame(family_renamed).name == "renamed"
        @test axis_names(input_frame(family_renamed)) == axis_names(frame)
        @test axis_units(input_frame(family_renamed)) == axis_units(frame)
    end

    @test_throws "explicit pixel_shape replacement" prepend_frame(wcs, Frame2D(name="new_input"), IdentityTransform{2}())
    @test_throws "explicit bounding_box replacement" prepend_frame(wcs, Frame2D(name="new_input"), IdentityTransform{2}(); pixel_shape=(30, 40))
    prepended = prepend_frame(wcs, Frame2D(name="new_input"), IdentityTransform{2}(); pixel_shape=(30, 40), bounding_box=BoundingBox(((0.0, 29.0), (0.0, 39.0))))
    @test available_frames(prepended) == ("new_input", "input", "middle", "output")
    @test prepended.pixel_shape == (30, 40)
    @test prepended.bounding_box.intervals == ((0.0, 29.0), (0.0, 39.0))
    @test pixel_to_world(prepended, pixel) == world
    cleared = prepend_frame(wcs, Frame2D(name="unbounded_input"), IdentityTransform{2}(); pixel_shape=nothing, bounding_box=nothing)
    @test cleared.pixel_shape === nothing
    @test cleared.bounding_box === nothing

    appended = append_frame(wcs, IdentityTransform{2}(), Frame2D(name="new_output"))
    @test available_frames(appended) == ("input", "middle", "output", "new_output")
    @test appended.pixel_shape == wcs.pixel_shape
    @test appended.bounding_box === box
    @test pixel_to_world(appended, pixel) == world

    # Edits chain by rebinding, and a failed edit leaves the WCS it was given unchanged.
    edited = set_transform(wcs, "input", "middle", IdentityTransform{2}())
    @test pixel_to_world(edited, pixel) == SVector(2.0, 6.0)
    edited = rename_frame(edited, "middle", "focal")
    edited = insert_frame(edited, "focal", IdentityTransform{2}(), Frame2D(name="temporary"))
    edited = remove_frame(edited, "temporary")
    @test available_frames(edited) == ("input", "focal", "output")
    @test_throws "pipeline compilation failed" set_transform(edited, "input", "focal", PipelineCompileFailure())

    allocated_forward(edited, coords) = @allocated pixel_to_world(edited, coords)
    allocated_backward(edited, coords) = @allocated world_to_pixel(edited, coords)
    pixel_to_world(renamed, pixel)
    world_to_pixel(renamed, world)
    @test @inferred(pixel_to_world(renamed, pixel)) == world
    @test @inferred(world_to_pixel(renamed, world)) ≈ pixel
    @test allocated_forward(renamed, pixel) == 0
    @test allocated_backward(renamed, world) == 0
end

@testset "Frame axis ordering" begin
    input = CoordinateFrame{3}(name="input", axes_names=("a", "b", "c"), axes_order=(2, 3, 1))
    middle = CoordinateFrame{3}(name="middle", axes_names=("p", "q", "r"), axes_order=(3, 1, 2))
    output = CoordinateFrame{3}(name="output", axes_names=("u", "v", "w"), axes_order=(2, 1, 3))
    first_transform = AffineTransform(@SMatrix [10.0 0.0 0.0; 0.0 100.0 0.0; 0.0 0.0 1000.0])
    wcs = WCS(Step(input, first_transform), Step(middle, IdentityTransform{3}()), Step(output);
        pixel_shape=(10, 20, 30), bounding_box=BoundingBox(((0.0, 9.0), (0.0, 19.0), (0.0, 29.0))))

    # As in Python GWCS, coordinates are in transform order and axes_order only labels them: frame axis i is at
    # transform position axes_order[i], so the accessors report the input's names as (c, a, b).
    @test axis_names(input) == ("c", "a", "b")
    @test pixel_axis_names(wcs) == ("c", "a", "b")
    @test world_axis_names(wcs) == ("v", "u", "w")
    pixels = SVector(1.0, 2.0, 3.0)
    world = SVector(10.0, 200.0, 3000.0)
    @test pixel_to_world(wcs, pixels) == world
    @test world_to_pixel(wcs, world) == pixels
    @test in_domain(wcs, pixels)
    @test in_image(wcs, world)
    @test !in_domain(wcs, SVector(1.0, 20.0, 3.0))
    @test transform(wcs, "input", "middle", pixels...) == world
    @test transform(wcs, "middle", "output", world...) == world
    @test transform(wcs, "output", "input", world...) == pixels
    @test axis_correlation_matrix(wcs) == [i == j for i in 1:3, j in 1:3]

    batch = [1.0 4.0; 2.0 5.0; 3.0 6.0]
    world_batch = pixel_to_world(wcs, batch)
    @test world_batch == [10.0 40.0; 200.0 500.0; 3000.0 6000.0]
    @test world_to_pixel(wcs, world_batch) == batch

    # "a" is transform input 2, so fixing it removes the second pixel_shape entry and bounding-box interval.
    sliced = fix_inputs(wcs, "a" => 1.0)
    @test input_frame(sliced).axes_names == ("b", "c")
    @test input_frame(sliced).axes_order == (2, 1)
    @test pixel_axis_names(sliced) == ("c", "b")
    @test sliced.pixel_shape == (10, 30)
    @test sliced.bounding_box.intervals == ((0.0, 9.0), (0.0, 29.0))
    @test pixel_to_world(sliced, 2.0, 3.0) == SVector(20.0, 100.0, 3000.0)
    @test world_to_pixel(sliced, SVector(20.0, 100.0, 3000.0)) == SVector(2.0, 3.0)
    @test fix_inputs(wcs, 2 => 1.0).steps == sliced.steps
    @test_throws "lies outside the WCS bounding box" fix_inputs(wcs, "a" => 20.0)

    allocated_forward(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    allocated_backward(wcs, coords) = @allocated world_to_pixel(wcs, coords)
    pixel_to_world(wcs, pixels)
    world_to_pixel(wcs, world)
    @test @inferred(pixel_to_world(wcs, pixels)) == world
    @test @inferred(world_to_pixel(wcs, world)) == pixels
    @test allocated_forward(wcs, pixels) == 0
    @test allocated_backward(wcs, world) == 0

    @test_throws DimensionMismatch CoordinateFrame{2}(name="bad", axes_order=(1,))
    @test_throws "entries must be integers" CoordinateFrame{2}(name="bad", axes_order=(true, 2))
    @test_throws "entries must be positive" CoordinateFrame{2}(name="bad", axes_order=(0, 1))
    @test_throws "entries must be unique" CoordinateFrame{2}(name="bad", axes_order=(1, 1))
    incomplete = CoordinateFrame{2}(name="incomplete", axes_order=(1, 3))
    @test_throws "must be a permutation of 1:2" WCS(Step(incomplete, IdentityTransform{2}()), Step(CoordinateFrame{2}(name="valid")))
end
