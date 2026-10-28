_display_type(value) = nameof(typeof(value))

function _show_properties(io::IO, value)
    print(io, _display_type(value), "(")
    for (index, field) in enumerate(fieldnames(typeof(value)))
        index > 1 && print(io, ", ")
        print(io, field, "=")
        show(io, getfield(value, field))
    end
    print(io, ")")
end

Base.show(io::IO, transform::AbstractTransform) = _show_properties(io, transform)

function Base.show(io::IO, transform::AnnotatedTransform)
    print(io, "AnnotatedTransform(")
    transform.name === nothing || (print(io, "name="); show(io, transform.name); print(io, ", "))
    show(io, transform.transform)
    transform.inputs === nothing || (print(io, ", inputs="); show(io, transform.inputs))
    transform.outputs === nothing || (print(io, ", outputs="); show(io, transform.outputs))
    print(io, ")")
end

function Base.show(io::IO, transform::ComposedTransform)
    print(io, "ComposedTransform(")
    join(io, transform.transforms, " | ")
    print(io, ")")
end

function Base.show(io::IO, transform::ParallelTransform)
    print(io, "ParallelTransform(")
    join(io, transform.transforms, " & ")
    print(io, ")")
end

function Base.show(io::IO, transform::WithInverse)
    print(io, "WithInverse(")
    show(io, transform.forward)
    print(io, ", inverse=")
    show(io, transform.backward)
    print(io, ")")
end

function Base.show(io::IO, frame::AbstractCoordinateFrame)
    print(io, _display_type(frame), "(")
    show(io, frame.name)
    print(io, ", axes=")
    show(io, axis_names(frame))
    print(io, ")")
end

function Base.show(io::IO, step::Step)
    print(io, "Step(")
    show(io, step.frame)
    if step.transform !== nothing
        print(io, ", ")
        show(io, step.transform)
    end
    print(io, ")")
end

function Base.show(io::IO, wcs::WCS{Nin,Nout}) where {Nin,Nout}
    print(io, "WCS(")
    show(io, wcs.name)
    print(io, ", ", Nin, " → ", Nout, ", ", length(wcs.steps), " steps)")
end

function _show_transform_tree(io::IO, transform::AbstractTransform, indent)
    print(io, " "^indent)
    show(io, transform)
end

function _show_transform_tree(io::IO, transform::AnnotatedTransform, indent)
    print(io, " "^indent, "AnnotatedTransform")
    transform.name === nothing || (print(io, " name="); show(io, transform.name))
    transform.inputs === nothing || (print(io, " inputs="); show(io, transform.inputs))
    transform.outputs === nothing || (print(io, " outputs="); show(io, transform.outputs))
    print(io, '\n')
    _show_transform_tree(io, transform.transform, indent + 2)
end

function _show_transform_tree(io::IO, transform::Union{ComposedTransform,ParallelTransform}, indent)
    operator = transform isa ComposedTransform ? "ComposedTransform (|)" : "ParallelTransform (&)"
    print(io, " "^indent, operator)
    for child in transform.transforms
        print(io, '\n')
        _show_transform_tree(io, child, indent + 2)
    end
end

function _show_transform_tree(io::IO, transform::WithInverse, indent)
    print(io, " "^indent, "WithInverse\n", " "^(indent + 2), "forward:\n")
    _show_transform_tree(io, transform.forward, indent + 4)
    print(io, '\n', " "^(indent + 2), "inverse:\n")
    _show_transform_tree(io, transform.backward, indent + 4)
end

function Base.show(io::IO, ::MIME"text/plain", step::Step)
    print(io, "Step\n  frame: ")
    _show_properties(io, step.frame)
    if step.transform !== nothing
        print(io, "\n  transform (", input_dimension(step.transform), " → ", output_dimension(step.transform), "):\n")
        _show_transform_tree(io, step.transform, 4)
    end
end

function Base.show(io::IO, ::MIME"text/plain", wcs::WCS{Nin,Nout}) where {Nin,Nout}
    print(io, "WCS ")
    show(io, wcs.name)
    print(io, " (", Nin, " → ", Nout, "), ", length(wcs.steps), " steps")
    wcs.pixel_shape === nothing || (print(io, "\n  pixel_shape: "); show(io, wcs.pixel_shape))
    wcs.bounding_box === nothing || (print(io, "\n  bounding_box: "); show(io, wcs.bounding_box))
    for (index, step) in enumerate(wcs.steps)
        print(io, "\n\n", index, ". frame: ")
        _show_properties(io, step.frame)
        if step.transform !== nothing
            print(io, "\n   transform (", input_dimension(step.transform), " → ", output_dimension(step.transform), "):\n")
            _show_transform_tree(io, step.transform, 5)
        end
    end
end
