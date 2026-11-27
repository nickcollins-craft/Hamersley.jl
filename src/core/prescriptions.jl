#=
Value prescriptions.

A boundary condition says *what* is prescribed (a trace, or the quantity conjugate to it),
*where* (a region and a selection of components), and *how much*: the value, as a
function of the reference position ``X`` and the time ``t``. This file defines the last of
these.

A prescribed value only ever enters the right-hand side: the constrained values of an
essential condition, or the load of a natural one. It never enters the tangent, so no
prescription can change whether a problem is linear. The one exception is a follower load
whose value depends on the current position ``x = X + u`` (a hydrostatic pressure on a
deforming body), and that dependence is declared, so the linearity and coupling analyses
see it.

Any prescription could be written as a function of ``(X, t)``; the other kinds exist
because the system can see inside them. A tabulated history declares its breakpoints
(where the step controller should land, since a kink smeared across a step loses accuracy
for rate-dependent materials) and that it is uniform in space (so that its shape can be
checked before anything is discretised).
=#


"""
    AbstractPrescription

Supertype of the ways a boundary value can be prescribed.
"""
abstract type AbstractPrescription end


"""
    Constant(value)

The value `value`, everywhere on the region and at all times. Plain numbers, tuples and
vectors given as a value are constants.

# Fields
- `value::V`: the value: a number for a scalar, a tuple or vector for several components.
"""
struct Constant{V} <: AbstractPrescription
    value::V
end


"""
    Prescribed(f; breakpoints = (), current_position = false)

The value ``f(X, t)``, a function of the reference position and the time. A function given
as a value is a `Prescribed` value without breakpoints.

With `current_position = true` the function is called as ``f(X, x, t)``, with ``x = X + u``
the current position. This is only allowed for follower loads (see
[`is_follower_quantity`](@ref)), whose dependence on the deformation is linearised with
them.

# Fields
- `f::F`: the function.
- `breakpoints::Tuple{Vararg{Float64}}`: the times at which ``f`` has a kink, where the step
  controller should land.
- `current_position::Bool`: whether ``f`` also depends on the current position.
"""
struct Prescribed{F} <: AbstractPrescription
    f::F
    breakpoints::Tuple{Vararg{Float64}}
    current_position::Bool
end

function Prescribed(f; breakpoints=(), current_position::Bool=false)
    # Breakpoints are kept sorted, so that the step controller can walk through them.
    sorted_breakpoints = Tuple(sort!(Float64[time for time in breakpoints]))
    prescription = Prescribed(f, sorted_breakpoints, current_position)
    return prescription
end


"""
    Tabulated(times, values)

A history that is piecewise linear in time and uniform in space: `values[i]` at
`times[i]`, interpolated linearly in between, and held at the first and last values
outside the table. A ramp is a table of two rows.

# Fields
- `times::Vector{Float64}`: the times of the rows, strictly increasing.
- `values::Vector{V}`: the value at each time.
"""
struct Tabulated{V} <: AbstractPrescription
    times::Vector{Float64}
    values::Vector{V}
end

function Tabulated(times, values)
    # We store the rows as vectors, whatever collections they were given in; the
    # verification checks that they are consistent.
    table = Tabulated(Float64[time for time in times], collect(values))
    return table
end


"""
    prescription(value)

Return the prescription that `value` names: a prescription is returned as it is, a
function is [`Prescribed`](@ref), and anything else is [`Constant`](@ref).
"""
prescription(value::AbstractPrescription) = value
prescription(f::Function) = Prescribed(f)
prescription(value) = Constant(value)


"""
    breakpoints(prescription)

Return the times at which `prescription` has a kink, sorted: the rows of a table, the
declared breakpoints of a function, and none for a constant.
"""
breakpoints(::Constant) = ()
breakpoints(prescription::Prescribed) = prescription.breakpoints
breakpoints(table::Tabulated) = Tuple(table.times)


"""
    depends_on_current_position(prescription)

Return `true` when `prescription` depends on the current position ``x = X + u``.
"""
depends_on_current_position(::AbstractPrescription) = false
depends_on_current_position(prescription::Prescribed) = prescription.current_position


"""
    prescribed_value(prescription, X, t; x = X)

Evaluate `prescription` at the reference position `X` and the time `t` (and, for a value
that depends on it, the current position `x`).
"""
prescribed_value(constant::Constant, X, t; x=X) = constant.value

function prescribed_value(prescription::Prescribed, X, t; x=X)
    # A value that depends on the current position is given it as a second argument.
    if prescription.current_position
        return prescription.f(X, x, t)
    else
        return prescription.f(X, t)
    end
end

function prescribed_value(table::Tabulated, X, t; x=X)
    times = table.times
    values = table.values

    # Outside the table the history is held at its first or last value.
    if t <= first(times)
        return first(values)
    elseif t >= last(times)
        return last(values)
    end

    # Inside, we interpolate linearly between the two rows that bracket t.
    upper = searchsortedfirst(times, t)
    lower = upper - 1
    θ = (t - times[lower])/(times[upper] - times[lower])
    value = interpolate_rows(values[lower], values[upper], θ)
    return value
end


"""
    interpolate_rows(a, b, θ)

Return ``(1 - θ) a + θ b`` for two rows of a table, componentwise for tuples.
"""
interpolate_rows(a::Number, b::Number, θ) = (1 - θ)*a + θ*b
interpolate_rows(a::Tuple, b::Tuple, θ) = map((ai, bi) -> (1 - θ)*ai + θ*bi, a, b)
interpolate_rows(a::AbstractVector, b::AbstractVector, θ) = @. (1 - θ)*a + θ*b


"""
    table_issues(table)

Return a message for every way in which `table` is malformed: no rows, a different number
of times and values, or times that do not strictly increase.
"""
function table_issues(table::Tabulated)
    messages = String[]

    # A table needs at least one row, and as many values as times.
    if isempty(table.times)
        push!(messages, "has no rows")
    end
    if length(table.times) != length(table.values)
        push!(
            messages,
            "has $(length(table.times)) times but $(length(table.values)) values",
        )
    end

    # Interpolation needs the times in strictly increasing order.
    if !all(i -> table.times[i] < table.times[i+1], 1:(length(table.times)-1))
        push!(messages, "has times that do not strictly increase")
    end
    return messages
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Value shapes -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractValueShape

Supertype of the shapes a prescribed value can be required to have.
"""
abstract type AbstractValueShape end


"""
    ScalarValue()

A single number.
"""
struct ScalarValue <: AbstractValueShape end


"""
    VectorValue(length)

A tuple or vector of numbers, of length `length` when that is known without the spatial
dimension, and of any length otherwise.

# Fields
- `length::Union{Int,Nothing}`: the required length, or `nothing`.
"""
struct VectorValue <: AbstractValueShape
    length::Union{Int,Nothing}
end


"""
    UnknownShape()

A shape that cannot be checked before the mesh is known (for example the value of a
tensor field).
"""
struct UnknownShape <: AbstractValueShape end


"""
    fits_shape(value, shape)

Return `true` when the constant `value` has the shape `shape`. A single number fits a
vector of length one, since one selected component is often written as a number.
"""
fits_shape(value, ::UnknownShape) = true
fits_shape(value, ::ScalarValue) = value isa Real
fits_shape(value, ::VectorValue) = false

function fits_shape(value::Union{Tuple,AbstractVector}, shape::VectorValue)
    # Every entry must be a number, and the length must match when it is known.
    entries_are_numbers = all(entry -> entry isa Real, value)
    length_matches = shape.length === nothing || length(value) == shape.length
    return entries_are_numbers && length_matches
end

fits_shape(value::Real, shape::VectorValue) = shape.length == 1


"""
    shape_label(shape)

Return a description of `shape` for messages.
"""
shape_label(::ScalarValue) = "a number"
shape_label(shape::VectorValue) =
    shape.length === nothing ? "a tuple of numbers" : "a tuple of $(shape.length) numbers"
shape_label(::UnknownShape) = "a value"


"""
    prescription_shape_issues(prescription, shape)

Return a message for every value of `prescription` that does not have the shape `shape`,
and for every way a table is malformed. A function is checked at compile time, when it can
be evaluated on the region.
"""
prescription_shape_issues(::Prescribed, ::AbstractValueShape) = String[]

function prescription_shape_issues(constant::Constant, shape::AbstractValueShape)
    messages = String[]
    if !fits_shape(constant.value, shape)
        push!(messages, "has value $(constant.value), but needs $(shape_label(shape))")
    end
    return messages
end

function prescription_shape_issues(table::Tabulated, shape::AbstractValueShape)
    # A table must be well formed, and every row must have the shape.
    messages = table_issues(table)
    for value in table.values
        if !fits_shape(value, shape)
            push!(messages, "has row $value, but needs $(shape_label(shape))")
        end
    end
    return messages
end
