#=
Pull-backs and push-forwards.

Under finite strains a quantity measured per unit current volume or area must be pulled
back before it can be integrated over the reference configuration, and a law written in
spatial terms (for example isotropic Fourier conduction, ``q = -k ∇_x T``) needs its inputs pushed
forward. Both are vertices in the graph, because both carry a derivative with respect to the
deformation gradient: the pull-back ``q_0 = J F^{-1} q`` couples a heat flux to the
displacement.

Users never name the transformation itself. The rule (Piola, covariant, Nanson, …) is a
property of the quantity's tag (see `transformation_rule` in `core/quantities.jl`). Under
small strains every pull-back and push-forward is the identity, so a formulation written
with them runs unchanged under either assumption.

Measures are pulled back too: an integral over the current volume becomes one over the
reference volume with a factor ``J``, and one over the current boundary picks up Nanson's
formula.
=#


"""
    PullBack(vertex; motion = nothing)

The pull-back of `vertex` to the reference configuration.

# Fields
- `vertex::N`: the vertex being pulled back.
- `motion::Union{Symbol,Nothing}`: the displacement field defining the motion, or
  `nothing` to use the motion declared by [`FiniteStrain`](@ref).
"""
struct PullBack{N<:AbstractVertex} <: AbstractVertex
    vertex::N
    motion::Union{Symbol,Nothing}
end

PullBack(vertex::AbstractVertex; motion::Union{Symbol,Nothing}=nothing) =
    PullBack(vertex, motion)


"""
    PushForward(vertex; motion = nothing)

The push-forward of `vertex` to the current configuration.

# Fields
- `vertex::N`: the vertex being pushed forward.
- `motion::Union{Symbol,Nothing}`: the displacement field defining the motion, or
  `nothing` to use the motion declared by [`FiniteStrain`](@ref).
"""
struct PushForward{N<:AbstractVertex} <: AbstractVertex
    vertex::N
    motion::Union{Symbol,Nothing}
end

PushForward(vertex::AbstractVertex; motion::Union{Symbol,Nothing}=nothing) =
    PushForward(vertex, motion)


const AbstractTransformVertex = Union{PullBack,PushForward}


# ---------------------------------------------------------------------------------------- #
# ---------------------------------- Graph bookkeeping ----------------------------------- #
# ---------------------------------------------------------------------------------------- #

# A transform carries the tag of its inner quantity's counterpart. An unknown inner
# quantity, or one without a counterpart in the target configuration, gives no tag.
quantity(transform::PullBack, fields::Tuple) =
    wrapped_quantity(pulled_back, transform, fields)
quantity(transform::PushForward, fields::Tuple) =
    wrapped_quantity(pushed_forward, transform, fields)


"""
    transform_motion(transform, assumption)

Return the displacement field whose motion `transform` uses under `assumption`, or
`nothing` under small strains, where the transform is the identity.
"""
transform_motion(::AbstractTransformVertex, ::SmallStrain) = nothing

function transform_motion(transform::AbstractTransformVertex, assumption::FiniteStrain)
    # A transform may name its own motion; otherwise it uses the problem's.
    if transform.motion === nothing
        return motion_field(assumption)
    else
        return transform.motion
    end
end


function vertex_fields(
    transform::AbstractTransformVertex,
    assumption::AbstractKinematicAssumption,
)
    # The transformed value depends on the inner vertex and, under finite strains, on F.
    fields = copy(vertex_fields(transform.vertex, assumption))
    motion = transform_motion(transform, assumption)
    if motion !== nothing && !(motion in fields)
        push!(fields, motion)
    end
    return fields
end


# Under small strains a transform is the identity; under finite strains it involves F and
# J, so it is never affine.
is_affine(transform::AbstractTransformVertex, assumption::SmallStrain) =
    is_affine(transform.vertex, assumption)
is_affine(::AbstractTransformVertex, ::FiniteStrain) = false

subvertices(transform::AbstractTransformVertex) = (transform.vertex,)


# ---------------------------------------------------------------------------------------- #
# ---------------------------------- Measure pull-backs ---------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    pulled_back_measure(measure)

Return the reference-configuration counterpart of `measure`. Under the total Lagrangian
approach every integral is ultimately taken over ``Ω_0`` or ``∂Ω_0``.
"""
pulled_back_measure(measure::AbstractMeasure) = measure
pulled_back_measure(::CurrentVolume) = ReferenceVolume()
pulled_back_measure(::CurrentBoundary) = ReferenceBoundary()


"""
    measure_transformation_rule(measure)

Return the rule by which `measure` is mapped onto its reference counterpart:
``\\mathrm{d}v = J \\, \\mathrm{d}v_0`` for volumes, and Nanson's formula for areas. Returns
`nothing` for measures that are already referential.
"""
measure_transformation_rule(::AbstractMeasure) = nothing
measure_transformation_rule(::CurrentVolume) = VolumeDensityRule()
measure_transformation_rule(::CurrentBoundary) = NansonRule()


"""
    measure_fields(measure, assumption)

Return a vector of the names of the fields that `measure` depends on under `assumption`.

A reference measure is fixed, so it depends on nothing. Under [`FiniteStrain`](@ref) a
current measure is pulled back with a factor that depends on the motion (``J`` for a volume,
``J F^{-T} n_0`` for an area), so an integral over it depends on the motion even when its
integrand does not. Under [`SmallStrain`](@ref) the configurations coincide, so no measure
depends on anything.
"""
measure_fields(::AbstractMeasure, ::AbstractKinematicAssumption) = Symbol[]

function measure_fields(measure::AbstractMeasure, assumption::FiniteStrain)
    # Only a measure that is pulled back with a rule picks up the motion.
    if measure_transformation_rule(measure) === nothing
        return Symbol[]
    else
        return [motion_field(assumption)]
    end
end


"""
    is_affine(measure, assumption)

Return `true` when integrating against `measure` keeps an affine integrand affine under
`assumption`.

Under [`FiniteStrain`](@ref) a current measure is not: ``∫ δT \\, c T \\, \\mathrm{d}v =
∫ δT \\, c T J \\, \\mathrm{d}v_0`` is a product of ``T`` and ``J(u)``, and even a constant
body force per unit current volume gives a residual in ``J(u)``.
"""
is_affine(::AbstractMeasure, ::AbstractKinematicAssumption) = true
is_affine(measure::AbstractMeasure, assumption::FiniteStrain) =
    isempty(measure_fields(measure, assumption))
