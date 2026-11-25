#=
Pull-backs and push-forwards.

Under finite strains a quantity measured per unit current volume or area must be pulled
back before it can be integrated over the reference configuration, and a law written in
spatial terms (isotropic Fourier conduction, ``q = -k ∇_x θ``) needs its inputs pushed
forward. Both are nodes in the graph, because both carry a derivative with respect to the
deformation gradient: the pull-back ``Q = J F^{-1} q`` couples a heat flux to the
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
    PullBack(node; motion = nothing)

The pull-back of `node` to the reference configuration.

# Fields
- `node::N`: the node being pulled back.
- `motion::Union{Symbol,Nothing}`: the displacement field defining the motion, or
  `nothing` to use the motion declared by [`FiniteStrain`](@ref).
"""
struct PullBack{N<:AbstractNode} <: AbstractNode
    node::N
    motion::Union{Symbol,Nothing}
end

PullBack(node::AbstractNode; motion::Union{Symbol,Nothing} = nothing) =
    PullBack(node, motion)


"""
    PushForward(node; motion = nothing)

The push-forward of `node` to the current configuration.

# Fields
- `node::N`: the node being pushed forward.
- `motion::Union{Symbol,Nothing}`: the displacement field defining the motion, or
  `nothing` to use the motion declared by [`FiniteStrain`](@ref).
"""
struct PushForward{N<:AbstractNode} <: AbstractNode
    node::N
    motion::Union{Symbol,Nothing}
end

PushForward(node::AbstractNode; motion::Union{Symbol,Nothing} = nothing) =
    PushForward(node, motion)


const AbstractTransformNode = Union{PullBack,PushForward}


# ---------------------------------------------------------------------------------------- #
# ---------------------------------- Graph bookkeeping ----------------------------------- #
# ---------------------------------------------------------------------------------------- #

function quantity(transform::PullBack, fields::Tuple)
    # An unknown inner quantity, or one without a reference counterpart, gives no tag.
    inner = quantity(transform.node, fields)
    if inner === nothing
        return nothing
    else
        return pulled_back(inner)
    end
end


function quantity(transform::PushForward, fields::Tuple)
    # An unknown inner quantity, or one without a current counterpart, gives no tag.
    inner = quantity(transform.node, fields)
    if inner === nothing
        return nothing
    else
        return pushed_forward(inner)
    end
end


"""
    transform_motion(transform, assumption)

Return the displacement field whose motion `transform` uses under `assumption`, or
`nothing` under small strains, where the transform is the identity.
"""
transform_motion(::AbstractTransformNode, ::SmallStrain) = nothing

function transform_motion(transform::AbstractTransformNode, assumption::FiniteStrain)
    # A transform may name its own motion; otherwise it uses the problem's.
    if transform.motion === nothing
        return motion_field(assumption)
    else
        return transform.motion
    end
end


function node_fields(
    transform::AbstractTransformNode,
    assumption::AbstractKinematicAssumption,
)
    # The transformed value depends on the inner node and, under finite strains, on F.
    fields = copy(node_fields(transform.node, assumption))
    motion = transform_motion(transform, assumption)
    if motion !== nothing && !(motion in fields)
        push!(fields, motion)
    end
    return fields
end


# Under small strains a transform is the identity; under finite strains it involves F and
# J, so it is never affine.
is_affine(transform::AbstractTransformNode, assumption::SmallStrain) =
    is_affine(transform.node, assumption)
is_affine(::AbstractTransformNode, ::FiniteStrain) = false

subnodes(transform::AbstractTransformNode) = (transform.node,)


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
``\\mathrm{d}v = J \\, \\mathrm{d}V_0`` for volumes, and Nanson's formula for areas. Returns
`nothing` for measures that are already referential.
"""
measure_transformation_rule(::AbstractMeasure) = nothing
measure_transformation_rule(::CurrentVolume) = VolumeDensityRule()
measure_transformation_rule(::CurrentBoundary) = NansonRule()
