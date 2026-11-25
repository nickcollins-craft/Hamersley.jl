#=
Errors shared across the package.

Hamersley is being built structure-first: many functions exist so that the shape of the
package is fixed before the numerics behind them are written. Those functions throw a
`NotYetImplemented` error, so that a reader (or a test) can tell "this is a deliberate
placeholder" apart from a genuine bug.
=#


"""
    NotYetImplemented(feature)

Exception thrown by functions whose interface is fixed but whose implementation has not
been written yet.

# Fields
- `feature::String`: a short description of the missing capability, phrased so that it
  completes the sentence "… is not yet implemented in Hamersley".
"""
struct NotYetImplemented <: Exception
    feature::String
end


function Base.showerror(io::IO, exception::NotYetImplemented)
    # We phrase the message so that the feature description reads as the subject.
    print(
        io,
        "NotYetImplemented: ",
        exception.feature,
        " is not yet implemented in Hamersley",
    )
    return nothing
end


"""
    not_yet_implemented(feature)

Throw a [`NotYetImplemented`](@ref) error describing `feature`.

This is the single entry point that placeholder functions use, so that every placeholder
fails in the same, recognisable way.
"""
function not_yet_implemented(feature::AbstractString)
    throw(NotYetImplemented(String(feature)))
end
