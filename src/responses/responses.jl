abstract type AbstractResponse end

response_quantity(response::AbstractResponse) = response.quantity
driving_quantity(response::AbstractResponse) = response.driving
response_measure(response::AbstractResponse) = response.measure
response_dependencies(response::AbstractResponse) = response.dependencies
response_law(response::AbstractResponse) = response.law

function evaluate(response::AbstractResponse, state, qp)
    throw(MethodError(evaluate, (response, state, qp)))
end

function linearize(response::AbstractResponse, state, qp)
    throw(MethodError(linearize, (response, state, qp)))
end

struct ResponseIssue
    message::String
end

Base.show(io::IO, issue::ResponseIssue) = print(io, issue.message)
