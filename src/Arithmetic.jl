# The type of the numbers used in the computation
T = QQFieldElem;

# This line is useful to fix the arithmetic of the numbers. T is the type we use, and F is a function that returns the numbers 0 and 1 to the same type as T.
function F(n::Union{Int64, UInt16})::T
    return T(n)
end 

Marks_type = Vector{Int64}