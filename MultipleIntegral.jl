module IterMultipleIntegrals
export integral, faster_integral, product_integral
using Base: product

# range iterators product (hypervolume) implementation
# consistent performance in high dimensions
function product_integral(f, intervals...; N=100)
    """
    Example: G = {x₁⁴ + x₂⁴ + x₃⁴ <= 1}, volume G = ∫dμG
    q(x...) = +((x.^4)...) <= 1
    multiple_integral(q, [(-1,1) for _ in 1:3]...; N=200)
    """
    n = length(intervals)
    if n == 0
        return f()
    end
    steps = map(x->x[2]-x[1], intervals) ./ N
    d = *(steps...)
    partition = product((range(a + step/2, b - step/2, length=N)
            for ((a, b), step) in zip(intervals, steps))...)

    # summing f on intervals iterator product: no mallocs
    sum(partition) do x
        f(x...) * d
    end
end

alwaystrue(args...) = true
# accurate recursive (iterated integral) implementation
# but lots of mallocs in cycles
function integral(f, (a, b); indomain=alwaystrue, N=300)
    s = 0
    step = (b-a)/N
    for x in range(a+step/2, b-step/2, length=(N-1))
        s += indomain(x) && f(x)
    end
    return s * step
end

function integral(f, intervals...; indomain=alwaystrue, N=100)
    s = 0
    (a, b), intervals... = intervals
    step = (b-a)/N
    for x in range(a+step/2, b-step/2, length=(N-1))
        fproj = ((args...,) -> f(x, args...))
        indomainproj = ((args...,) -> indomain(x, args...))
        s += integral(fproj, intervals...; indomain=indomainproj, N=N)
    end
    return s * step
end

# recursive implementation, faster in small cases, less accurate
function faster_integral(f, (a, b); indomain=alwaystrue, N=300)
    step = (b-a)/N
    sum(x-> indomain(x) && f(x), range(a, b, length=(N))) * step
end

function faster_integral(f, intervals...; indomain=alwaystrue, N=100)
    (a, b), intervals... = intervals
    step = (b-a)/N
    sum(range(a, b, length=(N))) do x
        fproj = ((args...,) -> f(x, args...))
        indomainproj = ((args...,) -> indomain(x, args...))
        faster_integral(fproj, intervals...; indomain=indomainproj, N=N)
    end * step
end

end
# f(x,y) = (1<=(x*y)<=3) && (x<=y<=(2*x)) && y^2
# Is = ((0.5,2),(1,3))
# # equivalent
# g(x,y) = y^2
# gdom(x,y) = (1<=(x*y)<=3) && (x<=y<=(2*x))

# q(x...) = +((x.^4)...) <= 1
# Iq(n) = [(-1,1) for _ in 1:n]
# Iq3 = Iq(3);
# volume(fq, [(-1,1) for _ in 1:3]...; N=200)
