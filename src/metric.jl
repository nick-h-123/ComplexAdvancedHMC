abstract type AbstractMetric end

_string_M⁻¹(mat::AbstractMatrix, n_chars::Int=32) = _string_M⁻¹(diag(mat), n_chars)
function _string_M⁻¹(vec::AbstractVector, n_chars::Int=32)
    s_vec = string(vec)
    l = length(s_vec)
    s_dots = " ...]"
    n_diag_chars = n_chars - length(s_dots)
    return s_vec[1:min(n_diag_chars,end)] * (l > n_diag_chars ? s_dots : "")
end

struct UnitEuclideanMetric{T,A<:Union{Tuple{Int},Tuple{Int,Int}}} <: AbstractMetric
    M⁻¹::UniformScaling{T}
    size::A
end

UnitEuclideanMetric(::Type{T}, sz) where {T} = UnitEuclideanMetric(UniformScaling{T}(one(T)), sz)
UnitEuclideanMetric(sz) = UnitEuclideanMetric(Float64, sz)
UnitEuclideanMetric(::Type{T}, dim::Int) where {T} = UnitEuclideanMetric(UniformScaling{T}(one(T)), (dim,))
UnitEuclideanMetric(dim::Int) = UnitEuclideanMetric(Float64, (dim,))

renew(ue::UnitEuclideanMetric, M⁻¹) = UnitEuclideanMetric(M⁻¹, ue.size)

Base.size(e::UnitEuclideanMetric) = e.size
Base.size(e::UnitEuclideanMetric, dim::Int) = e.size[dim]
Base.show(io::IO, uem::UnitEuclideanMetric) = print(io, "UnitEuclideanMetric($(_string_M⁻¹(ones(uem.size))))")

struct EuclideanMetric{T,A<:AbstractVecOrMat{T}} <: AbstractMetric
    # Diagnal of the inverse of the mass matrix
    M⁻¹     ::  A
    # Sqare root of the inverse of the mass matrix
    sqrtM⁻¹ ::  A
    # Pre-allocation for intermediate variables
    _temp   ::  A
end

function EuclideanMetric(M⁻¹::AbstractVecOrMat{T}) where {T<:AbstractFloat}
    return EuclideanMetric(M⁻¹, sqrt.(M⁻¹), similar(M⁻¹))
end
EuclideanMetric(::Type{T}, sz) where {T} = EuclideanMetric(ones(T, sz...))
EuclideanMetric(sz) = EuclideanMetric(Float64, sz)
EuclideanMetric(::Type{T}, dim::Int) where {T} = EuclideanMetric(ones(T, dim))
EuclideanMetric(dim::Int) = EuclideanMetric(Float64, dim)

renew(ue::EuclideanMetric, M⁻¹) = EuclideanMetric(M⁻¹)

Base.size(e::EuclideanMetric, dim...) = size(e.M⁻¹, dim...)
Base.show(io::IO, dem::EuclideanMetric) = print(io, "EuclideanMetric($(_string_M⁻¹(dem.M⁻¹)))")


struct HermitianMetric{T,A<:AbstractVecOrMat{T}} <: AbstractMetric
    # Diagnal of the inverse of the mass matrix
    M⁻¹     ::  A
    # Sqare root of the inverse of the mass matrix
    sqrtM⁻¹ ::  A
    # Pre-allocation for intermediate variables
    _temp   ::  A
end

function HermitianMetric(M⁻¹::AbstractVecOrMat{T}) where {T<:Complex}
    return HermitianMetric(M⁻¹, sqrt.(M⁻¹), similar(M⁻¹))
end
HermitianMetric(::Type{T}, sz) where {T} = HermitianMetric(ones(T, sz...))
HermitianMetric(sz) = HermitianMetric(Float64, sz)
HermitianMetric(::Type{T}, dim::Int) where {T} = HermitianMetric(ones(T, dim))
HermitianMetric(dim::Int) = HermitianMetric(Float64, dim)

renew(ue::HermitianMetric, M⁻¹) = HermitianMetric(M⁻¹)

Base.size(e::HermitianMetric, dim...) = size(e.M⁻¹, dim...)
Base.show(io::IO, dem::HermitianMetric) = print(io, "HermitianMetric($(_string_M⁻¹(dem.M⁻¹)))")

# getname functions
for T in (UnitEuclideanMetric, EuclideanMetric, HermitianMetric)
    @eval getname(::Type{<:$T}) = $T
end
getname(m::T) where {T<:AbstractMetric} = getname(T)

# `rand` functions for `metric` types.

function _rand(
    rng::Union{AbstractRNG, AbstractVector{<:AbstractRNG}},
    metric::UnitEuclideanMetric{T}
) where {T}
    r = randn(rng, T, size(metric)...)
    return r
end

function _rand(
    rng::Union{AbstractRNG, AbstractVector{<:AbstractRNG}},
    metric::EuclideanMetric{T}
) where {T}
    r = randn(rng, T, size(metric)...)
    r ./= metric.sqrtM⁻¹
    return r
end

function _rand(
    rng::Union{AbstractRNG, AbstractVector{<:AbstractRNG}},
    metric::HermitianMetric{T}
) where {T<:Complex}
    r = randn(rng, T, size(metric)...)
    if size(metric.M⁻¹) == 1
        r = 0.5*(r + conj(reverse(test)))
    elseif size(metric.M⁻¹) == 2
        r = 0.5*(r + conj(transpose(test)))
    end
    r ./= metric.sqrtM⁻¹
    return r
end

Base.rand(rng::AbstractRNG, metric::AbstractMetric) = _rand(rng, metric)    # this disambiguity is required by Random.rand
Base.rand(rng::AbstractVector{<:AbstractRNG}, metric::AbstractMetric) = _rand(rng, metric)
Base.rand(metric::AbstractMetric) = rand(GLOBAL_RNG, metric)
