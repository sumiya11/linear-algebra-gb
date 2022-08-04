

#------------------------------------------------------------------------------
####### Basis #######

# s-pair, a pair of polynomials
struct SPair
    # first generator as index from the basis array
    poly1::Int
    # second generator -//-
    poly2::Int
    # position of lcm(poly1, poly2) in hashtable
    lcm::Int
    # total degree of lcm
    deg::UInt
end

mutable struct Pairset
    pairs::Vector{SPair}
    # number of filled pairs,
    # Initially zero
    load::Int
end

function initialize_pairset(; initial_size=2^6) # TODO: why 64?
    pairs = Vector{SPair}(undef, initial_size)
    return Pairset(pairs, 0)
end

function Base.isempty(ps::Pairset)
    return ps.load == 0
end

function check_enlarge_pairset!(ps::Pairset, added::Int)
    sz = length(ps.pairs)
    # TODO: shrink pairset ?
    # @warn ps.load + added sz
    if ps.load + added >= sz
        # @error "resuze"  ps.load + added sz ps.load + added >= sz
        newsz = max(2 * sz, ps.load + added)
        resize!(ps.pairs, newsz)
    end
end

#------------------------------------------------------------------------------

mutable struct Basis{T}
    # vector of polynomials, each polynomial is a vector of monomials,
    # each monomial is represented with it's position in hashtable
    gens::Vector{Vector{Int}}
    # polynomial coefficients
    coeffs::Vector{Vector{T}}

    #= Keeping track of sizes   =#
    #=  ndone <= ntotal <= size =#
    # total allocated size,
    # size == length(gens) is always true
    size::Int
    # number of processed polynomials in `gens`
    # Iitially a zero
    ndone::Int
    # total number of polys filled in `gens`
    # (these will be handled during next update! call)
    # Iitially a zero
    ntotal::Int

    #= Keeping track of redundancy =#
    #= invariant =#
    #= length(lead) == length(nonred) == count(isred) == nlead =#
    # if element of the basis
    # is redundant
    isred::Vector{Int8}
    # positions of non-redundant elements in the basis
    nonred::Vector{Int}
    # division masks of leading monomials of
    # non redundant basis elements
    lead::Vector{UInt32}
    # number of filled elements in lead
    nlead::Int

    # characteristic of ground field
    ch::UInt64
end

function initialize_basis(ring::PolyRing, ngens::Int, ::Type{T}) where {T<:Coeff}
    #=
        always true
        length(gens) == length(coeffs) == length(isred) == size
    =#

    sz     = ngens * 2 # hmmm
    ndone  = 0
    ntotal = 0
    nlead  = 0

    gens   = Vector{Vector{Int}}(undef, sz)
    coeffs = Vector{Vector{T}}(undef, sz)
    isred  = zeros(Int8, sz)
    nonred = Vector{Int}(undef, sz)
    lead   = Vector{UInt32}(undef, sz)

    ch = ring.ch

    Basis(gens, coeffs, sz, ndone, ntotal, isred, nonred, lead, nlead, ch)
end

function initialize_basis(ring::PolyRing, hashedexps, coeffs::Vector{Vector{T}}) where {T<:Coeff}
    sz     = length(hashedexps) # hmmm
    ndone  = 0
    ntotal = 0
    nlead  = 0

    isred  = zeros(Int8, sz)
    nonred = Vector{Int}(undef, sz)
    lead   = Vector{UInt32}(undef, sz)

    ch = ring.ch

    Basis(hashedexps, coeffs, sz, ndone, ntotal, isred, nonred, lead, nlead, ch)
end

#------------------------------------------------------------------------------

function copy_basis_thorough(basis::Basis{T}) where {T}
    #  That cost a day of debugging ////
    gens   = Vector{Vector{Int}}(undef, basis.size)
    coeffs = Vector{Vector{T}}(undef, basis.size)
    @inbounds for i in 1:basis.ntotal
        gens[i] = Vector{Int}(undef, length(basis.gens[i]))
        coeffs[i] = Vector{T}(undef, length(basis.coeffs[i]))

        @inbounds for j in 1:length(basis.gens[i])
            gens[i][j] = basis.gens[i][j]
            coeffs[i][j] = basis.coeffs[i][j]
        end
    end
    isred  = copy(basis.isred)
    nonred = copy(basis.nonred)
    lead = copy(basis.lead)
    Basis(gens, coeffs, basis.size, basis.ndone,
            basis.ntotal, isred, nonred, lead,
            basis.nlead, basis.ch)
end

#------------------------------------------------------------------------------

function check_enlarge_basis!(basis::Basis{T}, added::Int) where {T}
    if basis.ndone + added >= basis.size
        basis.size = max(basis.size * 2, basis.ndone + added)
        resize!(basis.gens, basis.size)
        resize!(basis.coeffs, basis.size)
        resize!(basis.isred, basis.size)
        basis.isred[basis.ndone+1:end] .= 0
        resize!(basis.nonred, basis.size)
        resize!(basis.lead, basis.size)
    end
end

#------------------------------------------------------------------------------

# Normalize each element of the input basis
# by dividing it by leading coefficient
function normalize_basis!(basis::Basis{CoeffFF})
    cfs = basis.coeffs
    @inbounds for i in 1:basis.ntotal
        # mul = inv(cfs[i][1])
        # hack for now, TODODO
        if !isassigned(cfs, i)
            continue
        end
        ch = basis.ch
        mul = invmod(cfs[i][1], ch) % ch
        @inbounds for j in 2:length(cfs[i])
            # cfs[i][j] *= mul
            # TODO: faster division
            cfs[i][j] = cfs[i][j]*mul % ch
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end

# Normalize each element of the input basis
# by dividing it by leading coefficient
function normalize_basis!(basis::Basis{CoeffQQ})
    cfs = basis.coeffs
    @inbounds for i in 1:basis.ntotal
        # mul = inv(cfs[i][1])
        # hack for now, TODODO
        if !isassigned(cfs, i)
            continue
        end
        mul = inv(cfs[i][1])
        @inbounds for j in 2:length(cfs[i])
            cfs[i][j] *= mul
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end