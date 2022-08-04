
#------------------------------------------------------------------------------
####### Hashtable ######

# TODO: make immutable
#=
    stores hash of a monomial,
    corresponding divmask to speed up divisibility checks,
    index for position matrix (defaults to zero),
    and the todal degree
=#
mutable struct Hashvalue
    hash::UInt32

    #=
    =#
    divmask::UInt32

    idx::Int
    deg::Int
end

function copy_hashvalue(x::Hashvalue)
    Hashvalue(x.hash, x.divmask, x.idx, x.deg)
end


#=
    Hashtable implementing linear probing
    and designed to store and operate with monomials
=#
mutable struct MonomialHashtable
    exponents::Vector{ExponentVector}

    # maps exponent hash to its position in exponents array
    hashtable::Vector{Int}

    # stores hashes, division masks,
    # and other valuable info
    # for each hashtable enrty
    hashdata::Vector{Hashvalue}

    # values to hash exponents with, i.e
    # hash(e) = sum(hasher .* e)
    hasher::Vector{UInt32}

    #= Ring information =#
    # number of variables
    nvars::Int
    # raw length of exponent vector
    explen::Int
    # ring monomial ordering
    ord::Symbol

    #= Divisibility =#
    # divisor map to check divisibility faster
    divmap::Vector{UInt32}
    # variables selected for divmap
    divvars::Vector{Int}
    # count of divmap variables
    ndivvars::Int
    # bits per div variable
    ndivbits::Int

    size::Int
    # elements added
    load::Int
    #
    offset::Int
end

#------------------------------------------------------------------------------

function hashnextindex(h::UInt32, j::UInt32, mod::UInt32)
     (h + j) & mod + UInt32(1)
end

#------------------------------------------------------------------------------

# initialize and set fields for basis hashtable
function initialize_basis_hash_table(
        ring::PolyRing,
        rng::Random.AbstractRNG;
        initial_size::Int=2^16)

    exponents = Vector{Vector{UInt16}}(undef, initial_size)
    hashdata  = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(Int, initial_size)

    nvars = ring.nvars
    explen = ring.explen
    ord = ring.ord

    # initialize hashing vector
    hasher = zeros(UInt32, explen)
    for i in 1:explen
        # we don't want hash vector components to be zero
        while hasher[i] == 0
            hasher[i] = rand(rng, UInt32)
        end
    end

    # exponents[1:load] cover all stored exponents
    # , also exponents[1] is zeroed by default
    load = 1
    size = initial_size

    # exponents array starts from index offset,
    # We store buffer array at index 1
    offset = 2

    # initialize fast divisibility params
    charbit    = 8 # TODO ??
    int32bits  = charbit * sizeof(Int32)
    int32bits != 32 && error("Strange story with ints")
    ndivbits   = div(int32bits, nvars)
    # division mask stores at least 1 bit
    # per each of first charbit*sizeof(Int32) variables
    ndivbits == 0 && (ndivbits += 1)
    ndivvars = nvars < int32bits ? nvars : int32bits
    divvars  = Vector{Int}(undef, ndivvars)
    divmap   = Vector{UInt32}(undef, ndivvars * ndivbits)
    # count only first ndivvars variables for divisibility checks
    for i in 1:ndivvars
        divvars[i] = i
    end

    # first stored exponent used as buffer lately
    exponents[1] = zeros(UInt16, explen)

    MonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, explen, ord,
        divmap, divvars, ndivvars, ndivbits,
        size, load, offset)
end

function copy_hashtable(ht::MonomialHashtable)
    exps = Vector{Vector{UInt16}}(undef, ht.size)
    table = Vector{Int}(undef, ht.size)
    data = Vector{Hashvalue}(undef, ht.size)
    exps[1] = zeros(UInt16, ht.explen)

    @inbounds for i in 2:ht.load
        # TODO: ELIDE COPY
        exps[i] = copy(ht.exponents[i])
        table[i] = ht.hashtable[i]
        data[i] = copy_hashvalue(ht.hashdata[i])
    end

    MonomialHashtable(
        ht.exponents, ht.hashtable, ht.hashdata, ht.hasher,
        ht.nvars, ht.explen, ht.ord,
        ht.divmap, ht.divvars, ht.ndivvars, ht.ndivbits,
        ht.size, ht.load, ht.offset)
end

#------------------------------------------------------------------------------

# initialize hashtable either for `symbolic_preprocessing` or for `update` functions
# These are of the same purpose and structure as basis hashtable,
# but are more local oriented
function initialize_secondary_hash_table(basis_ht::MonomialHashtable)

    # 2^6 seems to be the best out of 2^5, 2^6, 2^7
    initial_size = 2^6

    exponents = Vector{Vector{UInt16}}(undef, initial_size)
    hashdata  = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(Int, initial_size)

    # preserve ring info
    explen = basis_ht.explen
    nvars  = basis_ht.nvars
    ord    = basis_ht.ord

    # preserve division info
    divmap   = basis_ht.divmap
    divvars  = basis_ht.divvars
    ndivvars = basis_ht.ndivvars
    ndivbits = basis_ht.ndivbits

    # preserve hasher
    hasher = basis_ht.hasher

    load = 1
    size = initial_size
    offset = 2

    exponents[1] = zeros(UInt16, explen)

    MonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, explen, ord,
        divmap, divvars, ndivvars, ndivbits,
        size, load, offset)
end

#------------------------------------------------------------------------------

# resizes (if needed) ht so that it can store `size` elements,
# and clears all previoud data
function reinitialize_hash_table!(ht::MonomialHashtable, size::Int)
    if size > ht.size
        while size > ht.size
            ht.size *= 2
        end
        # TODO
        resize!(ht.hashdata, ht.size)
        resize!(ht.exponents, ht.size)
    end
    ht.hashtable = zeros(Int, ht.size)
    ht.load = 1
end

# doubles the size of storage in `ht`,
# and rehashes all elements
function enlarge_hash_table!(ht::MonomialHashtable)
    ht.size *= 2
    resize!(ht.hashdata,  ht.size)
    resize!(ht.exponents, ht.size)

    ht.hashtable = zeros(Int, ht.size)

    mod = UInt32(ht.size - 1)
    for i in ht.offset:ht.load
        # hash for this elem is already computed
        he = ht.hashdata[i].hash
        hidx = he
        @inbounds for j in UInt32(1):UInt32(ht.size)
            hidx = hashnextindex(he, j, mod)
            ht.hashtable[hidx] != 0 && continue
            ht.hashtable[hidx] = i
            break
        end
    end
end