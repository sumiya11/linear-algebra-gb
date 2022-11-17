module f4hashtable;
% Monomial hashtable implementation.

%   struct Hashvalue;
%
%   stores hash of a monomial,
%   corresponding divmask to speed up divisibility checks,
%   index for position matrix (defaults to zero),
%   and the todal degree 
%
% the struct is mutable, getters and setters are defined below
asserted procedure hashtable_Hashvalue(
        hash: ExponentHash, 
        divmask: DivisionMask, 
        idx: Integer, 
        deg: Degree): Hashvalue;
    begin scalar v;
        v := dv_undef(4);
        % v[1] := hash;
        putv(v, 1, hash);
        putv(v, 2, divmask);
        putv(v, 3, idx);
        putv(v, 4, deg);
        return v
    end;

asserted procedure hashtable_hvget_hash(hv: Hashvalue): ExponentHash;
    getv(hv, 1); 

asserted procedure hashtable_hvget_divmask(hv: Hashvalue): DivisionMask;
    getv(hv, 2); 

asserted procedure hashtable_hvget_idx(hv: Hashvalue): Integer;
    getv(hv, 3); 

asserted procedure hashtable_hvget_deg(hv: Hashvalue): Degree;
    getv(hv, 4); 

asserted procedure hashtable_hvset_hash(hv: Hashvalue, x: ExponentHash): ExponentHash;
    putv(hv, 1, x); 

asserted procedure hashtable_hvset_divmask(hv: Hashvalue, x: DivisionMask): DivisionMask;
    putv(hv, 2, x); 

asserted procedure hashtable_hvset_idx(hv: Hashvalue, x: Integer): Integer;
    putv(hv, 3, x); 

asserted procedure hashtable_hvset_deg(hv: Hashvalue, x: Degree): Degree;
    putv(hv, 4, x); 

asserted procedure hashtable_copy_hashvalue(hv: Hashvalue): Hashvalue;
    hashtable_Hashvalue(
        hashtable_hvget_hash(hv),
        hashtable_hvget_divmask(hv),
        hashtable_hvget_idx(hv),
        hashtable_hv_deg(hv)
    );

%   struct MonomialHashtable
%   
% MonomialHashtable implementing linear probing
% and designed to store and operate with monomials.
% 
% Stores the following:
% . exponents - a vector of exponent vectors (monomials),
% . hashtable - a vector, maps exponent hash to its position
%               in the exponents array,
% . hashdata - a vector, stores hashes, division masks, 
%              and other valuable info for each valid `hashtable` enrty,
% . hasher - a vector, stores values to hash exponents with, i.e
%            hash(e) = sum(hasher .* e),
% . nvars - an integer, number of variables,
% . explen - an integer, raw length of exponent vector,
% . ord - ring monomial ordering (in some representation),
% . divmap - divisor map to check monom-by-monom divisibility faster,
% . divvars - variables selected for divmap,
% . ndivvars - count of divmap variables,
% . ndivbits - bits per div variables
% . size - total size (capacity) of hashtable
% . loadj - number of elements added (always <= size)
% . offset
asserted procedure hashtable_MonomialHashtable(
                exponents: Vector, hashtable: Vector, hashdata: Vector, 
                hasher: Vector, nvars: Integer, explen: Integer, ord: Any, 
                divmap: Vector, divvars: Vector, ndivvars: Integer, ndivbits: Integer, 
                size: Integer, loadj: Integer, offset: Integer): MonomialHashtable;
    begin scalar v;
        v := dv_undef(14);
        putv(v, 1, exponents);
        putv(v, 2, hashtable);
        putv(v, 3, hashdata);
        putv(v, 4, hasher);
        putv(v, 5, nvars);
        putv(v, 6, explen);
        putv(v, 7, ord);
        putv(v, 8, divmap);
        putv(v, 9, divvars);
        putv(v, 10, ndivvars);
        putv(v, 11, ndivbits);
        putv(v, 12, size);
        putv(v, 13, loadj);
        putv(v, 14, offset);
        return v
    end;

% Getters
asserted procedure hashtable_htget_exponents(ht: MonomialHashtable): Vector;
    getv(ht, 1);

asserted procedure hashtable_htget_hashtable(ht: MonomialHashtable): Vector;
    getv(ht, 2);

asserted procedure hashtable_htget_hashdata(ht: MonomialHashtable): Vector;
    getv(ht, 3);

asserted procedure hashtable_htget_hasher(ht: MonomialHashtable): Vector;
    getv(ht, 4);

asserted procedure hashtable_htget_nvars(ht: MonomialHashtable): Integer;
    getv(ht, 5);

asserted procedure hashtable_htget_explen(ht: MonomialHashtable): Integer;
    getv(ht, 6);

asserted procedure hashtable_htget_ord(ht: MonomialHashtable): Any;
    getv(ht, 7);

asserted procedure hashtable_htget_divmap(ht: MonomialHashtable): Vector;
    getv(ht, 8);

asserted procedure hashtable_htget_divvars(ht: MonomialHashtable): Vector;
    getv(ht, 9);

asserted procedure hashtable_htget_ndivvars(ht: MonomialHashtable): Integer;
    getv(ht, 10);

asserted procedure hashtable_htget_ndivbits(ht: MonomialHashtable): Integer;
    getv(ht, 11);

asserted procedure hashtable_htget_size(ht: MonomialHashtable): Integer;
    getv(ht, 12);

asserted procedure hashtable_htget_load(ht: MonomialHashtable): Integer;
    getv(ht, 13);

asserted procedure hashtable_htget_offset(ht: MonomialHashtable): Integer;
    getv(ht, 14);

% Setters
asserted procedure hashtable_htset_exponents(ht: MonomialHashtable, x: Vector): Vector;
    putv(ht, 1, x);

asserted procedure hashtable_htset_hashtable(ht: MonomialHashtable, x: Vector): Vector;
    putv(ht, 2, x);

asserted procedure hashtable_htset_hashdata(ht: MonomialHashtable, x: Vector): Vector;
    putv(ht, 3, x);

asserted procedure hashtable_htset_hasher(ht: MonomialHashtable, x: Vector): Vector;
    putv(ht, 4, x);

asserted procedure hashtable_htset_nvars(ht: MonomialHashtable, x: Integer): Integer;
    putv(ht, 5, x);

asserted procedure hashtable_htset_explen(ht: MonomialHashtable, x: Integer): Integer;
    putv(ht, 6, x);

asserted procedure hashtable_htset_ord(ht: MonomialHashtable, x: Any): Any;
    putv(ht, 7, x);

asserted procedure hashtable_htset_divmap(ht: MonomialHashtable, x: Vector): Vector;
    putv(ht, 8, x);

asserted procedure hashtable_htset_divvars(ht: MonomialHashtable, x: Vector): Vector;
    putv(ht, 9, x);

asserted procedure hashtable_htset_ndivvars(ht: MonomialHashtable, x: Integer): Integer;
    putv(ht, 10, x);

asserted procedure hashtable_htset_ndivbits(ht: MonomialHashtable, x: Integer): Integer;
    putv(ht, 11, x);

asserted procedure hashtable_htset_size(ht: MonomialHashtable, x: Integer): Integer;
    putv(ht, 12, x);

asserted procedure hashtable_htset_load(ht: MonomialHashtable, x: Integer): Integer;
    putv(ht, 13, x);

asserted procedure hashtable_htset_offset(ht: MonomialHashtable, x: Integer): Integer;
    putv(ht, 14, x);

% Julia: variable name "mod" changed to "modj"
% Returns the next look-up position in the table 
% (that is, implementing open addressing with quadratic probing) 
asserted procedure hashtable_hashnextindex(h: ExponentHash, j: ExponentHash, modj: ExponentHash): ExponentHash;
    land(h + j, modj) + 1;

% initialize and set fields for basis hashtable
asserted procedure hashtable_initialize_basis_hash_table(ring: PolyRing, initial_size: Integer): MonomialHashtable;
    begin scalar initial_size, exponents, hashdata, hashtable,
                    nvars, explen, ord, hasher, loadj, size, offset,
                    ndivbits, ndivvars, divvars, divmap; 

        exponents := dv_undef(initial_size);
        hashdata := dv_undef(initial_size);
        hashtable := dv_zeros(initial_size);
        
        nvars := io_prget_nvars(ring);
        explen := io_prget_explen(ring);
        ord := io_prget_ord(ring);
        
        % initialize hashing vector
        hasher := dv_zeros(explen);
        % fill hasher
        for i := 1:explen do
            % we don't want hash vector components to be zero
            while hasher[i] = 0 do
                hasher[i] := random(2^f4_sizeofInt32!*);
        
        % exponents array starts from index offset,
        % We store buffer array at index 1
        loadj := 1;
        size := initial_size;

        % exponents array starts from index offset,
        % We store buffer array at index 1
        offset := 2;

        % initialize fast divisibility params
        ASSERT(f4_sizeofInt32!* = 32);
        ndivbits := f4_sizeofInt32!* / nvars;
        % division mask stores at least 1 bit
        % per each of first f4_sizeofInt32!* variables
        if ndivbits = 0 then
            ndivbits := ndivbits + 1;
        ndivvars := if nvars < f4_sizeofInt32!* then
            nvars
        else
            f4_sizeofInt32!*;
        divvars := dv_undef(ndivvars);
        divmap := dv_zeros(ndivvars * ndivbits);
        % count only first ndivvars variables for divisibility checks
        for i := 1:ndivvars do
            divvars[i] := i;

        % first stored exponent used as buffer lately
        exponents[1] := dv_zeros(explen);

        return hashtable_MonomialHashtable(
                exponents, hashtable, hashdata, hasher,
                nvars, explen, ord, 
                divmap, divvars, ndivvars, ndivbits, 
                size, loadj, offset)
    end;

asserted procedure hashtable_copy_hashtable(ht: MonomialHashtable): MonomialHashtable;
    begin scalar loadj, sz, htexps, httable, htdata, exps, table, data;
        loadj := hashtable_htget_load(ht);
        sz := hashtable_htget_size(ht);
        htexps := hashtable_htget_exponents(ht);
        httable := hashtable_htget_hashtable(ht);
        htdata := hashtable_htget_hashdata(ht);

        exps := dv_undef(sz);
        table := dv_undef(sz);
        data := dv_undef(sz);

        for i := 2:loadj do <<
            putv(exps, i, getv(htexps, i));
            putv(table, i, getv(httable, i));
            putv(data, i, getv(htdata, i))
        >>;

        return hashtable_MonomialHashtable(
            htexps, httable, htdata, hashtable_htget_hasher(ht),
            hashtable_htget_nvars(ht), hashtable_htget_explen(ht), hashtable_htget_ord(ht),
            hashtable_htget_divmap(ht), hashtable_htget_divvars(ht), hashtable_htget_ndivvars(ht), 
            hashtable_htget_ndivbits(ht),
            hashtable_htget_size(ht), hashtable_htget_load(ht), hashtable_htget_offset(ht)
        )
    end;

% initialize hashtable either for `symbolic_preprocessing` or for `update` functions
% These are of the same purpose and structure as basis hashtable,
% but are more local oriented
asserted procedure hashtable_initialize_secondary_hash_table(ht: MonomialHashtable): MonomialHashtable;
    begin scalar initial_size, exponents, hashdata, hashtable,
                    explen, nvars, ord, divmap, divvars, ndivvars, ndivbits,
                    hasher, loadj, size, offset;
        initial_size := 2^6; % best size in Julia is 2^6
        
        exponents := dv_undef(initial_size);
        hashdata := dv_undef(initial_size);
        hashtable := dv_zeros(initial_size);

        % preserve ring info
        nvars := hashtable_htget_nvars(ht);
        explen := hashtable_htget_explen(ht);
        ord := hashtable_htget_ord(ht);

        % preserve division info
        divmap := hashtable_htget_divmap(ht);
        divvars := hashtable_htget_divvars(ht);
        ndivvars := hashtable_htget_ndivvars(ht);
        ndivbits := hashtable_htget_ndivbits(ht);

        % preserve hasher
        hasher := hashtable_htget_hasher(ht);

        loadj := 1;
        size := initial_size;
        offset := 2;

        putv(exponents, 1, dv_zeros(explen));
        
        return hashtable_MonomialHashtable(
                exponents, hashtable, hashdata, hasher,
                nvars, explen, ord,
                divmap, divvars, ndivvars, ndivbits,
                size, loadj, offset)
    end;

% resizes (if needed) ht so that it can store `size` elements,
% and clears all previoud data
asserted procedure hashtable_reinitialize_hash_table(ht: MonomialHashtable, size: Integer);
    begin scalar size, htsize;
        htsize := hashtable_htget_size(ht);
        if size > htsize then <<
            while size > htsize do
                htsize := htsize * 2;
            hashtable_htset_size(ht, htsize);
            hashtable_htset_hashdata(ht, dv_resize(hashtable_htget_hashdata(ht), htsize));
            hashtable_htset_exponents(ht, dv_resize(hashtable_htget_exponents(ht), htsize))
        >>;
        hashtable_htset_hashtable(ht, dv_zeros(htsize));
        hashtable_htset_load(ht, 1)
    end;

% doubles the size of storage in `ht`,
% and rehashes all elements
asserted procedure hashtable_enlarge_hash_table(ht: MonomialHashtable);
    begin scalar htsize, modj, he, hidx, htdata, httable, j;
        htsize := hashtable_htget_size(ht) * 2;
        hashtable_htset_size(ht, htsize);
        % 3 vectors of equal size inside ht: hashdata (add), exponents (add), hashtable (table)
        hashtable_htset_hashdata(ht, dv_resize(hashtable_htget_hashdata(ht), htsize));
        hashtable_htset_exponents(ht, dv_resize(hashtable_htget_exponents(ht), htsize));

        hashtable_htset_hashtable(ht, dv_zeros(htsize));

        htdata := hashtable_htget_hashdata(ht);
        httable := hashtable_htget_hashtable(ht);

        modj := htsize - 1;
        for i := hashtable_htget_offset(ht) : hashtable_htget_load(ht) do <<
            % hash for this elem is already computed
            he := hashtable_hvget_hash(htdata[i]);
            hidx := he;
            j := 0;
            repeat <<
                j := j + 1;
                ASSERT(j <= htsize);
                hidx := hashtable_hashnextindex(he, j, modj)
            >> until httable[hidx] = 0;
            httable[hidx] := i
        >>    
    end;

asserted procedure hashtable_insert_in_hash_table(ht: MonomialHashtable, e: ExponentVector): ExponentIdx;
    begin scalar he, hasher, htexplen, htsize, hidx, modj, i, httable, htdata, htexps,
                    vidx, present, ve, divmask, toreturn;
        % generate hash
        he := 0;

        hasher := hashtable_htget_hasher(ht);
        htexplen := hashtable_htget_explen(ht);
        for i := 1:htexplen do
            he := he + getv(hasher, i)*getv(e, i);

        % find new elem position in the table
        htsize := hashtable_htget_size(ht);
        hidx := he;
        modj := htsize - 1;
        i := 1;

        httable := hashtable_htget_hashtable(ht);
        htdata := hashtable_htget_hashdata(ht);
        htexps := hashtable_htget_exponents(ht);

        % Julia: can not return from cycle
        toreturn := 0;

    Restart:
        while i < htsize do << % begin
            hidx := hashtable_hashnextindex(he, i, modj);
            vidx := getv(httable, hidx);

            % if not free
            if vidx = 0 then
                go to Break;

            if hashtable_hvget_hash(getv(htdata, vidx)) neq he then <<
                % if not free and not same hash
                i := i + 1;
                go to Restart
            >>;
            
            present := getv(htexps, vidx);
            for j := 1:htexplen do <<
                if getv(present, j) neq getv(e, j) then <<
                    i := i + 1;
                    go to Restart
                >>
            >>;

            % already present in hashtable
            go to Return_
        >>;
    Break:
        
        % add its position to hashtable, and insert exponent to that position
        vidx := hashtable_htget_load(ht) + 1;
        putv(httable, hidx, vidx);

        putv(htexps, vidx, dv_similar(e));
        ve := getv(htexps, vidx);
        for i := 1 : dv_length(e) do
            putv(ve, i, getv(e, i));
        
        divmask := hashtable_generate_monomial_divmask(e, ht);
        putv(htdata, vidx, hashtable_Hashvalue(he, divmask, 0, f4_getvlast(e)));

        hashtable_htset_load(ht, hashtable_htget_load(ht) + 1);

    Return_:
        return vidx
    end;

% Having `ht` filled with monomials from input polys,
% computes ht.divmap and divmask for each of already stored monomials
asserted procedure hashtable_fill_divmask(ht: MonomialHashtable);
    begin scalar ndivvars, divvars, ndivbits, divmap, min_exp, max_exp,
                    htexps, htdata, e, ctr, steps, unmasked, divmask;
        ndivvars := hashtable_htget_ndivvars(ht);
        divvars := hashtable_htget_divvars(ht);
        ndivbits := hashtable_htget_ndivbits(ht);
        divmap := hashtable_htget_divmap(ht);

        min_exp := dv_undef(ndivvars);
        max_exp := dv_undef(ndivvars);

        htexps := hashtable_htget_exponents(ht);
        htdata := hashtable_htget_hashdata(ht);
        e := getv(htexps, hashtable_htget_offset(ht));

        for i := 1:ndivvars do <<
            putv(min_exp, i, getv(e, getv(divvars, i)));
            putv(max_exp, i, getv(e, getv(divvars, i)))
        >>;

        for i := hashtable_htget_offset(ht):hashtable_htget_load(ht) do <<
            e := getv(htexps, i);
            for j := 1:ndivvars do <<
                if getv(e, getv(divvars, j)) > getv(max_exp, j) then
                    putv(max_exp, j, getv(e, getv(divvars, j)));
                if getv(e, getv(divvars, j)) < getv(min_exp, j) then
                    putv(min_exp, j, getv(e, getv(divvars, j)))
            >>
        >>;

        ctr := 1;
        steps := 0;
        for i := 1:ndivvars do <<
            steps := (getv(max_exp, i) - getv(min_exp, i)) / ndivbits;
            if steps = 0 then
                steps := steps + 1;
            for j := 1:ndivbits do <<
                putv(divmap, ctr, steps);
                steps := steps + 1;
                ctr := ctr + 1
            >>
        >>;

        for vidx := hashtable_htget_offset(ht):hashtable_htget_load(ht) do <<
            unmasked := getv(htdata, vidx);
            e := getv(htexps, vidx);
            divmask := hashtable_generate_monomial_divmask(e, ht);
            putv(htdata, vidx, hashtable_Hashvalue(hashtable_hvget_hash(unmasked), divmask, 0, f4_getvlast(e)))
        >>

    end;

asserted procedure hashtable_generate_monomial_divmask(e: ExponentVector, ht: MonomialHashtable): DivisionMask;
    begin scalar divvars, divmap, ctr, res;
        divvars := hashtable_htget_divvars(ht);
        divmap := hashtable_htget_divmap(ht);

        ctr := 1;
        res := 0;

        for i := 1:hashtable_htget_ndivvars(ht) do <<
            for j := 1:hashtable_htget_ndivbits(ht) do <<
                if getv(e, getv(divvars, i)) >= getv(divmap, ctr) then
                    res := lor(res, lshift(1, ctr - 1));
                ctr := ctr + 1
            >>
        >>;

        return res
    end;

% h1 divisible by h2
asserted procedure hashtable_is_monom_divisible(h1: ExponentIdx, h2: ExponentIdx, ht: MonomialHashtable): Boolean;
    begin scalar htdata, e1, e2, result, htexps;
        htdata := hashtable_htget_hashdata(ht);

        if not (land(hashtable_hvget_divmask(getv(htdata, h2)), lnot(hashtable_hvget_divmask(getv(htdata, h1)))) = 0) then
            return nil;

        htexps := hashtable_htget_exponents(ht);
        e1 := getv(htexps, h1);
        e2 := getv(htexps, h2);

        % Julia
        result := t;
        for i := 1:hashtable_htget_explen(ht) do
            if getv(e1, i) < getv(e2, i) then <<
                result := nil;
                go to Return_
            >>;

    Return_:
        return result
    end;

asserted procedure hashtable_is_gcd_const(h1: ExponentIdx, h2: ExponentIdx, ht: MonomialHashtable): Boolean;
    begin scalar htexps, e1, e2, result;
        htexps := hashtable_htget_exponents(ht);
        e1 := getv(htexps, h1);
        e2 := getv(htexps, h2);

        result := t;
        for i := 1:hashtable_htget_explen(ht)-1 do
            if (not (getv(e1, i) = 0)) and (not (getv(e2, i) = 0)) then <<
                result := nil;
                go to Return_
            >>;

    Return_:
        return result
    end;

% compare pairwise divisibility of lcms from a[first:last] with lcm
asserted procedure hashtable_check_monomial_division_in_update(
                a: Vector, first: Integer, last: Integer,
                lcmj: ExponentIdx, ht: MonomialHashtable);
    begin scalar htdata, htexps, divmask, lcmexp, j, ea;
        htdata := hashtable_htget_hashdata(ht);
        htexps := hashtable_htget_exponents(ht);

        % pairs are sorted, we only need to check entries above starting point

        divmask := hashtable_hvget_divmask(getv(htdata, lcmj));
        lcmexp := getv(htexps, lcmj);
        
        j := first;
    Restart:
        while j <= last do <<
            % if bad lcm
            if getv(a, j) = 0 then <<
                j := j + 1;
                go to Restart
            >>;
            
            % fast division check
            if land(lnot(hashtable_hvget_divmask(getv(htdata, getv(a, j)))), divmask) neq 0 then <<
                j := j + 1;
                go to Restart
            >>;

            ea := getv(htexps, getv(a, j));
            for i := 1:hashtable_htget_explen(ht) do <<
                if getv(ea, i) < getv(lcmexp, i) then <<
                    j := j + 1;
                    go to Restart
                >>
            >>;
            % mark as redundant
            putv(a, j, 0)
        >>;

    end;

% add monomials from `poly` multiplied by exponent vector `etmp`
% with hash `htmp` to hashtable `symbol_ht`,
% and substitute hashes in row
asserted procedure hashtable_insert_multiplied_poly_in_hash_table(
            row: Vector, htmp: ExponentHash, etmp: ExponentVector, 
            poly: Vector, ht: MonomialHashtable, symbol_ht: MonomialHashtable): Vector;
    begin scalar len, explen, modj, bexps, bdata, btable, sexps, sdata, stable,
                    l, h, e, lastidx, enew, k, i, vidx, estored, sexpsnew, divmask;

        % length of poly to add
        len := dv_length(poly);
        explen := hashtable_htget_explen(ht);

        modj := hashtable_htget_size(symbol_ht) - 1;

        bexps := hashtable_htget_exponents(ht);
        bdata := hashtable_htget_hashdata(ht);
        btable := hashtable_htget_hashtable(ht);

        sexps := hashtable_htget_exponents(symbol_ht);
        sdata := hashtable_htget_hashdata(symbol_ht);
        stable := hashtable_htget_hashtable(symbol_ht);

        l := 1;
    Letsgo:
        while l <= len do begin
            % we iterate over all monoms of the given poly,
            % multiplying them by htmp/etmp,
            % and inserting into symbolic hashtable

            % hash is linear, so that
            % hash(e1 + e2) = hash(e1) + hash(e2)
            % We also assume that the hashing vector is shared same
            % between all created hashtables
            h := htmp + hashtable_hvget_hash(getv(bdata, getv(poly, l)));
            
            e := getv(bexps, getv(poly, l));

            lastidx := hashtable_htget_load(symbol_ht) + 1;
            enew := getv(sexps, 1);
            
            % multiplied monom exponent
            for j := 1:explen do 
                putv(enew, j, getv(etmp, j) + getv(e, j));

            % now insert into hashtable
            k := h;

            i := 1;

        Restart:
            while i <= hashtable_htget_size(symbol_ht) do <<
                k := hashtable_hashnextindex(h, i, modj);

                vidx := getv(stable, k);
                % if index is free
                if vidx = 0 then
                    go to Break;

                % if different exponent is stored here
                if hashtable_hvget_hash(getv(sdata, vidx)) neq h then <<
                    i := i + 1;
                    go to Restart
                >>;

                estored := getv(sexps, vidx);
                for j := 1:explen do <<
                    % hash collision, restarting the search
                    if getv(estored, j) neq getv(enew, j) then <<
                        i := i + 1;
                        go to Restart
                    >>
                >>;
                
                putv(row, l, vidx);
                l := l + 1;

                go to Letsgo
            >>;
        Break:

            % add multiplied exponent to hash table
            if null getv(sexps, lastidx) then
                putv(sexps, lastidx, dv_undef(explen));
            sexpsnew := getv(sexps, lastidx);
            % multiplied monom exponent
            for j := 1:explen do
                putv(sexpsnew, j, getv(enew, j));
            putv(stable, k, lastidx);

            divmask := hashtable_generate_monomial_divmask(enew, symbol_ht);
            putv(sdata, lastidx, hashtable_Hashvalue(h, divmask, 0, f4_getvlast(enew)));
            
            putv(row, l, lastidx);
            l := l + 1;
            hashtable_htset_load(symbol_ht, hashtable_htget_load(symbol_ht) + 1)
        end;

        return row
    end;

% Julia: invalid parameter "load";
% solution: switched from name "load" to "loadj"
%
% If symbolic hash table of the given size and load factor should be
% enlarged after adding the polynomial of length `added` 
asserted procedure hashtable_symbolic_ht_needscale(loadj: Integer, added: Integer, size: Integer): Boolean;
    1.4*(loadj + added) >= size;

asserted procedure hashtable_multiplied_poly_to_matrix_row(
                symbol_ht: MonomialHashtable, basis_ht: MonomialHashtable, 
                htmp: ExponentHash, etmp: ExponentVector, poly: MonomsVector): Vector;
    begin scalar row, htload;

        row := dv_similar(poly);
        htload := hashtable_htget_load(symbol_ht);
        while hashtable_symbolic_ht_needscale(htload, dv_length(poly), hashtable_htget_size(symbol_ht)) do
            hashtable_enlarge_hash_table(symbol_ht);
        
        return hashtable_insert_multiplied_poly_in_hash_table(row, htmp, etmp, poly, basis_ht, symbol_ht)
    end;

asserted procedure hashtable_insert_in_basis_hash_table_pivots(
                row: Vector, 
                ht: MonomialHashtable, 
                symbol_ht: MonomialHashtable, 
                col2hash: Vector);
    begin scalar sdata, sexps, modj, explen, bdata, bexps, bhash, 
                    l, hidx, h, lastidx, e, k, i, hm, pos, ehm;
        
        while hashtable_htget_size(ht) - hashtable_htget_load(ht) <= dv_length(row) do
            hashtable_enlarge_hash_table(ht);

        sdata := hashtable_htget_hashdata(symbol_ht);
        sexps := hashtable_htget_exponents(symbol_ht);

        modj := hashtable_htget_size(ht) - 1;
        explen := hashtable_htget_explen(ht);
        bdata := hashtable_htget_hashdata(ht);
        bexps := hashtable_htget_exponents(ht);
        bhash := hashtable_htget_hashtable(ht);

        l := 1;
    Letsgo:
        while l <= dv_length(row) do begin
            hidx := getv(col2hash, getv(row, l));

            % symbolic hash
            h := hashtable_hvget_hash(getv(sdata, hidx));

            lastidx := hashtable_htget_load(ht) + 1;
            putv(bexps, lastidx, getv(sexps, hidx));
            e := getv(bexps, lastidx);

            k := h;
            i := 1;
        Restart:
            while i <= hashtable_htget_size(ht) do <<
                k := hashtable_hashnextindex(h, i, modj);
                hm := getv(bhash, k);

                if hm = 0 then 
                    go to Break;

                if hashtable_hvget_hash(getv(bdata, hm)) neq h then
                    i := i + 1
                else <<
                    ehm := getv(bexps, hm);
                    for j := 1:explen do
                        if not (getv(e, j) = getv(ehm, j)) then <<
                            i := i + 1;
                            go to Restart
                        >>;
                    putv(row, l, hm);
                    l := l + 1;
                    go to Letsgo
                >>
            >>;
        Break:

            pos := lastidx;
            putv(bhash, k, pos);
            putv(row, l, pos);
            l := l + 1;

            putv(bdata, pos, hashtable_Hashvalue(h, hashtable_hvget_divmask(getv(sdata, hidx)), hashtable_hvget_idx(getv(sdata, hidx)), hashtable_hvget_deg(getv(sdata, hidx))));

            hashtable_htset_load(ht, hashtable_htget_load(ht) + 1)
        end

    end;

% computes lcm of he1 and he2 as exponent vectors from ht1
% and inserts it in ht2
asserted procedure hashtable_get_lcm(he1: ExponentIdx, he2: ExponentIdx,
                                    ht1: MonomialHashtable, ht2: MonomialHashtable);
    begin scalar e1, e2, etmp, htexps;

        htexps := hashtable_htget_exponents(ht1);
        
        e1 := getv(htexps, he1);
        e2 := getv(htexps, he2);
        etmp := getv(htexps, 1);

        putv(etmp, dv_length(etmp), 0);
        for i := 1:hashtable_htget_explen(ht1)-1 do <<
            putv(etmp, i, max(getv(e1, i), getv(e2, i)));
            putv(etmp, dv_length(etmp), f4_getvlast(etmp) + getv(etmp, i))
        >>;

        return hashtable_insert_in_hash_table(ht2, etmp)
    end;

endmodule; % end of hashtable module

end;