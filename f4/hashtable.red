
%
%   struct Hashvalue;
%
%   stores hash of a monomial,
%   corresponding divmask to speed up divisibility checks,
%   index for position matrix (defaults to zero),
%   and the todal degree 
%
asserted procedure hashtable_Hashvalue(
        hash: ExponentHash, 
        divmask: DivisionMask, 
        idx: Integer, 
        deg: Degree): Hashvalue;
    {'hv, hash, divmask, idx, deg};

asserted procedure hashtable_hvget_hash(hv: Hashvalue): ExponentHash;
    cadr hv; 

asserted procedure hashtable_hvget_divmask(hv: Hashvalue): DivisionMask;
    caddr hv; 

asserted procedure hashtable_hvget_idx(hv: Hashvalue): Integer;
    cadddr hv; 

asserted procedure hashtable_hvget_deg(hv: Hashvalue): Degree;
    caddddr hv; 

asserted procedure hashtable_copy_hashvalue(hv: Hashvalue): Hashvalue;
    hashtable_Hashvalue(
        hashtable_hvget_hash(hv),
        hashtable_hvget_divmask(hv),
        hashtable_hvget_idx(hv),
        hashtable_hv_deg(hv)
    );

%--------------------------------------------------------------------------------------------------

# Returns the next look-up position in the table 
# (that is, implementing open addressing with quadratic probing) 
asserted procedure hashtable_hashnextindex(h: ExponentHash, j: ExponentHash, mod: ExponentHash): ExponentHash;
    land((h + j), mod) + 1;
end

%--------------------------------------------------------------------------------------------------

%   struct MonomialHashtable
%
% {'ht, exponents, hashtable, hashdata, hasher,
%   nvars, explen, ord, 
%   divmap, divvars, ndivvars, ndivbits,
%   size, load, offset}
%
%
asserted procedure hashtable_MonomialHashtable(
                exponents, hashtable, hashdata, hasher,nvars, explen, ord, 
                divmap, divvars, ndivvars, ndivbits, 
                size, load, offset): MonomialHashtable;
    {'ht, exponents, hashtable, hashdata, hasher,
        nvars, explen, ord, 
        divmap, divvars, ndivvars, ndivbits, 
        size, load, offset};

% initialize and set fields for basis hashtable
asserted procedure hashtable_initialize_basis_hashtable(ring: PolyRing): MonomialHashtable;
    begin scalar initial_size, exponents, hashdata, hashtable;
            integer nvars, explen; 
        initial_size := 16;

        exponents := dv_undef(initial_size);
        hashdata := dv_undef(initial_size);
        hashtable := dv_undef(initial_size);
        
        nvars := io_pr_nvars(ring);
        explen := io_pr_explen(ring);
        ord := io_pr_explen(ring);
        
        % initialize hashing vector
        hasher := dv_zeros(nvars);
        % fill hasher
        for i := 1:explen do
            % we don't want hash vector components to be zero
            while getv(hasher, i) = 0 do
                putv(hasher, i, random(sizeofInt32));
        
        % exponents array starts from index offset,
        % We store buffer array at index 1
        load := 1;
        size := initial_size;
        offset := 2;

        % initialize fast divisibility params
        ASSERT(sizeofInt32 = 32);
        ndivbits := sizeofInt32 / nvars;
        % division mask stores at least 1 bit
        % per each of first sizeofInt32 variables
        if ndivbits = 0 then
            ndivbits := ndivbits + 1;
        ndivvars := if nvars < int32bits then
            nvars
        else
            int32bits;
        divvars := dv_undef(ndivvars);
        divmap := dv_undef(ndivvars * ndivbits);
        % count only first ndivvars variables for divisibility checks
        for i := 1:ndivvars do
            putv(divvars, i, i);

        putv(exponents, 1, dv_zeros(explen));

        return hashtable_MonomialHashtable(
                exponents, hashtable, hashdata, hasher,
                nvars, explen, ord, 
                divmap, divvars, ndivvars, ndivbits, 
                size, load, offset)
    end;

asserted procedure hashtable_htget_exponents(ht: Hashtable): DynamicVector;
    car cdrn(ht, 1);

asserted procedure hashtable_htget_hashtable(ht: Hashtable): DynamicVector;
    car cdrn(ht, 2);

asserted procedure hashtable_htget_hashdata(ht: Hashtable): DynamicVector;
    car cdrn(ht, 3);

asserted procedure hashtable_htget_hasher(ht: Hashtable): DynamicVector;
    car cdrn(ht, 4);

asserted procedure hashtable_htget_nvars(ht: Hashtable): Integer;
    car cdrn(ht, 5);

asserted procedure hashtable_htget_explen(ht: Hashtable): Integer;
    car cdrn(ht, 6);

asserted procedure hashtable_htget_ord(ht: Hashtable);
    car cdrn(ht, 7);

asserted procedure hashtable_htget_divmap(ht: Hashtable);
    car cdrn(ht, 8);

asserted procedure hashtable_htget_divvars(ht: Hashtable);
    car cdrn(ht, 9);

asserted procedure hashtable_htget_ndivvars(ht: Hashtable);
    car cdrn(ht, 10);

asserted procedure hashtable_htget_ndivbits(ht: Hashtable);
    car cdrn(ht, 11);

asserted procedure hashtable_htget_size(ht: Hashtable): Integer;
    car cdrn(ht, 12);

asserted procedure hashtable_htget_load(ht: Hashtable): Integer;
    car cdrn(ht, 13);

asserted procedure hashtable_htget_offset(ht: Hashtable): Integer;
    car cdrn(ht, 14);

asserted procedure hashtable_htset_exponents(ht: Hashtable, x): DynamicVector;
    car cdrn(ht, 1) := x;

asserted procedure hashtable_htset_hashtable(ht: Hashtable, x): DynamicVector;
    car cdrn(ht, 2) := x;

asserted procedure hashtable_htset_hashdata(ht: Hashtable, x): DynamicVector;
    car cdrn(ht, 3) := x;

asserted procedure hashtable_htset_hasher(ht: Hashtable, x): DynamicVector;
    car cdrn(ht, 4) := x;

asserted procedure hashtable_htset_nvars(ht: Hashtable, x): Integer;
    car cdrn(ht, 5) := x;

asserted procedure hashtable_htset_explen(ht: Hashtable, x): Integer;
    car cdrn(ht, 6) := x;

asserted procedure hashtable_htset_ord(ht: Hashtable, x);
    car cdrn(ht, 7) := x;

asserted procedure hashtable_htset_divmap(ht: Hashtable, x);
    car cdrn(ht, 8) := x;

asserted procedure hashtable_htset_divvars(ht: Hashtable, x);
    car cdrn(ht, 9) := x;

asserted procedure hashtable_htset_ndivvars(ht: Hashtable, x);
    car cdrn(ht, 10) := x;

asserted procedure hashtable_htset_ndivbits(ht: Hashtable, x);
    car cdrn(ht, 11) := x;

asserted procedure hashtable_htset_size(ht: Hashtable, x): Integer;
    car cdrn(ht, 12) := x;

asserted procedure hashtable_htset_load(ht: Hashtable, x): Integer;
    car cdrn(ht, 13) := x;

asserted procedure hashtable_htset_offset(ht: Hashtable, x): Integer;
    car cdrn(ht, 14) := x;

asserted procedure hashtable_copy_hash_table(ht: MonomialHashtable): MonomialHashtable;
    begin scalar sz, load;
        load := hashtable_htget_load(ht);
        sz := hashtable_htget_size(ht);
        htexps := hashtable_htget_exponents(ht);
        httable := hashtable_htget_hashtable(ht);
        htdata := hashtable_htget_hashdata(ht);

        exps := dv_undef(sz);
        table := dv_undef(sz);
        data := dv_undef(sz);

        for i := 2:load do <<
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

%--------------------------------------------------------------------------------------------------

asserted procedure hashtable_initialize_secondary_hash_table(ht: MonomialHashtable): MonomialHashtable;
    begin scalar initial_size;
        initial_size := 8; % Julia 2^6
        
        exponents := dv_undef(sz);
        hashdata := dv_undef(sz);
        hashtable := dv_undef(sz);

        nvars := hashtable_htget_nvars(ht);
        explen := hashtable_htget_explen(ht);
        ord := hashtable_htget_ord(ht);

        divmap := hashtable_htget_divmap(ht);
        divvars := hashtable_htget_divvars(ht);
        ndivvars := hashtable_htget_ndivvars(ht);
        ndivbits := hashtable_htget_ndivbits(ht);

        hasher := hashtable_htget_hasher(ht);
        load := 1;
        size := initial_size;
        offset := 2;

        putv(exponents, 1, dv_zeros(explen));
        
        return  hashtable_MonomialHashtable(
                exponents, hashtable, hashdata, hasher,
                nvars, explen, ord,
                divmap, divvars, ndivvars, ndivbits,
                size, load, offset)
    end;

%--------------------------------------------------------------------------------------------------

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
    begin scalar x;
        htsize := hashtable_htget_size(ht) * 2;
        hashtable_htset_size(ht, htsize);
        hashtable_htset_hashdata(ht, dv_resize(hashtable_htget_hashdata(ht), htsize));
        hashtable_htset_exponents(ht, dv_resize(hashtable_htget_exponents(ht), htsize))

        hashtable_htset_hashtable(ht, dv_zeros(htsize));

        htdata := hashtable_htget_hashdata(ht);
        httable := hashtable_htget_hashtable(ht);

        mod = htsize - 1;
        for i := hashtable_htget_offset(ht):hashtable_htget_load(ht) do <<
            % hash for this elem is already computed
            he := hashtable_hvget_hash(getv(htdata, i));
            hidx := he;
            for j := 1:htsize do <<
                hidx := hashtable_hashnextindex(he, j, mod);
                if getv(httable, hidx) = 0 then <<
                    setv(httable, hidx, i);
                    j := htsize
                >>
            >>
        >>
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure hashtable_insert_in_hash_table(ht: MonomialHashtable, e::ExponentVector);
    begin scalar he;
        he := 0;
        hasher := hashtable_htget_hasher(ht);
        htexplen := hashtable_htget_explen(ht);
        for i := 1:htexplen 
            he := he + getv(hasher, i)*getv(e, i);
        
        % find new elem position in the table
        htsize := hashtable_htget_size(ht);
        hidx := he;
        mod = htsize - 1;
        i = 1;

        httable := hashtable_htget_hashtable(ht);
        htdata := hashtable_htget_hashdata(ht);
        htexps := hashtable_htget_exponents(ht);

    Restart:
        while i < htsize do <<
            hidx := hashtable_hashnextindex(he, i, mod);
            vidx := getv(httable, hidx);
            % if not free
            if not (vidx = 0) then <<
                if not (hashtable_hvget_hash(getv(htdata, vidx)) = he) then
                    i := i + 1
                else <<
                    present := getv(htexps, vidx);
                    for j in 1:htexplen do <<
                        if not (getv(present, j) = getv(e, j)) then <<
                            i := i + 1;
                            go to Restart
                        >>
                    >>
                >>
            >>;

            return vidx
        >>;

        vidx := hashtable_htget_load(ht) + 1;
        putv(httable, hidx, vidx);

        putv(htexps, vidx, copy(e));
        ve := getv(htexps, vidx);
        for i := 1:length(e) do
            putv(ve, i, getv(e, i));
        
        divmask := hashtable_generate_monomial_divmask(e, ht);
        putv(htdata, vidx, hashtable_Hashvalue(he, divmask, 0, getvlast(e)));

        hashtable_htset_load(ht, hashtable_htget_load(ht) + 1);

        return vidx
    end;

%--------------------------------------------------------------------------------------------------

# Having `ht` filled with monomials from input polys,
# computes ht.divmap and divmask for each of already stored monomials
asserted procedure hashtable_fill_divmask(ht: MonomialHashtable);
    begin scalar x;
        ndivvars := hashtable_htget_ndivvars(ht);
        divvars := hashtable_htget_divvars(ht);
        ndivbits := hashtable_htget_ndivvars(ht);
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
        for i := 1:ndivvars <<
            steps := (getv(max_exp, i) - getv(min_exp, i)) / ndivbits;
            if steps = 0 then
                steps := steps + 1;
            for j := 1:ndivbits do <<
                putv(divmap, ctr, steps);
                steps := steps + 1;
                ctr := ctr + 1
            >>
        >>;

        for vidx := hashtable_htset_offset(ht):hashtable_htget_load(ht) do <<
            unmasked := getv(htdata, vidx);
            e := getv(htexps, vidx);
            divmask := hashtable_generate_monomial_divmask(e, ht);
            putv(htdata, vidx, hashtable_Hashvalue(hashtable_hvget_hash(unmasked), divmask, 0, getvlast(e)))
        >>;

    end;

asserted procedure hashtable_generate_monomial_divmask(e: ExponentVector, ht: MonomialHashtable): DivisionMask;
    begin scalar ctr;
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

%--------------------------------------------------------------------------------------------------

asserted procedure hashtable_is_monom_divisible(h1: ExponentIdx, h2: ExponentIdx, ht: MonomialHashtable);
    begin scalar e1;
        htdata := hashtable_htget_hashdata(ht);

        % if fast_division_check

        htexps := hashtable_htget_exponents(ht);
        e1 := getv(htexps, h1);
        e2 := getv(htexps, h2);

        toreturn := t;
        for i := 1:hashtable_htget_explen(ht) do <<
            if getv(e1, i) < getv(e2, i) then
                toreturn := nil
        >>;

        return toreturn
    end;

asserted procedure hashtable_is_gcd_const(h1: ExponentIdx, h2: ExponentIdx, ht: MonomialHashtable);
    begin scalar e1;
        htexps := hashtable_htget_exponents(ht);
        e1 := getv(htexps, h1);
        e2 := getv(htexps, h2);

        toreturn := t;
        for i := 1:hashtable_htget_explen(ht) do <<
            if (not (getv(e1, i) = 0)) and (not (getv(e2, i) = 0)) then
                toreturn := nil
        >>;

        return toreturn
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure hashtable_check_monomial_division_in_update(
                a: Vector, first: Integer, last: Integer,
                lcmj: ExponentIdx, ht: MonomialHashtable);
    begin scalar divmask;
        htdata := hashtable_htget_hashdata(ht);
        htexps := hashtable_htget_exponents(ht);

        divmask := hashtable_hvget_divmask(getv(htdata, lcmj));
        lcmexp := getv(htexps, lcmj);
        
        j := first;
    Restart:
        while j <= last do <<
            if getv(a, j) = 0 then
                j := j + 1
            else <<
                % if fast_division_check
                ea := getv(htexps, getv(a, j));
                for j := 1:hashtable_htget_explen(ht) do <<
                    if getv(ea, i) < getv(lcmexp, i) then <<
                          j := j + 1;
                          go to Restart
                    >>
                >>
            >>;
            putv(a, j, 0)
        >>
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure hashtable_insert_multiplied_poly_in_hash_table(
            row: Vector, htmp: ExponentHash, etmp: ExponentVector, 
            poly: Vector, ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar len;
        len := length(poly);
        explen := hashtable_htget_explen(ht);

        mod := hashtable_htget_size(symbol_ht) - 1;

        bexps := hashtable_htget_exponents(ht);
        bdata := hashtable_htget_hashdata(ht);

        sexps := hashtable_htget_exponents(symbol_ht);
        sdata := hashtable_htget_hashdata(symbol_ht);
        stable := hashtable_htget_hashtable(symbol_ht);

        l := 1;
    Letsgo:
        while l <= len do <<
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

            envew := getv(sexps, 1);
            
            for j := 1:explen do 
                putv(enew, getv(etmp, j) + getv(e, j));

            k := h;
            i := 1;

        Restart:
            while i <= hashtable_htget_size(symbol_ht) do <<
                k = hashtable_hashnextindex(h, i, mod);
                vidx := getv(httable, k);
                if vidx = 0 then
                    i := hashtable_htget_size(symbol_ht)
                else <<
                    if not (hashtable_hvget_hash(getv(htdata, vidx)) = h) then
                        i := i + 1
                    else <<
                        estored := getv(sexps, vidx);
                        for j := 1:explen do <<
                            if not (getv(estored, j) = getv(enew, j)) then <<
                                i := i + 1;
                                go to Restart
                            >>
                        >>
                    >>
                >>;
                putv(row, l, vidx);
                l := l + 1;
                go to Letsgo
            >>;
            if null getv(sexps, lastidx) then
                putv(sexps, lastidx, dv_undef(explen));
            sexpsnew := getv(sexps, lastidx);
            for j := 1:explen do
                putv(sexpsnew, j, getv(enew, j));
            putv(stable, k, lastidx);
            divmask := hashtable_generate_monomial_divmask(enew, symbol_ht);
            putv(sdata, lastidx, hashtable_Hashvalue(h, divmask, 0, getvlast(enew)));
            putv(row, l, lastidx);
            l := l + 1;
            hashtable_htset_load(symbol_ht, hashtable_htget_load(symbol_ht) + 1);
        >>;

        return row
    end;

% Julia: invalid parameter "load";
% solution: switched from name "load" to "l"
asserted procedure hashtable_symbolic_ht_needscale(l: Integer, added: Integer, size: Integer);
    1.4*(l + added) >= size;

asserted procedure hashtable_multiplied_poly_to_matrix_row(
                symbolic_ht: MonomialHashtable, basis_ht: MonomialHashtable, 
                htmp: ExponentHash, etmp: ExponentVector, poly: MonomsVector);
    begin scalar row;
        row := dv_similar(poly);
        htload := hashtable_htget_load(symbol_ht);
        while hashtable_symbolic_ht_needscale(htload, length(poly), hashtable_htget_size(symbolic_ht)) do
            hashtable_enlarge_hash_table(symbolic_ht);
        
        return hashtable_insert_multiplied_poly_in_hash_table(row, htmp, etmp, poly, basis_htm symbolic_ht)
    end;

%--------------------------------------------------------------------------------------------------


%--------------------------------------------------------------------------------------------------

% computes lcm of he1 and he2 as exponent vectors from ht1
% and inserts it in ht2
asserted procedure hashtable_get_lcm(he1: ExponentIdx, he2: ExponentIdx,
                                    ht1: MonomialHashtable, ht2: MonomialHashtable);
    begin scalar e1;
        htexps := hashtable_htget_exponents(ht1);
        e1 := getv(htexps, h1);
        e2 := getv(htexps, h2);
        etmp := getv(htexps, 1);

        putv(etmp, length(etmp), 0);
        for i := 1:hashtable_htget_explen(ht1)-1 do <<
            putv(etmp, i, max(getv(e1, i), getv(e2, i)));
            putv(etmp, length(etmp), getvlast(etmp) + getv(etmp, i))
        >>;

        return hashtable_insert_in_hash_table(ht2, etmp)
    end;

%--------------------------------------------------------------------------------------------------
