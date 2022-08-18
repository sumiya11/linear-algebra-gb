module f4f4;
% The F4 implementation.

%--------------------------------------------------------------------------------------------------

% Julia: the comment in Julia is not correct.
asserted procedure f4_reduction(ring: PolyRing, basis: Basis, matrixj: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    <<
        matrix_convert_hashes_to_columns(matrixj, symbol_ht);

        sorting_sort_matrix_rows_decreasing(matrixj); % for pivots,  AB part
        sorting_sort_matrix_rows_increasing(matrixj); % for reduced, CD part

        matrix_linear_algebra(ring, matrixj, basis);

        matrix_convert_matrix_rows_to_basis_elements(matrixj, basis, ht, symbol_ht)
    >>;

%--------------------------------------------------------------------------------------------------

% Initializes Basis and MonomialHashtable structures,
% fills input data from exponents and coeffs
%
% Hashtable initial size is set to tablesize
asserted procedure f4_initialize_structures(ring: PolyRing, exponents: Vector, 
                                            coeffs: Vector, tablesize: Integer);
    begin scalar basis, basis_ht;
        % basis for storing basis elements,
        % hashtable for hashing monomials occuring in the basis
        basis := basis_initialize_basis(ring, dv_length(exponents));
        basis_ht := hashtable_initialize_basis_hash_table(ring, tablesize);
        
        if f4_debug() then <<
            prin2t {"initialize_structures.."};
            prin2t basis;
            prin2t basis_ht
        >>;

        % filling the basis and hashtable with the given inputs
        basis_fill_data(basis, basis_ht, exponents, coeffs);
        
        if f4_debug() then <<
            prin2t {"initialize_structures: filled"};
            prin2t basis;
            prin2t basis_ht
        >>;

        hashtable_fill_divmask(basis_ht);

        if f4_debug() then <<
            prin2t {"initialize_structures: divmasks"};
            prin2t hashtable_htget_hashdata(basis_ht)
        >>;

        sorting_sort_gens_by_lead_increasing(basis, basis_ht);

        if f4_debug() then <<
            prin2t {"initialize_structures: after sort by lead"};
            prin2t basis
        >>;

        basis_normalize_basis(ring, basis);
        
        if f4_debug() then <<
            prin2t {"initialize_structures: after normalize"};
            prin2t basis
        >>;

        return basis . basis_ht
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure f4_reducegb(ring: PolyRing, basis: Basis, matrixj: MacaulayMatrix,
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar exponents, etmp, uprows, nonred, gens, up2coef, symdata, nrows,
                    k, lead, htdata, ntotal, i;
        
        exponents := hashtable_htget_exponents(ht);
        etmp := getv(exponents, 1);
        for i := 1:upbv(etmp) do
            putv(etmp, i, 0);
        %  etmp is now set to zero, and has zero hash

        matrix_reinitialize_matrix(matrixj, basis_bget_nlead(basis));
        uprows := matrix_mget_uprows(matrixj);

        nonred := basis_bget_nonred(basis);
        gens := basis_bget_gens(basis);
        up2coef := matrix_mget_up2coef(matrixj);
        symdata := hashtable_htget_hashdata(symbol_ht);

        for i := 1:basis_bget_nlead(basis) do <<
            nrows := matrix_mget_nrows(matrixj) + 1;
            matrix_mset_nrows(matrixj, nrows);

            putv(uprows, nrows, hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, 0, etmp, getv(gens, getv(nonred, i))));

            putv(up2coef, nrows, getv(nonred, i));
            % set lead index as 1
            hashtable_hvset_idx(getv(symdata, getv(getv(uprows, nrows), 1)), 1)
        >>;

        % needed for correct counting in symbol
        matrix_mset_ncols(matrixj, matrix_mget_nrows(matrixj));
        matrix_mset_nup(matrixj, matrix_mget_nrows(matrixj));

        f4_symbolic_preprocessing(basis, matrixj, ht, symbol_ht);
        for i := hashtable_htget_offset(symbol_ht):hashtable_htget_load(symbol_ht) do
            hashtable_hvset_idx(getv(symdata, i), 1);

        matrix_convert_hashes_to_columns(matrixj, symbol_ht);
        matrix_mset_ncols(matrixj, matrix_mget_nleft(matrixj) + matrix_mget_nright(matrixj));

        sorting_sort_matrix_rows_decreasing(matrixj);

        matrix_interreduce_matrix_rows(ring, matrixj, basis);

        matrix_convert_matrix_rows_to_basis_elements(matrixj, basis, ht, symbol_ht);
        % no longer need in two hashtables

        basis_bset_ntotal(basis, matrix_mget_npivots(matrixj) + basis_bget_ndone(basis));
        basis_bset_ndone(basis, matrix_mget_npivots(matrixj));

        % we may have added some multiples of reduced basis polynomials
        % from the matrixj, so we get rid of them.
        k := 0;
        i := 1;
        gens := basis_bget_gens(basis);
        nonred := basis_bget_nonred(basis);
        ntotal := basis_bget_ntotal(basis);
        lead := basis_bget_lead(basis);
        htdata := hashtable_htget_hashdata(ht);
    Letsgo:
        while i <= basis_bget_ndone(basis) do <<
            for j := 1:k do <<
                if hashtable_is_monom_divisible(
                        getv(getv(gens, ntotal - i + 1), 1),
                        getv(getv(gens, getv(nonred, j)), 1),
                        ht) then <<

                    i := i + 1;
                    go to Letsgo
                >>
            >>;
            k := k + 1;
            putv(nonred, k, ntotal - i + 1);
            putv(lead, k, hashtable_hvget_divmask(getv(htdata, getv(getv(gens, getv(nonred, k)), 1))));
            i := i + 1
        >>;
        basis_bset_nlead(basis, k)

    end;

%--------------------------------------------------------------------------------------------------

asserted procedure f4_find_multiplied_reducer(basis: Basis, matrixj: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable,
                                vidx: Integer);
    begin scalar symexps, symdata, e, etmp, divmask, blen, leaddiv, nlead, gens, 
                    nonred, htexps, explen, uprows, up2coef, i, rpoly, rexp,
                    nup, h;
        symexps := hashtable_htget_exponents(symbol_ht);
        symdata := hashtable_htget_hashdata(symbol_ht);

        e := getv(symexps, vidx);
        etmp := getv(symexps, 1);
        divmask := hashtable_hvget_divmask(getv(symdata, vidx));

        blen := basis_bget_ndone(basis);
        leaddiv := basis_bget_lead(basis);
        nlead := basis_bget_nlead(basis);
        
        gens := basis_bget_gens(basis);
        htexps := hashtable_htget_exponents(ht);
        nonred := basis_bget_nonred(basis);
        explen := hashtable_htget_explen(ht);

        uprows := matrix_mget_uprows(matrixj);
        up2coef := matrix_mget_up2coef(matrixj);

        % searching for a poly from whose leading monom
        % divides the given exponent e
        i := 1;
    Letsgo:
        while i <= nlead and not (land(getv(leaddiv, i), lnot(divmask)) = 0) do
            i := i + 1;

        % here found polynomial from basis with leading monom
        % dividing symbol_ht.exponents[vidx]
        if i <= nlead then <<
            rpoly := getv(gens, getv(nonred, i));
            rexp := getv(htexps, getv(rpoly, 1));

            for j := 1:explen do <<
                % if it actually does not divide and divmask lies
                if getv(e, j) < getv(rexp, j) then <<
                    i := i + 1;
                    go to Letsgo
                >>;
                putv(etmp, j, getv(e, j) - getv(rexp, j))
            >>;
            % now etmp = e // rexp in terms of monomias,
            % hash is linear
            h := hashtable_hvget_hash(getv(symdata, vidx)) - hashtable_hvget_hash(getv(symdata, getv(rpoly, 1)));
            
            nup := matrix_mget_nup(matrixj);
            putv(uprows, nup + 1, hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, h, etmp, rpoly));
            putv(up2coef, nup + 1, getv(nonred, i));

            % upsize matrixj
            hashtable_hvset_idx(getv(symdata, vidx), 2);
            matrix_mset_nup(matrixj, nup + 1);
            i := i + 1
        >>

    end;

%--------------------------------------------------------------------------------------------------

% Julia: comment in Julia is not correct, we don't copy it here
asserted procedure f4_symbolic_preprocessing(basis: Basis, matrixj: MacaulayMatrix, 
                                            ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar symbol_load, nrr, onrr, i, symdata;

        % check matrixj sizes (we want to omit this I guess)

        symbol_load := hashtable_htget_load(symbol_ht);

        nrr := matrix_mget_ncols(matrixj);
        onrr := matrix_mget_ncols(matrixj);

        while matrix_mget_size(matrixj) <= nrr + symbol_load do <<
            % Julia: Ask Thomas
            matrix_mset_size(matrixj, matrix_mget_size(matrixj)*2);
            matrix_mset_uprows(matrixj, dv_resize(matrix_mget_uprows(matrixj), matrix_mget_size(matrixj)));
            matrix_mset_up2coef(matrixj, dv_resize(matrix_mget_up2coef(matrixj), matrix_mget_size(matrixj)))
        >>;

        % for each lcm present in symbolic_ht set on select stage
        i := hashtable_htget_offset(symbol_ht);
        % First round, we add multiplied polynomials which divide  =#
        % a monomial exponent from selected spairs
        while i <= symbol_load do <<
            symdata := hashtable_htget_hashdata(symbol_ht);
            % not a reducer already
            if hashtable_hvget_idx(getv(symdata, i)) = 0 then <<
                hashtable_hvset_idx(getv(symdata, i), 1);
                matrix_mset_ncols(matrixj, matrix_mget_ncols(matrixj) + 1);
                f4_find_multiplied_reducer(basis, matrixj, ht, symbol_ht, i)
            >>;
            i := i + 1
        >>;

        % Second round, we add multiplied polynomials which divide  =#
        % lcm added on previous for loop
        while i <= hashtable_htget_load(symbol_ht) do <<
            if matrix_mget_size(matrixj) = matrix_mget_nup(matrixj) then <<
                matrix_mset_size(matrixj, matrix_mget_size(matrixj) * 2);
                matrix_mset_uprows(matrixj, dv_resize(matrix_mget_uprows(matrixj), matrix_mget_size(matrixj)));
                matrix_mset_up2coef(matrixj, dv_resize(matrix_mget_up2coef(matrixj), matrix_mget_size(matrixj)))
            >>;

            symdata := hashtable_htget_hashdata(symbol_ht);
            hashtable_hvset_idx(getv(symdata, i), 1);
            matrix_mset_ncols(matrixj, matrix_mget_ncols(matrixj) + 1);
            f4_find_multiplied_reducer(basis, matrixj, ht, symbol_ht, i);
            i := i + 1
        >>;

        % shrink matrixj sizes, set constants
        matrix_mset_uprows(matrixj, dv_resize(matrix_mget_uprows(matrixj), matrix_mget_nup(matrixj)));

        matrix_mset_nrows(matrixj, matrix_mget_nrows(matrixj) + matrix_mget_nup(matrixj) - onrr);
        matrix_mset_nlow(matrixj, matrix_mget_nrows(matrixj) - matrix_mget_nup(matrixj));
        matrix_mset_size(matrixj, matrix_mget_nrows(matrixj))

    end;

%--------------------------------------------------------------------------------------------------

asserted procedure f4_select_normal(pairset: Pairset, basis: Basis, matrixj: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar ps, min_idx, npairs, uprows, up2coef, lowrows,
                    gens, htexps, htdata, symdata, etmp, i,
                    loadj, lcmj, j, bgens, prev, poly, vidx, eidx, elcm,
                    htmp, min_deg, explen, low2coef;
        
        sorting_sort_pairset_by_degree(pairset, 1, basis_psget_load(pairset) - 1);

        ps := basis_psget_pairs(pairset);
        min_deg := hashtable_hvget_deg(getv(ps, 1));
        min_idx := 0;

        while min_idx < basis_psget_load(pairset) and hashtable_hvget_deg(getv(ps, min_idx+1)) = min_deg do
            min_idx := min_idx + 1;

        % number of selected pairs
        npairs := min_idx;

        if f4_debug() then
            prin2t {"Selected", npairs, "pairs"};

        sorting_sort_pairset_by_lcm(pairset, npairs, ht);

        matrix_reinitialize_matrix(matrixj, npairs);

        uprows := matrix_mget_uprows(matrixj);
        up2coef := matrix_mget_up2coef(matrixj);
        low2coef := matrix_mget_low2coef(matrixj);
        lowrows := matrix_mget_lowrows(matrixj);

        % polynomials from pairs in order (p11, p12)(p21, p21)
        % (future rows of the matrixj)
        gens := dv_undef(2 * npairs);
        
        htexps := hashtable_htget_exponents(ht);
        htdata := hashtable_htget_hashdata(ht);

        explen := hashtable_htget_explen(ht);

        symdata := hashtable_htget_hashdata(symbol_ht);

        etmp := getv(htexps, 1);
        i := 1;
        while i <= npairs do <<
            matrix_mset_ncols(matrixj, matrix_mget_ncols(matrixj) + 1);
            loadj := 1;
            lcmj := basis_spget_lcmj(getv(ps, i));
            j := i;

            % we collect all generators with same lcm into gens
            while j <= npairs and basis_spget_lcmj(getv(ps, j)) = lcmj do <<
                putv(gens, loadj, basis_spget_poly1(getv(ps, j)));
                loadj := loadj + 1;
                putv(gens, loadj, basis_spget_poly2(getv(ps, j)));
                loadj := loadj + 1;
                j := j + 1
            >>;
            loadj := loadj - 1;

            % sort by number in the basis (by=identity)
            sorting_sort_generators_by_position(gens, loadj);

            % now we collect reducers, and reduced

            % Julia: basis.gens
            bgens := basis_bget_gens(basis);

            % first generator index in groebner basis
            prev := getv(gens, 1);
            % first generator in hash table
            poly := getv(bgens, prev);
            % first generator lead monomial index in hash data
            vidx := getv(poly, 1);

            % first generator exponent
            eidx := getv(htexps, vidx);
            % exponent of lcm corresponding to first generator
            elcm := getv(htexps, lcmj);
            for u := 1:explen do
                putv(etmp, u, getv(elcm, u) - getv(eidx, u));
            % now etmp contents complement to eidx in elcm

            % hash of complement
            htmp := hashtable_hvget_hash(getv(htdata, lcmj)) - hashtable_hvget_hash(getv(htdata, vidx));

            % add row as a reducer
            matrix_mset_nup(matrixj, matrix_mget_nup(matrixj) + 1);
            putv(uprows, matrix_mget_nup(matrixj), hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, htmp, etmp, poly));
            % map upper row to index in basis
            putv(up2coef, matrix_mget_nup(matrixj), prev);

            %  mark lcm column as reducer in symbolic hashtable
            hashtable_hvset_idx(getv(symdata, getv(getv(uprows, matrix_mget_nup(matrixj)), 1)), 2);
            % increase number of rows set
            matrix_mset_nrows(matrixj, matrix_mget_nrows(matrixj) + 1);

            % over all polys with same lcm,
            % add them to the lower part of matrixj
            for k := 1:loadj do <<
                % duplicate generator,
                % we can do so as long as generators are sorted
                if not (getv(gens, k) = prev) then <<
                    % if the table was reallocated
                    elcm := getv(htexps, lcmj);

                    prev := getv(gens, k);
                    poly := getv(bgens, prev);
                    vidx := getv(poly, 1);
                    %  leading monom idx
                    eidx := getv(htexps, vidx);
                    for u := 1:explen do
                        putv(etmp, u, getv(elcm, u) - getv(eidx, u));
                    
                    htmp := hashtable_hvget_hash(getv(htdata, lcmj)) - hashtable_hvget_hash(getv(htdata, vidx));

                    %  add row to be reduced
                    matrix_mset_nlow(matrixj, matrix_mget_nlow(matrixj) + 1);
                    putv(lowrows, matrix_mget_nlow(matrixj), hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, htmp, etmp, poly));
                    % map lower row to index in basis
                    putv(low2coef, matrix_mget_nlow(matrixj), prev);

                    hashtable_hvset_idx(getv(getv(lowrows, matrix_mget_nlow(matrixj)), 1), 2);

                    matrix_mset_nrows(matrixj, matrix_mget_nrows(matrixj) + 1)
                >>
            >>;

            i := j
        >>;

        matrix_mset_lowrows(matrixj, dv_resize(matrix_mget_lowrows(matrixj), matrix_mget_nrows(matrixj) - matrix_mget_ncols(matrixj)));

        % remove selected parirs from pairset
        for i := 1:basis_psget_load(pairset) - npairs do
            putv(ps, i, getv(ps, i + npairs));
        basis_psset_load(pairset, basis_psget_load(pairset) - npairs)
    end;

%--------------------------------------------------------------------------------------------------

% Computes the Groebner basis inplace in `basis`
%
% Input ivariants:
%   - ring is set, and ring fields ~ ht fields
%   - divmasks in ht are set
%   - basis is filled so that
%   - basis.ntotal = actual number of elements
%   - basis.ndone  = 0
%   - basis.nlead  = 0
%   - basis contains no zero polynomials (!)
% Output invariants:
%   - basis.ndone == basis.ntotal == basis.nlead
%   - basis.gens and basis.coeffs of size basis.ndone
%   - basis elements are sorted increasingly wrt ordering on lead elements
%   - divmasks in basis are filled and coincide to divmasks in hashtable
%
% Julia: in Reduce, there is no `tracer`, `linalg`, and `rng` arguments
asserted procedure f4_f4(ring: PolyRing, basis: Basis, ht: MonomialHashtable, reduced: Bool);
    begin scalar pairset, matrixj, update_ht, symbol_ht, plcm, d;

        ASSERT(io_prget_ord(ring) = hashtable_htget_ord(ht));
        ASSERT(io_prget_nvars(ring) = hashtable_htget_nvars(ht));
        ASSERT(io_prget_explen(ring) = hashtable_htget_explen(ht));

        ASSERT(basis_bget_ndone(basis) = 0);

        if f4_debug() then <<
            prin2t "f4: Input";
            prin2t {"Ring: ", ring};
            prin2t {"Basis: ", basis};
            prin2t {"HT: ", ht}
        >>;

        % Julia: in Reduce, Pairset initialization is moved here
        pairset := basis_initialize_pairset();

        % matrixj storing coefficients in rows
        % wrt columns representing the current monomial basis
        matrixj := matrix_initialize_matrix(ring);

        % initialize hash tables for update and symbolic preprocessing steps
        update_ht := hashtable_initialize_secondary_hash_table(ht);
        symbol_ht := hashtable_initialize_secondary_hash_table(ht);
        
        if f4_debug() then <<
            prin2t "f4: After initialize";
            prin2t {"Pairset: ", pairset};
            prin2t {"Matrix: ", matrixj};
            prin2t {"UpdapteHT: ", update_ht}
        >>;

        plcm := dv_undef(0);
        plcm := basis_update(pairset, basis, ht, update_ht, plcm);
        
        if f4_debug() then <<
            prin2t "After first update";
            prin2t {"Pairset: ", pairset};
            prin2t {"Basis: ", basis};
            prin2t {"HT: ", ht}
        >>;

        % while there are pairs to be reduced
        d := 0;
        while not (basis_psget_load(pairset) = 0) do <<
            d := d + 1;
            if f4_debug() then <<
                prin2t {"f4: iteration: ", d};
                prin2t {"f4: available pairs:", basis_psget_load(pairset)}
            >>;

            % selects pairs for reduction from pairset following normal strategy
            % (minimal lcm degrees are selected),
            % and puts these into the matrixj rows
            f4_select_normal(pairset, basis, matrixj, ht, symbol_ht);

            if f4_debug() then <<
                prin2t {"f4: after select: ", pairset};
                prin2t {"f4: matrix:", matrixj}
            >>;

            f4_symbolic_preprocessing(basis, matrixj, ht, symbol_ht);

            if f4_debug() then <<
                prin2t {"f4: after symbolic: ", matrixj}
            >>;

            % reduces polys and obtains new potential basis elements
            f4_reduction(ring, basis, matrixj, ht, symbol_ht);

            if f4_debug() then <<
                prin2t {"f4: after reduction: ", matrixj}
            >>;

            % update the current basis with polynomials produced from reduction,
            % does not copy,
            % checks for redundancy
            plcm := basis_update(pairset, basis, ht, update_ht, plcm);

            matrixj := matrix_initialize_matrix(ring);
            symbol_ht := hash_initialize_secondary_hashtable(ht);
        >>;

        % remove redundant elements
        basis_filter_redundant(basis);
        
        if f4_debug() then <<
            prin2t {"f4: after filter redundant: ", basis}
        >>;

        if reduced then
            f4_reducegb(ring, basis, matrixj, ht, symbol_ht);
        
        basis_standardize_basis(ring, basis, ht, hashtable_htget_ord(ht));

        if f4_debug() then <<
            prin2t {"f4: after standardize: ", basis}
        >>;

    end;

endmodule; % end of module f4;

end; % end of file