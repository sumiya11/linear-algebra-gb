module f4f4;
% The F4 implementation. The heart of this package.
% This file corresponds to file f4/f4.jl in Groebner.jl

revision('f4f4, "$Id$");

copyright('f4f4, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%    * Redistributions of source code must retain the relevant
%      copyright notice, this list of conditions and the following
%      disclaimer.
%    * Redistributions in binary form must reproduce the above
%      copyright notice, this list of conditions and the following
%      disclaimer in the documentation and/or other materials provided
%      with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

% Sort rows of the given matrix and call matrix reduction
asserted procedure f4_reduction(ring: PolyRing, basis: Basis, matrixj: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    <<
        matrix_convert_hashes_to_columns(matrixj, symbol_ht);

        sorting_sort_matrix_rows_decreasing(matrixj); % for pivots,  AB part
        sorting_sort_matrix_rows_increasing(matrixj); % for reduced, CD part

        matrix_linear_algebra(ring, matrixj, basis);

        matrix_convert_matrix_rows_to_basis_elements(matrixj, basis, ht, symbol_ht)
    >>;

% Initializes Basis and MonomialHashtable structures,
% fills input data from exponents and coeffs
%
% MonomialHashtable initial size is set to tablesize
asserted procedure f4_initialize_structures(ring: PolyRing, exponents: Vector, 
                                            coeffs: Vector, tablesize: Integer);
    begin scalar basis, basis_ht;
        % basis for storing basis elements,
        % hashtable for hashing monomials occuring in the basis
        basis := basis_initialize_basis(ring, dv_length(exponents));
        basis_ht := hashtable_initialize_basis_hash_table(ring, tablesize);

        % filling the basis and hashtable with the given inputs
        basis_fill_data(basis, basis_ht, exponents, coeffs);

        hashtable_fill_divmask(basis_ht);

        sorting_sort_gens_by_lead_increasing(basis, basis_ht);

        basis_normalize_basis(ring, basis);

        return basis . basis_ht
    end;

% Initializes Basis and MonomialHashtable structures,
% fills input data from exponents and coeffs
%
% MonomialHashtable initial size is set to tablesize
asserted procedure f4_initialize_structures_no_normalize(ring: PolyRing, exponents: Vector, 
                                            coeffs_qq: Vector, coeffs_ff: Vector, tablesize: Integer);
    begin scalar basis, basis_ht;
        % basis for storing basis elements,
        % hashtable for hashing monomials occuring in the basis
        basis := basis_initialize_basis(ring, dv_length(exponents));
        basis_ht := hashtable_initialize_basis_hash_table(ring, tablesize);
        
        % filling the basis and hashtable with the given inputs
        basis_fill_data(basis, basis_ht, exponents, coeffs_ff);

        % every monomial in hashtable is associated with its divmask
        % to perform divisions faster. Filling those
        hashtable_fill_divmask(basis_ht);
        
        % sort input, smaller leading terms first
        sorting_sort_gens_by_lead_increasing2(basis, basis_ht, coeffs_qq);

        return basis . basis_ht
    end;

% Initializes Basis and MonomialHashtable structures,
% fills input data from exponents and coeffs
%
% MonomialHashtable initial size is set to tablesize
asserted procedure f4_initialize_structures_ff(ring: PolyRing, exponents: Vector, 
                                            coeffs: Vector, tablesize: Integer);
    begin scalar coeffs_ff;
        coeffs_ff := dv_undef(dv_length(coeffs));
        for i := 1:dv_length(coeffs_ff) do
            putv(coeffs_ff, i, dv_undef(dv_length(coeffs[i])));
        return f4_initialize_structures_no_normalize(ring, exponents, coeffs, coeffs_ff, tablesize)
    end;

% "Interreduce" the given basis inplace
asserted procedure f4_reducegb(ring: PolyRing, basis: Basis, matrixj: MacaulayMatrix,
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar exponents, etmp, uprows, nonred, gens, up2coef, symdata, nrows,
                    k, lead, htdata, ntotal, i;
        
        exponents := hashtable_htget_exponents(ht);
        etmp := exponents[1];
        
        for i := 1:dv_length(etmp) do
            etmp[i] := 0;
        %  etmp is now set to zero, and has zero hash

        matrix_reinitialize_matrix(matrixj, basis_bget_nlead(basis));
        uprows := matrix_mget_uprows(matrixj);

        nonred := basis_bget_nonred(basis);
        gens := basis_bget_gens(basis);
        up2coef := matrix_mget_up2coef(matrixj);
        symdata := hashtable_htget_hashdata(symbol_ht);

        % add all non redundant elements from basis
        % as matrix upper rows
        for i := 1:basis_bget_nlead(basis) do <<
            nrows := matrix_mget_nrows(matrixj) + 1;
            matrix_mset_nrows(matrixj, nrows);

            putv(matrix_mget_uprows(matrixj), nrows, hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, 0, etmp, getv(gens, getv(nonred, i))));

            putv(matrix_mget_up2coef(matrixj), nrows, getv(nonred, i));
            % set lead index as 1
            hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), getv(getv(matrix_mget_uprows(matrixj), nrows), 1)), 1)
        >>;

        % needed for correct counting in symbol
        matrix_mset_ncols(matrixj, matrix_mget_nrows(matrixj));
        matrix_mset_nup(matrixj, matrix_mget_nrows(matrixj));

        f4_symbolic_preprocessing(basis, matrixj, ht, symbol_ht);
        for i := hashtable_htget_offset(symbol_ht):hashtable_htget_load(symbol_ht) do
            hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), i), 1);

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
            putv(lead, k, hashtable_hvget_divmask(getv(hashtable_htget_hashdata(ht), getv(getv(gens, getv(nonred, k)), 1))));
            i := i + 1
        >>;
        basis_bset_nlead(basis, k)

    end;

asserted procedure f4_find_multiplied_reducer(basis: Basis, matrixj: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable,
                                vidx: Integer);
    begin scalar symdata, e, etmp, divmask, blen, leaddiv, gens, 
                    nonred, explen, uprows, up2coef, i, rpoly, rexp,
                    nup, h;

        e := getv(hashtable_htget_exponents(symbol_ht), vidx);
        etmp := getv(hashtable_htget_exponents(ht), 1);
        divmask := hashtable_hvget_divmask(getv(hashtable_htget_hashdata(symbol_ht), vidx));

        blen := basis_bget_ndone(basis);
        leaddiv := basis_bget_lead(basis);
        
        gens := basis_bget_gens(basis);
        nonred := basis_bget_nonred(basis);
        explen := hashtable_htget_explen(ht);

        uprows := matrix_mget_uprows(matrixj);
        up2coef := matrix_mget_up2coef(matrixj);

        % searching for a poly from whose leading monom
        % divides the given exponent e
        i := 1;
    Letsgo:
        while i <= basis_bget_nlead(basis) and not (land(getv(leaddiv, i), lnot(divmask)) = 0) do
            i := i + 1;
        
        % here found polynomial from basis with leading monom
        % dividing symbol_ht.exponents[vidx]
        if i <= basis_bget_nlead(basis) then <<
            rpoly := getv(gens, getv(nonred, i));
            rexp := getv(hashtable_htget_exponents(ht), getv(rpoly, 1));

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
            h := hashtable_hvget_hash(getv(hashtable_htget_hashdata(symbol_ht), vidx)) - hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), getv(rpoly, 1)));
            
            nup := matrix_mget_nup(matrixj);
            putv(matrix_mget_uprows(matrixj), nup + 1, hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, h, etmp, rpoly));
            putv(matrix_mget_up2coef(matrixj), nup + 1, getv(nonred, i));
            
            symdata := hashtable_htget_hashdata(symbol_ht);

            % upsize matrixj
            hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), vidx), 2);
            matrix_mset_nup(matrixj, nup + 1);
            i := i + 1
        >>

    end;

% Julia: comment in Julia is not correct, we don't copy it here
asserted procedure f4_symbolic_preprocessing(basis: Basis, matrixj: MacaulayMatrix, 
                                            ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar symbol_load, nrr, onrr, i;

        % check matrixj sizes (we want to omit this I guess)

        symbol_load := hashtable_htget_load(symbol_ht);

        nrr := matrix_mget_ncols(matrixj);
        onrr := matrix_mget_ncols(matrixj);

        while matrix_mget_size(matrixj) <= nrr + symbol_load do <<
            matrix_mset_size(matrixj, matrix_mget_size(matrixj)*2);
            matrix_mset_uprows(matrixj, dv_resize(matrix_mget_uprows(matrixj), matrix_mget_size(matrixj)));
            matrix_mset_up2coef(matrixj, dv_resize(matrix_mget_up2coef(matrixj), matrix_mget_size(matrixj)))
        >>;

        % for each lcm present in symbolic_ht set on select stage
        i := hashtable_htget_offset(symbol_ht);
        % First round, we add multiplied polynomials which divide  =#
        % a monomial exponent from selected spairs
        while i <= symbol_load do <<
            % not a reducer already
            if hashtable_hvget_idx(getv(hashtable_htget_hashdata(symbol_ht), i)) = 0 then <<
                hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), i), 1);
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
            
            hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), i), 1);
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

% for history
% smacro procedure symdata(symbol_ht);
%     hashtable_htget_hashdata(symbol_ht);

asserted procedure f4_select_normal(pairset: Pairset, basis: Basis, matrixj: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar ps, min_idx, npairs,
                    gens, htexps, htdata, etmp, i, symdata,
                    loadj, lcmj, j, bgens, prev, poly, vidx, eidx, elcm,
                    htmp, min_deg, explen;
        
        sorting_sort_pairset_by_degree(pairset, 1, basis_psget_load(pairset) - 1);

        ps := basis_psget_pairs(pairset);
        min_deg := hashtable_hvget_deg(getv(ps, 1));
        min_idx := 0;

        while min_idx < basis_psget_load(pairset) and hashtable_hvget_deg(getv(ps, min_idx+1)) = min_deg do
            min_idx := min_idx + 1;

        % number of selected pairs
        npairs := min_idx;
        sorting_sort_pairset_by_lcm(pairset, npairs, ht);

        matrix_reinitialize_matrix(matrixj, npairs);

        % polynomials from pairs in order (p11, p12)(p21, p21)
        % (future rows of the matrixj)
        gens := dv_undef(2 * npairs);
        
        htexps := hashtable_htget_exponents(ht);
        htdata := hashtable_htget_hashdata(ht);

        explen := hashtable_htget_explen(ht);

        symdata := hashtable_htget_hashdata(symbol_ht);

        etmp := getv(hashtable_htget_exponents(ht), 1);
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
            eidx := getv(hashtable_htget_exponents(ht), vidx);
            % exponent of lcm corresponding to first generator
            elcm := getv(hashtable_htget_exponents(ht), lcmj);
            for u := 1:explen do
                putv(etmp, u, getv(elcm, u) - getv(eidx, u));
            % now etmp contents complement to eidx in elcm

            % hash of complement
            htmp := hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), lcmj)) - hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), vidx));

            % add row as a reducer
            matrix_mset_nup(matrixj, matrix_mget_nup(matrixj) + 1);

            putv(matrix_mget_uprows(matrixj), matrix_mget_nup(matrixj), hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, htmp, etmp, poly));

            symdata := hashtable_htget_hashdata(symbol_ht);

            % map upper row to index in basis
            putv(matrix_mget_up2coef(matrixj), matrix_mget_nup(matrixj), prev);

            %  mark lcm column as reducer in symbolic hashtable
            hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), getv(getv(matrix_mget_uprows(matrixj), matrix_mget_nup(matrixj)), 1)), 2);
            
            % increase number of rows set
            matrix_mset_nrows(matrixj, matrix_mget_nrows(matrixj) + 1);

            % over all polys with same lcm,
            % add them to the lower part of matrixj
            for k := 1:loadj do <<
                % duplicate generator,
                % we can do so as long as generators are sorted
                if getv(gens, k) neq prev then <<
                    % if the table was reallocated
                    elcm := getv(hashtable_htget_exponents(ht), lcmj);

                    prev := getv(gens, k);
                    poly := getv(bgens, prev);
                    vidx := getv(poly, 1);
                    %  leading monom idx
                    eidx := getv(hashtable_htget_exponents(ht), vidx);
                    for u := 1:explen do
                        putv(etmp, u, getv(elcm, u) - getv(eidx, u));
                    
                    htmp := hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), lcmj)) - hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), vidx));

                    %  add row to be reduced
                    matrix_mset_nlow(matrixj, matrix_mget_nlow(matrixj) + 1);
                    putv(matrix_mget_lowrows(matrixj), matrix_mget_nlow(matrixj), hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, htmp, etmp, poly));
                    % map lower row to index in basis
                    putv(matrix_mget_low2coef(matrixj), matrix_mget_nlow(matrixj), prev);

                    hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), getv(getv(matrix_mget_lowrows(matrixj), matrix_mget_nlow(matrixj)), 1)), 2);

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

asserted procedure f4_select_isgroebner(pairset: Pairset, basis: Basis, matrixj: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar ps, npairs,
                    gens, htexps, htdata, etmp, i, symdata,
                    loadj, lcmj, j, bgens, prev, poly, vidx, eidx, elcm,
                    htmp, explen;
        
        % sorting_sort_pairset_by_degree(pairset, 1, basis_psget_load(pairset) - 1);

        ps := basis_psget_pairs(pairset);
        npairs := basis_psget_load(pairset);

        sorting_sort_pairset_by_lcm(pairset, npairs, ht);

        matrix_reinitialize_matrix(matrixj, npairs);

        % polynomials from pairs in order (p11, p12)(p21, p21)
        % (future rows of the matrixj)
        gens := dv_undef(2 * npairs);
        
        htexps := hashtable_htget_exponents(ht);
        htdata := hashtable_htget_hashdata(ht);

        explen := hashtable_htget_explen(ht);

        symdata := hashtable_htget_hashdata(symbol_ht);

        etmp := getv(hashtable_htget_exponents(ht), 1);
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
            eidx := getv(hashtable_htget_exponents(ht), vidx);
            % exponent of lcm corresponding to first generator
            elcm := getv(hashtable_htget_exponents(ht), lcmj);
            for u := 1:explen do
                putv(etmp, u, getv(elcm, u) - getv(eidx, u));
            % now etmp contents complement to eidx in elcm

            % hash of complement
            htmp := hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), lcmj)) - hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), vidx));

            % add row as a reducer
            matrix_mset_nup(matrixj, matrix_mget_nup(matrixj) + 1);

            putv(matrix_mget_uprows(matrixj), matrix_mget_nup(matrixj), hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, htmp, etmp, poly));

            symdata := hashtable_htget_hashdata(symbol_ht);

            % map upper row to index in basis
            putv(matrix_mget_up2coef(matrixj), matrix_mget_nup(matrixj), prev);

            %  mark lcm column as reducer in symbolic hashtable
            hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), getv(getv(matrix_mget_uprows(matrixj), matrix_mget_nup(matrixj)), 1)), 2);
            
            % increase number of rows set
            matrix_mset_nrows(matrixj, matrix_mget_nrows(matrixj) + 1);

            % over all polys with same lcm,
            % add them to the lower part of matrixj
            for k := 1:loadj do <<
                % duplicate generator,
                % we can do so as long as generators are sorted
                if getv(gens, k) neq prev then <<
                    % if the table was reallocated
                    elcm := getv(hashtable_htget_exponents(ht), lcmj);

                    prev := getv(gens, k);
                    poly := getv(bgens, prev);
                    vidx := getv(poly, 1);
                    %  leading monom idx
                    eidx := getv(hashtable_htget_exponents(ht), vidx);
                    for u := 1:explen do
                        putv(etmp, u, getv(elcm, u) - getv(eidx, u));
                    
                    htmp := hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), lcmj)) - hashtable_hvget_hash(getv(hashtable_htget_hashdata(ht), vidx));

                    %  add row to be reduced
                    matrix_mset_nlow(matrixj, matrix_mget_nlow(matrixj) + 1);
                    putv(matrix_mget_lowrows(matrixj), matrix_mget_nlow(matrixj), hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, htmp, etmp, poly));
                    % map lower row to index in basis
                    putv(matrix_mget_low2coef(matrixj), matrix_mget_nlow(matrixj), prev);

                    hashtable_hvset_idx(getv(hashtable_htget_hashdata(symbol_ht), getv(getv(matrix_mget_lowrows(matrixj), matrix_mget_nlow(matrixj)), 1)), 2);

                    matrix_mset_nrows(matrixj, matrix_mget_nrows(matrixj) + 1)
                >>
            >>;

            i := j
        >>;

        matrix_mset_lowrows(matrixj, dv_resize(matrix_mget_lowrows(matrixj), matrix_mget_nrows(matrixj) - matrix_mget_ncols(matrixj)))
    end;

asserted procedure f4_select_tobereduced(basis: Basis, tobereduced: Basis, matrixj: MacaulayMatrix, symbol_ht: MonomialHashtable, ht: MonomialHashtable);
    begin scalar etmp, gen, h;
        
        % prepare to load all elems from tobereduced
        % into low rows
        matrix_reinitialize_matrix(matrixj, max(basis_bget_ntotal(basis), basis_bget_ntotal(tobereduced)));
        matrix_mset_lowrows(matrixj, dv_resize(matrix_mget_lowrows(matrixj), basis_bget_ntotal(tobereduced)));

        etmp := dv_zeros(hashtable_htget_explen(ht));

        for i := 1:basis_bget_ntotal(tobereduced) do <<
            matrix_mset_nrows(matrixj, matrix_mget_nrows(matrixj) + 1);
            gen := getv(basis_bget_gens(tobereduced), i);
            h := 0;
            putv(matrix_mget_lowrows(matrixj), matrix_mget_nrows(matrixj), hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, h, etmp, gen));
            putv(matrix_mget_low2coef(matrixj), matrix_mget_nrows(matrixj), i)
        >>;

        basis_bset_nlead(basis, basis_bget_ntotal(basis));
        basis_bset_ndone(basis, basis_bget_ntotal(basis));
        basis_bset_isred(basis, dv_zeros(dv_length(basis_bget_isred(basis))));
        for i := 1:basis_bget_nlead(basis) do <<
            putv(basis_bget_nonred(basis), i, i);
            putv(basis_bget_lead(basis), i, hashtable_hvget_divmask(getv(hashtable_htget_hashdata(ht), getv(getv(basis_bget_gens(basis), i), 1))))
        >>
    end;

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
asserted procedure f4_f4(ring: PolyRing, basis: Basis, ht: MonomialHashtable, reduced: Boolean);
    begin scalar pairset, matrixj, update_ht, symbol_ht, plcm, d;

        ASSERT(io_prget_ord(ring) = hashtable_htget_ord(ht));
        ASSERT(io_prget_nvars(ring) = hashtable_htget_nvars(ht));
        ASSERT(io_prget_explen(ring) = hashtable_htget_explen(ht));

        ASSERT(basis_bget_ndone(basis) = 0);

        % in Reduce, Pairset initialization is moved here
        pairset := basis_initialize_pairset();

        % matrixj storing coefficients in rows
        % wrt columns representing the current monomial basis
        matrixj := matrix_initialize_matrix(ring);

        % initialize hash tables for update and symbolic preprocessing steps
        update_ht := hashtable_initialize_secondary_hash_table(ht);
        symbol_ht := hashtable_initialize_secondary_hash_table(ht);

        plcm := dv_undef(0);
        plcm := basis_update(pairset, basis, ht, update_ht, plcm);

        % while there are pairs to be reduced
        d := 0;
        while basis_psget_load(pairset) neq 0 do <<
            d := d + 1;

            % selects pairs for reduction from pairset following normal strategy
            % (minimal lcm degrees are selected),
            % and puts these into the matrixj rows
            f4_select_normal(pairset, basis, matrixj, ht, symbol_ht);

            f4_symbolic_preprocessing(basis, matrixj, ht, symbol_ht);

            % reduces polys and obtains new potential basis elements
            f4_reduction(ring, basis, matrixj, ht, symbol_ht);

            % update the current basis with polynomials produced from reduction,
            % does not copy,
            % checks for redundancy
            plcm := basis_update(pairset, basis, ht, update_ht, plcm);

            matrixj := matrix_initialize_matrix(ring);
            symbol_ht := hashtable_initialize_secondary_hash_table(ht)
        >>;

        % remove redundant elements
        basis_filter_redundant(basis);

        if reduced then
            f4_reducegb(ring, basis, matrixj, ht, symbol_ht);
        
        basis_standardize_basis(ring, basis, ht, hashtable_htget_ord(ht));

    end;

endmodule; % end of module f4;

end; % end of file