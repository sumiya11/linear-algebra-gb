module f4;

%--------------------------------------------------------------------------------------------------

asserted procedure f4_reduction(ring: PolyRing, basis: Basis, matrix: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    <<
        matrix_convert_hashes_to_columns(matrix, symbol_ht);

        sorting_sort_matrix_rows_decreasing(matrix);
        sorting_sort_matrix_rows_increasing(matrix);

        matrix_linear_algebra(ring, matrix, basis);

        matrix_convert_matrix_rows_to_basis_elements(matrix, basis, ht, symbol_ht)
    >>;

%--------------------------------------------------------------------------------------------------

asserted procedure f4_initialize_structures(ring: PolyRing, exponents: Vector, 
                                            coeffs: Vector, tablesize: Integer);
    begin scalar basis, basis_ht;
        % basis for storing basis elements,
        % hashtable for hashing monomials occuring in the basis
        basis := basis_initialize_basis(ring);
        basis_ht := hashtable_initialize_basis_hashtable(ring);
        
        % filling the basis and hashtable with the given inputs
        basis_fill_Data(basis, basis_ht, exponents, coeffs);
        
        hashtable_fill_divmask(basis_ht);

        sorting_sort_gens_by_lead_increasing(basis, basis_ht);

        basis_normalize_basis(ring, basis);

        return basis . ht
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure f4_reducegb(ring: PolyRing, basis: Basis, matrix: MacaulayMatrix,
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar etmp;
        
        exponents := hashtable_htget_exponents(ht);
        etmp := getv(exponents, 1);
        for i := 1:upbv(etmp) do
            putv(etmp, i, 0);
        
        matrix_reinitialize_matrix(matrix, basis_bget_nlead(basis));
        uprows := matrix_mget_uprows(matrix);

        nonred := basis_bget_nonred(basis);
        gens := basis_bget_gens(basis);
        up2coef := matrix_mget_up2coef(matrix);
        symdata := hashtable_htget_hashdata(symbol_ht);

        for i := 1:basis_bget_nlead(basis) do <<
            nrows := matrix_mget_nrows(matrix) + 1;
            matrix_mset_nrows(matrix, nrows);

            putv(nrows, hashtable_multiplied_poly_to_matrix_row(symbol_ht, ht, 0, etmp, getv(gens, getv(nonred, i))));

            putv(up2coef, nrows, getv(nonred, i));
            hashtable_hvget_idx(getv(symdata, getv(getv(uprows, nrows), 1)), 1)
        >>;

        matrix_mset_ncols(matrix, matrix_mget_nrows(matrix));
        matrix_mset_nup(matrix, matrix_mget_nrows(matrix));

        f4_symbolic_preprocessing(basis, matrix, ht, symbol_ht);
        for i := hashtable_htget_offset(symbol_ht):hashtable_htget_load(symbol_ht) do
            hashtable_hvset_idx(getv(symdata, i), 1);

        matrix_convert_hashes_to_columns(matrix, symbol_ht);
        matrix_mset_ncols(matrix, matrix_mget_nleft(matrix) + matrix_mget_nright);

        % Julia
        sorting_sort_matrix_rows_decreasing(matrix);

        % Julia
        matrix_interreduce_matrix_rows(ring, matrix, basis);

        matrix_convert_matrix_rows_to_basis_elements(matrix, basis, ht, symbol_ht);

        basis_bset_ntotal(basis, matrix_mget_npivots(matrix) + basis_bset_ndone(basis));
        basis_bset_ndone(basis, matrix_mget_npivots(matrix));
        
        k := 0;
        i := 1;
        gens := basis_bget_gens(basis);
        nonred := basis_bget_nonred(basis);
        ntotal := basis_bget_ntotal(basis);
        lead := basis_bget_lead(basis);
        htdata := hashtable_htget_hashdata(ht);
    Letsgo:
        while i <= basis_bget_ndone(basis) do <<
            for j := 1:k do
                if hashtable_is_monom_divisible(
                        getv(getv(gens, ntotal - i + 1), 1),
                        getv(getv(gens, getv(nonred, j)), 1),
                        ht) then <<
                    i := i + 1;
                    go to Letsgo
                >>
            >>;
            k := k + 1;
            putv(nonred, k), ntotal - i + 1);
            putv(lead, k, hashtable_hvget_divmask(getv(htdata, getv(getv(gens, getv(nonred, k)), 1)));
            i := i + 1
        >>
        basis_bset_nlead(basis, k)

    end;

%--------------------------------------------------------------------------------------------------

asserted procedure f4_find_multiplied_reducer(basis: Basis, matrix: MacaulayMatrix, 
                                ht: MonomialHashtable, symbol_ht: MonomialHashtable,
                                vidx: Integer)
    begin scalar e;
        symexps := hashtable_htget_exponents(symbol_ht);
        symdata := hashtable_htget_hashdata(symbol_ht);

        e := getv(symexps, vidx);
        etmp := getv(exponents, 1);
        divmask := hashtable_hvget_divmask(getv(symdata, vidx));

        blen := basis_bget_ndone(basis);
        leaddiv := basis_bget_lead(basis);

        i := 1;
    Letsgo:
        

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
asserted procedure f4_f4(ring: PolyRing, basis: Basis, ht: MonomialHashtable);
    begin scalar pairset, matrix, update_htm symbol_ht, plcm;

        ASSERT(io_prget_ord(ring) = hashtable_htget_ord(ht));
        ASSERT(io_prget_nvars(ring) = hashtable_htget_nvars(ht));
        ASSERT(io_prget_explen(ring) = hashtable_htget_explen(ht));

        ASSERT(basis_bget_ndone(basis) = 0);

        pairset := basis_initialize_pairset();

        % matrix storing coefficients in rows
        % wrt columns representing the current monomial basis
        matrix := matrix_initialize_matrix(ring);

        % initialize hash tables for update and symbolic preprocessing steps
        update_ht := hash_initialize_secondary_hashtable(ht);
        symbol_ht := hash_initialize_secondary_hashtable(ht);
        
        plcm := dv_undef(0);
        core_update(pairset, basis, ht, update_ht, plcm);
        
        d := 0;
        while not pairs_isempty(pairset) do <<
            d := d + 1;
            prin2t {"Iter", d};

            core_select_normal(pairset, basis, matrix, ht, symbol_ht);

            core_symbolic_preprocessing(basis, matrix, ht, symbol_ht);

            core_reduction(ring, basis, matrix, ht, symbol_ht)

            core_update(pairset, basis, ht, update_ht, plcm);

            matrix := matrix_initialize_matrix(ring);
            symbol_ht := hash_initialize_secondary_hashtable(ht);
        >>;

        core_filter_redundant(basis);
        core_reducegb(ring, basis, matrix, ht, symbol_ht);
        core_standardize_basis(ring, basis, ht);

    end;

endmodule; % end of module f4;

end; % end of file