module f4;

asserted procedure f4_initializeStructures(ring, polys);
    begin;
        basis := basis_initializeBasis(ring);
        ht := hash_initializeHashtable(ring);
        basis_fillData(basis, ht, polys);
        % hash_fillDivmask(ht);
        sort_sortByLeadInc(basis, ht);
        basis_normalize(ring, basis);
        return basis . ht
    end;

asserted procedure f4_groebner1(ring, basis, ht);
    begin;
        matrix := matrix_initialize_matrix(ring);
        update_ht := hash_initialize_secondary_hashtable(ht);
        symbol_ht := hash_initialize_secondary_hashtable(ht);
        pairset := pairs_initialize_pairset();
        
        plcm := mkvect(0);
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