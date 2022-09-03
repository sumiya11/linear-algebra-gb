module f4isgroebner;

asserted procedure isgroebner_isgroebner_f4(ring: PolyRing, basis: Basis, ht: MonomialHashtable): Boolean;
    begin scalar matrixj, symbol_ht, update_ht, pairset, plcm;
        matrixj := matrix_initialize_matrix(ring);
        symbol_ht := hashtable_initialize_secondary_hash_table(ht);
        update_ht := hashtable_initialize_secondary_hash_table(ht);

        pairset := basis_initialize_pairset();

        plcm := dv_undef(0);
        basis_update(pairset, basis, ht, update_ht, plcm);

        if basis_psget_load(pairset) = 0 then
            return t;

        f4_select_isgroebner(pairset, basis, matrixj, symbol_ht);

        f4_symbolic_preprocessing(basis, matrixj, ht, symbol_ht);

        matrix_convert_hashes_to_columns(matrixj, symbol_ht);

        sorting_sort_matrix_rows_increasing(matrixj);
        sorting_sort_matrix_rows_decreasing(matrixj);

        return matrix_linear_algebra_isgroebner(ring, matrixj, basis)
    end;

endmodule; % end of f4isgroebner

end; % eof