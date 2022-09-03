module f4normalform;

asserted procedure normalform_normal_form_f4(ring: PolyRing, basis: Basis, ht: MonomialHashtable, tobereduced: Basis): Basis;
    begin scalar matrixj, symbol_ht;
        matrixj := matrix_initialize_matrix(ring);
        symbol_ht := hashtable_initialize_secondary_hash_table(ht);

        f4_select_tobereduced(basis, tobereduced, matrixj, symbol_ht, ht);

        f4_symbolic_preprocessing(basis, matrixj, ht, symbol_ht);

        matrix_convert_hashes_to_columns(matrixj, symbol_ht);

        sorting_sort_matrix_rows_decreasing(matrixj);

        matrix_linear_algebra_nf(ring, matrixj, tobereduced, basis);

        matrix_convert_nf_rows_to_basis_elements(matrixj, tobereduced, ht, symbol_ht);

        return tobereduced
    end;

endmodule; % end of f4normalform

end; % eof