module groebner;

% List of Poly --> List of Poly
asserted procedure groebner_groebner(ring, exps, coeffs): List;
    begin;
        basis . ht := f4_initialize_structures(ring, exps, coeffs);
        f4_f4(ring, basis, ht);
        gbexps := f4_hash_to_exponents(basis, ht);
        return gbexps . basis_bget_coeffs(basis)
    end;

endmodule; % end of groebner module

end; % of file