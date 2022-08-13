module groebner;
% Wrapper for f4

% Computes the Groebner basis.
%   . ring - a polynomial ring 
%   . exps - a list of polynomials' terms
%   . coeffs - a list of polynomials' coefficients
asserted procedure groebner_groebner(ring: PolyRing, exps: List, coeffs: List): List;
    begin scalar tablesize, basis, ht, gbexps;
        % Julia: for now, we want only 16 elements for debugging
        tablesize := 16;
        
        basis . ht := f4_initialize_structures(ring, exps, coeffs, tablesize);
        
        f4_f4(ring, basis, ht);
        
        gbexps := f4_hash_to_exponents(basis, ht);
        return gbexps . basis_bget_coeffs(basis)
    end;

endmodule; % end of groebner module

end; % of file