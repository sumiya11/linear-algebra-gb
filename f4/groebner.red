module groebner;
% Groebner basis computation. Contains a wrapper for f4 algorithm

%--------------------------------------------------------------------------------------------------

% Select the optimal size of hashtable
asserted procedure groebner_select_tablesize(ring: PolyRing, exps: Vector): Integer;
    begin scalar nvars;
        % Julia: for now, we want only 16 elements for debugging purposes
        nvars := io_pr_nvars(ring);
        return 16
    end;

%--------------------------------------------------------------------------------------------------

% Computes the Groebner basis.
%   . ring - a polynomial ring 
%   . exps - a list of polynomials' terms (not hashed)
%   . coeffs - a list of polynomials' coefficients
asserted procedure groebner_groebner(ring: PolyRing, exps: List, coeffs: List): List;
    begin scalar tablesize, basis, ht, gbexps;
        
        % select hashtable size
        tablesize := groebner_select_tablesize(ring, exps);
        
        basis . ht := f4_initialize_structures(ring, exps, coeffs, tablesize);
        
        f4_f4(ring, basis, ht, nil);
        
        gbexps := f4_hash_to_exponents(basis, ht);
        return gbexps . basis_bget_coeffs(basis)
    end;

endmodule; % end of groebner module

end; % of file