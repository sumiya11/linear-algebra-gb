module f4groebner;
% Groebner basis computation. Contains a wrapper for f4 algorithm

%--------------------------------------------------------------------------------------------------

% Select the optimal size of hashtable
asserted procedure groebner_select_tablesize(ring: PolyRing, exps: Vector): Integer;
    begin scalar nvars;
        % Julia: for now, we want only 8 elements for debugging purposes
        nvars := io_prget_nvars(ring);
        return 8
    end;

%--------------------------------------------------------------------------------------------------

% Computes the Groebner basis.
%   . ring - a polynomial ring 
%   . exps - a list of polynomials' terms (not hashed)
%   . coeffs - a list of polynomials' coefficients
asserted procedure groebner_groebner(ring: PolyRing, exps: Vector, coeffs: Vector): List;
    begin scalar tablesize, basis, ht, gbexps;
        
        % select hashtable size
        tablesize := groebner_select_tablesize(ring, exps);
        
        if f4_debug() then
            prin2t {"Hashtable size: ", tablesize};

        basis . ht := f4_initialize_structures(ring, exps, coeffs, tablesize);
        
        f4_f4(ring, basis, ht, nil);
        
        gbexps := basis_hash_to_exponents(basis, ht);

        if f4_debug() then <<
            prin2t {"After f4: ", basis};
            prin2t {"exps : ", gbexps}
        >>;

        return gbexps . basis_bget_coeffs(basis)
    end;

endmodule; % end of groebner module

end; % of file