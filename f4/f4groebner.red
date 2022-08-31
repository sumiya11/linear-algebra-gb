module f4groebner;
% Groebner basis computation. Contains a wrapper for f4 algorithm

%--------------------------------------------------------------------------------------------------

% Select the optimal size of hashtable
asserted procedure groebner_select_tablesize(ring: PolyRing, exps: Vector): Integer;
    begin scalar nvars;
        % Julia: for now, we want only 2^14 elements for debugging purposes
        nvars := io_prget_nvars(ring);
        return 2^14  % 2^10 - 2^16 in Julia
    end;

asserted procedure groebner_cleanup_gens(ring, gens_ff, prime);
    <<
        io_prset_ch(ring, prime);
        basis_normalize_basis(ring, gens_ff)
    >>;

%--------------------------------------------------------------------------------------------------

% Computes the Groebner basis.
%   . ring - a polynomial ring 
%   . exps - a list of polynomials' terms (not hashed)
%   . coeffs - a list of polynomials' coefficients
asserted procedure groebner_groebner1(ring: PolyRing, exps: Vector, coeffs: Vector): DottedPair;
    begin scalar tablesize, basis, ht, gbexps;
        
        % select hashtable size
        tablesize := groebner_select_tablesize(ring, exps);

        basis . ht := f4_initialize_structures(ring, exps, coeffs, tablesize);
        
        f4_f4(ring, basis, ht, t);
        
        gbexps := basis_hash_to_exponents(basis, ht);

        return gbexps . basis_bget_coeffs(basis)
    end;

asserted procedure groebner_groebner2(ring: PolyRing, exps: Vector, coeffs: Vector): DottedPair;
    begin scalar tablesize, gens_temp_ff, ht, coeffaccum,
                coeffs_zz, primetracker, i, gens_ff, prime;
        
        % select hashtable size
        tablesize := groebner_select_tablesize(ring, exps);

        gens_temp_ff . ht := f4_initialize_structures_ff(ring, exps, coeffs, tablesize);
        gens_ff := basis_copy_basis_thorough(gens_temp_ff);

        % now hashtable is filled correctly,
        % and gens_temp_ff exponents are correct and in correct order.
        % gens_temp_ff coefficients are filled with random stuff and
        % gens_temp_ff.ch is 0

        % to store integer and rational coefficients of groebner basis
        coeffaccum := coeffs_CoeffAccum();

        % scale coefficients of input to integers
        coeffs_zz := coeffs_scale_denominators(coeffs);

        % keeps track of used prime numbers
        primetracker := lucky_PrimeTracker(coeffs_zz);

        i := 1;

        % copy basis so that we initial exponents dont get lost
        gens_ff := basis_copy_basis_thorough(gens_temp_ff);

        prime := lucky_nextluckyprime(primetracker);

        % perform reduction and store result in gens_ff
        coeffs_reduce_modulo(coeffs_zz, basis_bget_coeffs(gens_ff), prime);

        % do some things to ensure generators are correct
        % Julia: normalize, somehow
        groebner_cleanup_gens(ring, gens_ff, prime);

        % compute groebner basis in finite field
        % ring.ch will be nonzero !!!
        f4_f4(ring, gens_ff, ht, t);

        % reconstruct into integers
        coeffs_reconstruct_crt(coeffaccum, primetracker, basis_bget_coeffs(gens_ff), prime);

        % reconstruct into rationals
        coeffs_reconstruct_modulo(coeffaccum, primetracker);

        if correctness_correctness_check(coeffaccum, primetracker, ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht) then
            correct := t;

        gap := 1;
        primegaps := {1,1,1,1,1};
        primemult := 2;

        while not correct do <<
            
            if primegaps then
                gap := pop(primegaps)
            else
                gap := gap*primemult;

            for j := 1:gap do <<
                i := i + 1;

                # copy basis so that initial exponents dont get lost
                gens_ff := basis_copy_basis_thorough(gens_temp_ff);
                prime := lucky_nextluckyprime(primetracker);
                % perform reduction and store result in gens_ff
                coeffs_reduce_modulo(coeffs_zz, basis_bget_coeffs(gens_ff), prime);
                % do some things to ensure generators are correct
                groebner_cleanup_gens(ring, gens_ff, prime);

                %  compute groebner basis in finite field
                f4(ring, gens_ff, ht, t);

                modular_reconstruct_crt(coeffaccum, primetracker, basis_bget_coeffs(gens_ff), prime)
            >>;
            
            % reconstruct into rationals
            modular_reconstruct_modulo(coeffaccum, primetracker);

            if correctness_correctness_check(coeffaccum, primetracker, ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht) then
                correct := t
        >>;

        gb_exps := basis_hash_to_exponents(gens_ff, ht);
        return gb_exps . coeffs_caget_gb_coeffs_qq(coeffaccum)
    end;

endmodule; % end of groebner module

end; % of file