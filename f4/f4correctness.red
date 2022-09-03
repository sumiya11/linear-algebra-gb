module f4correctness;

asserted procedure correctness_correctness_check(coeffaccum, primetracker, ring, exps, coeffs, coeffs_zz, gens_temp_ff, gb_ff, ht): Boolean;
    begin
        
        % first we check coefficients only
        if not correctness_heuristic_correctness_check(coeffs_caget_gb_coeffs_qq(coeffaccum), lucky_ptget_modulo(primetracker)) then
            return nil;

        if not correctness_randomized_correctness_check(coeffaccum, ring, coeffs_zz, gens_temp_ff, gb_ff, primetracker, ht) then
            return nil;

        return t
    end;

asserted procedure correctness_threshold_in_heuristic(sznum: Integer, szden: Integer, szmod: Integer): Boolean;
    1.30 * (sznum + sznum) >= szmod;

asserted procedure correctness_ilog2(x: Integer);
    begin integer i;
        i := 1;
        while 2^i <= x do i := i + 1;
        return i-1
    end;

asserted procedure correctness_sizeinbase(modulo: Integer): Integer;
    correctness_ilog2(abs(modulo));

asserted procedure correctness_heuristic_correctness_check(gbcoeffs_qq: Vector, modulo: Integer): Boolean;
    begin scalar lnm, result, n, d;
        
        lnm := correctness_sizeinbase(modulo);

        result := t;
        for i := 1:dv_length(gbcoeffs_qq) do <<
            for j := 1:dv_length(getv(gbcoeffs_qq, i)) do <<
                n := numr(gbcoeffs_qq[i, j]);
                d := denr(gbcoeffs_qq[i, j]);
                if correctness_threshold_in_heuristic(correctness_sizeinbase(n), correctness_sizeinbase(d), lnm) then <<
                    result := nil;
                    go to Return_
                >>
            >>
        >>;
    
    Return_:
        return result
    end;

asserted procedure correctness_randomized_correctness_check(coeffaccum, ring, coeffs_zz, gens_temp_ff, gb_ff, primetracker, ht): Boolean;
    begin scalar gens_ff_copy, gb_coeffs_zz, gb_ff_copy, result, goodprime;
        result := t;
        
        goodprime := lucky_nextgoodprime(primetracker);

        coeffs_reduce_modulo(coeffs_zz, basis_bget_coeffs(gens_temp_ff), goodprime);
        gens_ff_copy := basis_copy_basis_thorough(gens_temp_ff);
        groebner_cleanup_gens(ring, gens_ff_copy, goodprime);

        gb_coeffs_zz := coeffs_scale_denominators(coeffs_caget_gb_coeffs_qq(coeffaccum));
        gb_ff_copy := basis_copy_basis_thorough(gb_ff);
        coeffs_reduce_modulo(gb_coeffs_zz, basis_bget_coeffs(gb_ff_copy), goodprime);
        groebner_cleanup_gens(ring, gb_ff_copy, goodprime);

        % prin2t {"Prime in correctness:", goodprime};

        % prin2t {"Before normal form:"};
        % prin2t {"\tRing:", ring};
        % prin2t {"\tGens:", gens_ff_copy};
        % prin2t {"\tGB:", gb_ff_copy};

        normalform_normal_form_f4(ring, gb_ff_copy, ht, gens_ff_copy);

        % prin2t {"After normal form:", gens_ff_copy};

        for i := 1:basis_bget_ndone(gens_ff_copy) do <<
            if not dv_isempty(getv(basis_bget_coeffs(gens_ff_copy), i)) then <<
                result := nil;
                go to Return_
            >>
        >>;

        if not isgroebner_isgroebner_f4(ring, gb_ff_copy, ht) then <<
            result := nil;
            go to Return_
        >>;

    Return_:
        return result
    end;


endmodule; % end of module

end; % eof