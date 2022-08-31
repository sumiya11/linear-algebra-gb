module f4coeffs;
% Polynomial coefficients manipulations

% Julia: no CoeffBuffer

asserted procedure coeffs_CoeffAccum(): CoeffAccum;
    coeffs_CoeffAccum1(dv_undef(0), dv_undef(0), dv_undef(0));

asserted procedure coeffs_CoeffAccum1(gb_coeffs_zz: Vector, prev_gb_coeffs_zz: Vector, gb_coeffs_qq: Vector): CoeffAccum;
    begin scalar v;
        v := dv_undef(3);
        putv(v, 1, gb_coeffs_zz);
        putv(v, 2, prev_gb_coeffs_zz);
        putv(v, 3, gb_coeffs_qq);
        return v;
    end;

asserted procedure coeffs_caget_gb_coeffs_zz(ca: CoeffAccum): Vector;
    getv(ca, 1);

asserted procedure coeffs_caget_prev_gb_coeffs_zz(ca: CoeffAccum): Vector;
    getv(ca, 2);

asserted procedure coeffs_caget_gb_coeffs_qq(ca: CoeffAccum): Vector;
    getv(ca, 3);

asserted procedure coeffs_caset_gb_coeffs_zz(ca: CoeffAccum, x): Vector;
    putv(ca, 1, x);

asserted procedure coeffs_caset_prev_gb_coeffs_zz(ca: CoeffAccum, x): Vector;
    putv(ca, 2, x);

asserted procedure coeffs_caset_gb_coeffs_qq(ca: CoeffAccum, x): Vector;
    putv(ca, 3, x);

%--------------------------------------------------------------------------------------------------

asserted procedure coeffs_common_denominator(coeffs: Vector): Integer;
    begin scalar den; 
        den := 1;
        for i := 1:dv_length(coeffs) do
            den := lcm(den, denr coeffs[i]);
        return den
    end;

asserted procedure coeffs_scale_denominators1(coeffs_qq: Vector, coeffs_zz: Vector): Vector;
    begin scalar den, num, buf;

        ASSERT(dv_length(coeffs_qq) = dv_length(coeffs_zz));
        
        for i := 1:dv_length(coeffs_qq) do <<
            ASSERT(dv_length(coeffs_qq[i]) = dv_length(coeffs_zz[i]));

            den := coeffs_common_denominator(coeffs_qq[i]);
            for j := 1:dv_length(coeffs_qq[i]) do <<
                num := numr(coeffs_qq[i, j]);
                buf := den / denr(coeffs_qq[i, j]);
                coeffs_zz[i, j] := num * buf
            >>
        >>;

        return coeffs_zz
    end;

asserted procedure coeffs_scale_denominators(coeffs_qq: Vector): Vector;
    begin scalar coeffs_zz;
        coeffs_zz := dv_undef(dv_length(coeffs_qq));
        for i := 1:dv_length(coeffs_qq) do
            putv(coeffs_zz, i, dv_zeros(dv_length(coeffs_qq[i])));
        return coeffs_scale_denominators1(coeffs_qq, coeffs_zz)
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure coeffs_reduce_modulo(coeffs_zz: Vector, coeffs_ff: Vector, prime);
    begin scalar p, buf, c, cfs_zz_i;
        p := prime;
        for i := 1:dv_length(coeffs_zz) do <<
            cfs_zz_i := coeffs_zz[i];
            for j := 1:dv_length(cfs_zz_i) do <<
                c := cfs_zz_i[j];
                if c < 0 then <<
                    buf := c / p;
                    if remainder(c, p) < 0 then
                        buf := buf - 1;
                    buf := -buf;
                    buf := buf*prime;
                    c := c + buf
                >>;
                ASSERT(c >= 0);
                buf := remainder(c, p);
                coeffs_ff[i, j] := buf
            >>
        >>;
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure coeffs_resize_accum(coeffaccum: CoeffAccum, gb_coeffs);
    <<
        coeffs_caset_gb_coeffs_zz(coeffaccum, dv_resize(coeffs_caget_gb_coeffs_zz(coeffaccum), dv_length(gb_coeffs)));
        coeffs_caset_prev_gb_coeffs_zz(coeffaccum, dv_resize(coeffs_caget_prev_gb_coeffs_zz(coeffaccum), dv_length(gb_coeffs)));
        coeffs_caset_gb_coeffs_qq(coeffaccum, dv_resize(coeffs_caget_gb_coeffs_qq(coeffaccum), dv_length(gb_coeffs)));
        for i := 1:dv_length(gb_coeffs) do <<
            putv(coeffs_caget_gb_coeffs_zz(coeffaccum), i, dv_zeros(dv_length(getv(gb_coeffs, i))));
            putv(coeffs_caget_prev_gb_coeffs_zz(coeffaccum), i, dv_zeros(dv_length(getv(gb_coeffs, i))));
            putv(coeffs_caget_gb_coeffs_qq(coeffaccum), i, dv_zeros_sq(dv_length(getv(gb_coeffs, i))))
        >>
    >>;

asserted procedure coeffs_reconstruct_trivial_crt(coeffaccum: CoeffAccum, gb_coeffs_ff);
    begin scalar gb_coeffs_zz;
        gb_coeffs_zz := coeffs_caget_gb_coeffs_zz(coeffaccum);
        for i := 1:dv_length(gb_coeffs_ff) do <<
            for j := 1:dv_length(gb_coeffs_ff[i]) do <<
                ASSERT(gb_coeffs_ff[i, j] > 0);
                gb_coeffs_zz[i, j] := gb_coeffs_ff[i, j]
            >>
        >>
    end;

asserted procedure coeffs_reconstruct_crt(coeffaccum: CoeffAccum, primetracker: PrimeTracker, gb_coeffs_ff, ch: Integer);
    begin scalar buf, n1, n2, M, bigch, invm1, invm2, ca, cf;
        if dv_isempty(coeffs_caget_gb_coeffs_qq(coeffaccum)) then <<
            coeffs_resize_accum(coeffaccum, gb_coeffs_ff);
            coeffs_reconstruct_trivial_crt(coeffaccum, gb_coeffs_ff)
        >> else <<
            gb_coeffs_zz := coeffs_caget_gb_coeffs_zz(coeffaccum);
            prev_gb_coeffs_zz := coeffs_caget_prev_gb_coeffs_zz(coeffaccum);

            % copy to previous gb coeffs
            for i := 1:dv_length(gb_coeffs_zz) do
                for j := 1:dv_length(gb_coeffs_zz[i]) do
                    prev_gb_coeffs_zz[i, j] := gb_coeffs_zz[i, j];
            
            bigch := ch;
            M := lucky_ptget_modulo(primetracker) * ch;
            
            % Julia
            {buf, invm1, invm2} := modular_gcdext(lucky_ptget_modulo(primetracker), bigch);

            for i := 1:dv_length(gb_coeffs_ff) do <<
                for j := 1:dv_length(gb_coeffs_ff[i]) do <<
                    ca := gb_coeffs_zz[i, j];
                    cf := gb_coeffs_ff[i, j];
                    buf := modular_CRT(M, ca, invm1, cf, invm2, lucky_ptget_modulo(primetracker), bigch);
                    ASSERT(buf > 0);
                    gb_coeffs_zz[i, j] := buf
                >>
            >>
        >>;
        lucky_updatemodulo(primetracker)
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure coeffs_reconstruct_modulo(coeffaccum: CoeffAccum, primetracker: PrimeTracker): Boolean;
    begin scalar modulo, bnd, gb_coeffs_zz, gb_coeffs_qq, result, buf;

        modulo := lucky_ptget_modulo(primetracker);

        bnd := modular_rational_reconstruction_bound(modulo);

        gb_coeffs_zz := coeffs_caget_gb_coeffs_zz(coeffaccum);
        gb_coeffs_qq := coeffs_caget_gb_coeffs_qq(coeffaccum);

        result := t;
        for i := 1:dv_length(gb_coeffs_zz) do <<
            for j := 1:dv_length(gb_coeffs_zz[i]) do <<
                cz := gb_coeffs_zz[i, j];
                cq := gb_coeffs_qq[i, j];
                % num := numr(cq);
                % den := denr(cq);
                success . buf := modular_rational_reconstruction(bnd, cz, modulo);
                gb_coeffs_qq[i, j] := buf;
                if not success then <<
                    result := nil;
                    go to Return_
                >>
            >>
        >>;

    Return_:
        return result
    end;


endmodule; % end of f4coeffs

end; % eof