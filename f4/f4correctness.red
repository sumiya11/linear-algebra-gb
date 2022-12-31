module f4correctness;
% Checking correctness of Groebner basis over rationals.
% This file corresponds to file gb/correctness.jl in Groebner.jl

revision('f4correctness, "$Id$");

copyright('f4correctness, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%    * Redistributions of source code must retain the relevant
%      copyright notice, this list of conditions and the following
%      disclaimer.
%    * Redistributions in binary form must reproduce the above
%      copyright notice, this list of conditions and the following
%      disclaimer in the documentation and/or other materials provided
%      with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

% Calls all existing correctness checks for a reconstructed basis stored in `coeffaccum`
asserted procedure correctness_correctness_check(coeffaccum, primetracker, ring, exps, coeffs, coeffs_zz, gens_temp_ff, gb_ff, ht): Boolean;
    begin
        % first we check coefficients only
        if not correctness_heuristic_correctness_check(coeffs_caget_gb_coeffs_qq(coeffaccum), lucky_ptget_modulo(primetracker)) then
            return nil;

        if not correctness_randomized_correctness_check(coeffaccum, ring, coeffs_zz, gens_temp_ff, gb_ff, primetracker, ht) then
            return nil;

        return t
    end;

% Gives a threshold that determines if the rational number with: 
%  - `sznum` bitsize of numerator,
%  - `szden` bitsize of denominator,
%  - `szmod` bitsize of modulo
% is unlikely to be a coefficient in a correct groebner basis.  
asserted procedure correctness_threshold_in_heuristic(sznum: Integer, szden: Integer, szmod: Integer): Boolean;
    1.30 * (sznum + sznum) >= szmod;

% Integer logarithm of x base 2
asserted procedure correctness_ilog2(x: Integer);
    begin integer i;
        i := 1;
        while 2^i <= x do i := i + 1;
        return i-1
    end;

% bitsize of modulo
asserted procedure correctness_sizeinbase2(modulo: Integer): Integer;
    correctness_ilog2(abs(modulo));

% Heuristic correctness check: checks that all coefficients of the reconstructed basis
% agree with the threshold in `correctness_threshold_in_heuristic`
asserted procedure correctness_heuristic_correctness_check(gbcoeffs_qq: Vector, modulo: Integer): Boolean;
    begin scalar lnm, result, n, d;
        
        lnm := correctness_sizeinbase2(modulo);

        result := t;
        for i := 1:dv_length(gbcoeffs_qq) do <<
            for j := 1:dv_length(getv(gbcoeffs_qq, i)) do <<
                n := numr(gbcoeffs_qq[i, j]);
                d := denr(gbcoeffs_qq[i, j]);
                if correctness_threshold_in_heuristic(correctness_sizeinbase2(n), correctness_sizeinbase2(d), lnm) then <<
                    result := nil;
                    go to Return_
                >>
            >>
        >>;
    
    Return_:
        return result
    end;

% Randomized correctness check: checks that the reconstructed groebner basis
% is also a groebner basis modulo a random prime number 
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
        normalform_normal_form_f4(ring, gb_ff_copy, ht, gens_ff_copy);

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