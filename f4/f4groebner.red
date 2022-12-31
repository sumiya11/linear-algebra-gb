module f4groebner;
% Groebner basis computation. Contains a wrapper for f4 algorithm.
% This file corresponds to file gb/groebner.jl in Groebner.jl

revision('f4groebner, "$Id$");

copyright('f4groebner, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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

% Selectы the optimal size of hashtable
asserted procedure groebner_select_tablesize(ring: PolyRing, exps: Vector): Integer;
    2^15; % Hardcoded to 2^15, which proved to be good in almost all cases.

asserted procedure groebner_cleanup_gens(ring, gens_ff, prime);
    <<
        io_prset_ch(ring, prime);
        basis_normalize_basis(ring, gens_ff)
    >>;

% Computes the Groebner basis over Z/Zp or Q using direct computation 
%   . ring - polynomial ring (see f4io.red),
%   . exps - list of polynomials' terms,
%   . coeffs - a list of polynomials' coefficients.
asserted procedure groebner_groebner1(ring: PolyRing, exps: Vector, coeffs: Vector): DottedPair;
    begin scalar tablesize, basis, ht, gbexps;
        
        % select hashtable size
        tablesize := groebner_select_tablesize(ring, exps);

        basis . ht := f4_initialize_structures(ring, exps, coeffs, tablesize);
        
        f4_f4(ring, basis, ht, t);
        
        gbexps := basis_hash_to_exponents(basis, ht);

        return gbexps . basis_bget_coeffs(basis)
    end;

% Computes the Groebner basis over Q using modular algorithmю
% The strategy is incremental. It is implemented roughly in the following way:
%
% k = 1
% while !(correctly reconstructed)
%   k = k*2
%   select a batch of "register-size" prime numbers p1..pk
%   compute a batch of finite field groebner bases gb1..gbk
%   reconstruct gb1..gbk to gb_zz modulo prod(p1..pk) with CRT
%   reconstruct gb_zz to gb_qq with rational reconstruction
%   if the basis gb_qq is correct, then break
% end
% return gb_qq
%
asserted procedure groebner_groebner2(ring: PolyRing, exps: Vector, coeffs: Vector): DottedPair;
    begin scalar tablesize, gens_temp_ff, ht, coeffaccum, correct,
                coeffs_zz, primetracker, i, gens_ff, prime, gap,
                primegaps, primemult, gb_exps, success;
        
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
        % ring.ch will be nonzero !
        f4_f4(ring, gens_ff, ht, t);
        
        % reconstruct into integers
        coeffs_reconstruct_crt(coeffaccum, primetracker, basis_bget_coeffs(gens_ff), prime);

        % reconstruct into rationals
        success := coeffs_reconstruct_modulo(coeffaccum, primetracker);

        if success and correctness_correctness_check(coeffaccum, primetracker, ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht) then
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
                
                % copy basis so that initial exponents dont get lost
                gens_ff := basis_copy_basis_thorough(gens_temp_ff);
                prime := lucky_nextluckyprime(primetracker);
                % perform reduction and store result in gens_ff
                coeffs_reduce_modulo(coeffs_zz, basis_bget_coeffs(gens_ff), prime);
                % do some things to ensure generators are correct
                groebner_cleanup_gens(ring, gens_ff, prime);

                %  compute groebner basis in finite field
                f4_f4(ring, gens_ff, ht, t);

                coeffs_reconstruct_crt(coeffaccum, primetracker, basis_bget_coeffs(gens_ff), prime)
            >>;
            
            % reconstruct into rationals
            success := coeffs_reconstruct_modulo(coeffaccum, primetracker);

            if success and correctness_correctness_check(coeffaccum, primetracker, ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht) then
                correct := t
        >>;

        gb_exps := basis_hash_to_exponents(gens_ff, ht);
        return gb_exps . coeffs_caget_gb_coeffs_qq(coeffaccum)
    end;

endmodule; % end of groebner module

end; % of file