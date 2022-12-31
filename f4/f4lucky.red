module f4lucky;
% Lucky prime numbers.
% This file corresponds to file gb/lucky.jl in Groebner.jl

revision('f4lucky, "$Id$");

copyright('f4lucky, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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

% The helper module to keep track of prime numbers used in F4 modular computations.
% Provides the structure `Primetracker` to store intermediate "lucky" and "good" primes.
% We call the prime "lucky" if it is used directly for modular computation.
% We call the ptime "good" if it is used to verify the correctness of
% the modular computation.

% Dev comment: maybe think of our naming convention for "good" prime ? :)
% Strictly speaking, "lucky" prime is a conventional name, while "good" prime is not.

% These constants are set in f4constants.red.
% Lucky and good primes at any moment should not exceed these constants
fluid '(f4_largest!-small!-prime!*
        f4_largest!-small!-modulus!*);

% struct PrimeTracker
%
% PrimeTracker keeps track of primes for modular computation in f4.
asserted procedure lucky_PrimeTracker(coeffs: Vector): PrimeTracker;
    begin scalar firstprime;
        firstprime := isqrt(f4_largest!-small!-modulus!*) / 2;
        firstprime := nextprime(firstprime);
        return lucky_PrimeTracker1(coeffs, firstprime, nextprime(firstprime / 2), dv_undef(0), 1)
    end;

asserted procedure lucky_PrimeTracker1(coeffs: Vector, luckyprime: Integer, goodprime: Integer, primes: Vector, modulo: Integer): PrimeTracker;
    begin scalar v;
        v := dv_undef(5);
        putv(v, 1, coeffs);
        putv(v, 2, luckyprime);
        putv(v, 3, goodprime);
        putv(v, 4, primes);
        putv(v, 5, modulo);
        return v
    end;

asserted procedure lucky_ptget_coeffs(pt: PrimeTracker): Vector;
    getv(pt, 1);

asserted procedure lucky_ptget_luckyprime(pt: PrimeTracker): Integer;
    getv(pt, 2);

asserted procedure lucky_ptget_goodprime(pt: PrimeTracker): Integer;
    getv(pt, 3);

asserted procedure lucky_ptget_primes(pt: PrimeTracker): Vector;
    getv(pt, 4);

asserted procedure lucky_ptget_modulo(pt: PrimeTracker): Integer;
    getv(pt, 5);

asserted procedure lucky_ptset_coeffs(pt: PrimeTracker, x): Vector;
    putv(pt, 1, x);

asserted procedure lucky_ptset_luckyprime(pt: PrimeTracker, x): Integer;
    putv(pt, 2, x);

asserted procedure lucky_ptset_goodprime(pt: PrimeTracker, x): Integer;
    putv(pt, 3, x);

asserted procedure lucky_ptset_primes(pt: PrimeTracker, x): Vector;
    putv(pt, 4, x);

asserted procedure lucky_ptset_modulo(pt: PrimeTracker, x): Integer;
    putv(pt, 5, x);

% Check that a given prime is lucky!
asserted procedure lucky_isluckyprime(tracker: PrimeTracker, prime: Integer): Boolean;
    begin scalar p, coeffs, poly, c, result;
        result := t;
        p := prime;
        coeffs := lucky_ptget_coeffs(tracker);
        for i := 1:dv_length(coeffs) do <<
            poly := getv(coeffs, i);
            for j := 1:dv_length(poly) do <<
                c := getv(poly, j);
                if mod(c, p) = 0 then <<
                    result := nil;
                    go to Return_ 
                >>
            >>
        >>;
    Return_:
        return result
    end;

% Returns the next lucky prime and updates the internal state of PrimeTracker
asserted procedure lucky_nextluckyprime(tracker: PrimeTracker): Integer;
    begin scalar newsz, prime;
        prime := lucky_ptget_luckyprime(tracker);

        while not lucky_isluckyprime(tracker, prime) do
            % Julia: nextprime(5) is 5
            % Reduce: nextprime(5) is 7
            prime := nextprime(prime);
        
        lucky_ptset_luckyprime(tracker, nextprime(prime));

        % Julia: push
        newsz := dv_length(lucky_ptget_primes(tracker)) + 1;
        lucky_ptset_primes(tracker, dv_resize(lucky_ptget_primes(tracker), newsz));
        putv(lucky_ptget_primes(tracker), newsz, prime);

        return prime
    end;

% Updates the product of all primes in PrimeTracker
asserted procedure lucky_updatemodulo(tracker: PrimeTracker);
    lucky_ptset_modulo(tracker, lucky_ptget_modulo(tracker)*f4_getvlast(lucky_ptget_primes(tracker)));

% Returns the next good prime and updates the internal state of PrimeTracker
asserted procedure lucky_nextgoodprime(tracker: PrimeTracker): Integer;
    begin scalar newsz, prime;
        prime := lucky_ptget_luckyprime(tracker);

        while not lucky_isluckyprime(tracker, prime) do
            % Julia: nextprime(5) is 5
            % Reduce: nextprime(5) is 7
            prime := nextprime(prime);
        
        lucky_ptset_goodprime(tracker, nextprime(prime));

        % Julia: push
        newsz := dv_length(lucky_ptget_primes(tracker)) + 1;
        lucky_ptset_primes(tracker, dv_resize(lucky_ptget_primes(tracker), newsz));
        putv(lucky_ptget_primes(tracker), newsz, prime);

        return prime
    end;

endmodule; % end of lucky

end; % end of file