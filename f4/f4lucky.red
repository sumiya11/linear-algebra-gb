module f4lucky;
% Lucky prime numbers.

% Lucky primes are prime numbers used for modular computation.
% The sequence of lucky primes starts with FIRST_LUCKY_PRIME!*
fluid '(FIRST_LUCKY_PRIME!*);
FIRST_LUCKY_PRIME!* := 2^21+17;

% Good primes are prime numbers used for checking correctness.
% The sequence of good primes starts with FIRST_GOOD_PRIME!*
fluid '(FIRST_GOOD_PRIME!*);
FIRST_GOOD_PRIME!* := 2^20+7;

% There are TODO primes between
% FIRST_LUCKY_PRIME and FIRST_GOOD_PRIME

asserted procedure lucky_PrimeTracker(coeffs: Vector): PrimeTracker;
    lucky_PrimeTracker1(coeffs, FIRST_LUCKY_PRIME!*, FIRST_GOOD_PRIME!*, dv_undef(0), 1);

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

asserted procedure lucky_nextluckyprime(tracker: PrimeTracker): Integer;
    begin scalar p, coeffs, poly, c, result, newsz;
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

asserted procedure lucky_updatemodulo(tracker: PrimeTracker);
    lucky_ptset_modulo(tracker, lucky_ptget_modulo(tracker)*f4_getvlast(lucky_ptget_primes(tracker)));

asserted procedure lucky_nextgoodprime(tracker: PrimeTracker): Integer;
    begin scalar p, coeffs, poly, c, result, newsz;
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