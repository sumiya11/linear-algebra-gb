module f4constants;

fluid '(f4_largest!-small!-modulus!*);
fluid '(f4_largest!-small!-prime!*);

% on1 'assert;

if 'sixty!-four memq lispsystem!* and 'csl memq lispsystem!* then
   f4_largest!-small!-modulus!* := 2^(64-5)                  % 2^59
else if 'sixty!-four memq lispsystem!* and 'psl memq lispsystem!* then
   f4_largest!-small!-modulus!* := 2^(64-8)                  % 2^56
else
   f4_largest!-small!-modulus!* := largest!-small!-modulus;  % 2^23

#if (memq 'psl lispsystem!*)
   load!-package 'compiler;
   ASSERT( intp(f4_largest!-small!-modulus!* - 1) );
#elif (memq 'csl lispsystem!*)
   ASSERT( eq!-safe(f4_largest!-small!-modulus!* - 1) );
#endif

asserted procedure previousprime(p: Integer): Integer;
   % Returns the next prime number smaller than p.
   begin integer pp;
      if eqn(p, 3) then
         return 2;
      if eqn(p, 2) or eqn(p, 1) or eqn(p, 0) then
         return -2;
      pp := if evenp p then p - 1 else p - 2;
      while not primep pp do
         pp := pp - 2;
      return pp
   end;

f4_largest!-small!-prime!* := previousprime(f4_largest!-small!-modulus!*);

endmodule;

end;
