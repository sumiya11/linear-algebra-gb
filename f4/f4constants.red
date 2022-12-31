module f4constants;
% Platform-dependent detection of a largest possible reduction modulo.

revision('f4constants, "$Id$");

copyright('f4constants, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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

fluid '(f4_largest!-small!-modulus!*);
fluid '(f4_largest!-small!-prime!*);

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

endmodule;  % end of f4constants

end; % of file
