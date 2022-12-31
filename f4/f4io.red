module f4io;
% Input-output converions of polynomials in F4.
% This file corresponds to file io.jl in Groebner.jl

revision('f4io, "$Id$");

copyright('f4io, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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

% Julia: Note that the contents of this file differ significantly from 
%   the one from the Julia implementation.
%   This is caused by the fact that conversion between SQ and our polynomial type
%   differs from Julia (Julia has no concept of SQ).   

% struct PolyRing - stores information about polynomial ring
% 
% A PolyRing object contains the following:
%   nvars - number of variables
%   explen - length of exponent vectors
%   ord - current sort mode (should coincide with the one in vdpsortmode!*)
%   ch - characteristic of the ground field (0 for Q, > 0 for Z/Zp)
asserted procedure io_PolyRing(nvars: Integer, explen: Integer, ord, ch: Integer): PolyRing;
    begin scalar v;
        v := dv_undef(4);
        putv(v, 1, nvars);
        putv(v, 2, explen);
        putv(v, 3, ord);
        putv(v, 4, ch);
        return v
    end;

asserted procedure io_prget_nvars(pr: PolyRing): Integer;
    getv(pr, 1);

asserted procedure io_prget_explen(pr: PolyRing): Integer;
    getv(pr, 2);

asserted procedure io_prget_ord(pr: PolyRing);
    getv(pr, 3);

asserted procedure io_prget_ch(pr: PolyRing): Integer;
    getv(pr, 4);

asserted procedure io_prset_nvars(pr: PolyRing, x);
    putv(pr, 1, x);

asserted procedure io_prset_explen(pr: PolyRing, x);
    putv(pr, 2, x);

asserted procedure io_prset_ord(pr: PolyRing, x);
    putv(pr, 3, x);

asserted procedure io_prset_ch(pr: PolyRing, x): Integer;
    putv(pr, 4, x);

asserted procedure io_convert_lp_to_poly(inputLPs: List): List; 
    for each lp in inputLPs collect poly_sq2poly(simp lp);

asserted procedure io_convert_poly_to_lp(inputPolys: List): List;
    for each poly in inputPolys collect poly_poly2lp(poly);

% Converts a list of `Polynomial`s to
%    {polynomial ring, vector of vectors of monons, vector of vectors of coeffs}
asserted procedure io_convert_poly_to_internal(polys: List): List;
    begin scalar n, gens, coeffs, ring, i, gensi, coeffsi, 
                    gensij, explen, gensij_list, vars;
        
        vars := cdr global!-dipvars!*;

        n := length(polys);
        gens := dv_undef(n);
        coeffs := dv_undef(n);
        ring := io_PolyRing(length(vars), length(vars) + 1, vdpsortmode!*, 0);

        explen := io_prget_explen(ring);

        i := 1;
        for each poly in polys do <<
            gensi := poly_getTerms(poly);
            coeffsi := poly_getCoeffs(poly);

            gensi := f4_list2vector(gensi);
            coeffsi := f4_list2vector(coeffsi);
            putv(gens, i, gensi);
            putv(coeffs, i, coeffsi);
            for j := 1 : dv_length(gensi) do <<
                gensij_list := getv(gensi, j);
                gensij := dv_undef(explen);
                putv(gensij, explen, pop(gensij_list));
                for k := 1 : explen-1 do
                    putv(gensij, k, pop(gensij_list));
                putv(gensi, j, gensij)
            >>;

            i := i + 1
        >>;
        
        return {ring, gens, coeffs}
    end;

% Converts {polynomial ring, vector of vectors of monons, vector of vectors of coeffs}
% to a list of `Polynomial`s
asserted procedure io_convert_internal_to_poly(ring: PolyRing, bexps: Vector, bcoeffs: Vector): List;
    begin scalar anssq, bexpsi, bcoeffsi;
        for i := 1:dv_length(bexps) do <<
            bexpsi := getv(bexps, i);
            bexpsi := for j := 1:dv_length(bexpsi) collect
                (f4_getvlast(getv(bexpsi, j)) . for k := 1:io_prget_explen(ring)-1 collect 
                    getv(getv(bexpsi, j), k));
            bcoeffsi := for j := 1:dv_length(getv(bcoeffs, i)) collect
                getv(getv(bcoeffs, i), j);
            
            push(poly_Polynomial(bexpsi, bcoeffsi), anssq)
        >>;
        return reversip(anssq)
    end;

endmodule; % end of io module

end; % end of file