module f4io;
% Input-output converions of polynomials in F4.

% Julia: Note that the contents of this file differ significantlyfrom 
%   the one from the Julia implementation.
%   This is caused by the fact that conversion between SQ and our polynomial type
%   differs from Julia (Julia has no concept of SQ).   
% *Still, The interface provided by the file is the same as in Julia:*
%   

%--------------------------------------------------------------------------------------------------

% struct PolyRing - stores information about polynomial ring
% 
% A PolyRing object contains the following:
%   nvars - number of variables
%   explen - length of exponent vectors
%   ord - sort mode

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

%--------------------------------------------------------------------------------------------------

% list of lp --> {ring, vector of vectors of monons, vector of vectors of coeffs}
asserted procedure io_convert_to_internal(polys: List, vars: List, ord: Any): List;
    begin scalar n, gens, coeffs, ring, i, gensi, coeffsi, 
                    gensij, explen, gensij_list;
        
        if !*f4debug then
            prin2t "convert_to_internal..";
        
        vars := cdr vars;

        n := length(polys);
        gens := dv_undef(n);
        coeffs := dv_undef(n);
        ring := io_PolyRing(length(vars), length(vars) + 1, ord, 0);
        
        if !*f4debug then
            prin2t {"convert_to_internal: ring:", ring};

        % Julia!!
        % poly_initRing('list . vars, ord);

        explen := io_prget_explen(ring);

        i := 1;
        for each poly in polys do <<
            
            poly := poly_sq2poly(simp poly);
            gensi := poly_getTerms(poly);
            coeffsi := poly_getCoeffs(poly);

            if !*f4debug then <<
                prin2t {"poly", i};
                prin2t {gensi, coeffsi}
            >>;

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

            sorting_sort_input_to_change_ordering(gensi, coeffsi, 'revgradlex);

            i := i + 1
        >>;
        
        return {ring, gens, coeffs}
    end;

% {ring, vector of vectors of monons, vector of vectors of coeffs} --> list of lp
asserted procedure io_convert_to_output(ring: PolyRing, bexps: Vector, bcoeffs: Vector): List;
    begin scalar anssq, bexpsi, bcoeffsi;
        
        if !*f4debug then
            prin2t "convert_to_output..";

        for i := 1:dv_length(bexps) do <<
            bexpsi := getv(bexps, i);
            bexpsi := for j := 1:dv_length(bexpsi) collect
                (f4_getvlast(getv(bexpsi, j)) . for k := 1:io_prget_explen(ring)-1 collect 
                    getv(getv(bexpsi, j), k));
            bcoeffsi := for j := 1:dv_length(getv(bcoeffs, i)) collect
                getv(getv(bcoeffs, i), j);
            
            if !*f4debug then <<
                prin2t {"poly", i};
                prin2t {bexpsi, bcoeffsi}
            >>;

            push(poly_poly2lp(poly_Polynomial(bexpsi, bcoeffsi)), anssq)
        >>;
        return reversip(anssq)
    end;

endmodule; % end of io module

end; % end of file
