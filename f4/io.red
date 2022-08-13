module io;
% Input-output of polynomial in F4.

% Julia: Note that the contents of this file differ significantlyfrom 
%   the one from the Julia implementation.
%   This is caused by the fact that conversion between SQ and our polynomial type
%   differs from Julia (Julia has no concept of SQ).   
% *Still, The interface provided by the file is the same as in Julia:*
%   

% global!-dipvars!* is a list of kernels, 
% which currently serve as variables in polynomial ring.
% vdpsortmode!* is a "sort mode" identifier (e.g., lex, revgradlex)  
fluid '(global!-dipvars!*);
fluid '(vdpsortmode!*);

load!-package 'f5;

%--------------------------------------------------------------------------------------------------

% struct PolyRing - stores information about polynomial ring
% 
% A PolyRing object contains the following:
%   nvars - number of variables
%   explen - length of exponent vectors
%   ord - sort mode

asserted procedure io_PolyRing(nvars: Integer, explen: Integer, ord): PolyRing;
    {'pr, nvars, explen, ord};

asserted procedure io_prget_nvars(pr: PolyRing): Integer;
    car cdrn(pr, 1);

asserted procedure io_prget_explen(pr: PolyRing): Integer;
    car cdrn(pr, 2);

asserted procedure io_prget_ord(pr: PolyRing);
    car cdrn(pr, 3);

asserted procedure io_prset_nvars(pr: PolyRing, x);
    car cdrn(pr, 1) := x;

asserted procedure io_prset_explen(pr: PolyRing, x);
    car cdrn(pr, 2) := x;

asserted procedure io_prset_ord(pr: PolyRing, x);
    car cdrn(pr, 3) := x;

%--------------------------------------------------------------------------------------------------

asserted procedure io_init_globals(vars: List, ord);
    <<  
        torder({'list . vars, ord})
        % global!-dipvars!* := vars;
        % vdpsortmode!* := ord
    >>;

%--------------------------------------------------------------------------------------------------

% list of lp --> list of Poly1
asserted procedure io_convert_to_internal(polys: List): List;
    begin scalar x;
        sqpolys := for each poly in polys collect 
                        poly_sq2poly simp poly;
        gens := dv_undef(length(polys));
        coeffs := dv_undef(length(polys));
        for i := 1:length(sqpolys) do <<
            poly := pop(sqpolys);
            gensi := dv_undef(length(poly));
            termsi := poly_getTerms(poly);
            for j := 1:length(poly) do
                putv(gensi, j, f4putin_list2vector(pop(termsi)));
            putv(gens, i, gensi);
            putv(coeffs, i, f4putin_list2vector(poly_getCoeffs(poly)))
        >>;
        ring := io_PolyRing(length(global!-dipvars!*), 
                            length(global!-dipvars!*)+1, 
                            vdpsortmode!*);
        return {ring, gens, coeffs}
    end;

% list of Poly1 --> list of lp
asserted procedure io_convert_to_output(ring: PolyRing, bexps: Vector, bcoeffs: Vector): List;
    begin scalar ans, expsv, expsvv, coeffsv;
        for i := 1:upbv(bexps) do <<
            expsv := f4putin_vector2list(getv(bexps, i));
            expsvv := expsv;
            while expsvv do <<
                car expsvv := f4putin_vector2list(car expsvv);
                expsvv := cdr expsvv
            >>;
            coeffsv := f4putin_vector2list(getv(coeffs, i));
            push(poly_poly2sq(poly_Polynomial(expv, coeffsv)), ans)
        >>;
        return reversip(ans)
    end;

endmodule; % end of io module

end; % end of file
