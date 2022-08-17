module f4putin;
% f4putin. Main module. The F4 algorithm for computing Groebner bases.
% Contains the AM entry point.

% This package is migrated from Julia package 
% Groebner.jl maintained by Gowda S. and Demin A. 
%   https://github.com/sumiya11/Groebner.jl/

% DEVELOPER NOTE:
% Anywhere in this package sources, comments that start with the word "Julia"
% explain some of our design decisions during the migration of the package to Reduce.

% The f4putin module provides the implementation of the Faugere's F4 algorithm
%   https://doi.org/10.1016/S0022-4049(99)00005-5
%
% The interface contains the operator `f4` with the following signature
%     `f4(polynomials: List, vars: List, sortmode: Any)`
% Where 
%   . `polynomials` is a list of expressions, the ideal generators;
%   . `vars` is a list of Reduce kernels, the polynomial ring variables;
%   . `sortmode` is a sorting mode for monomial comparison ('lex, 'degrevlex)
% The `f4` function returns a list of expressions, which is the Groebner basis 
% of the ideal; each list element is a Lisp Prefix.

% Julia:
%   The file f4putin.red mirrors the file interface.jl from the Julia implementation.
%       (links are not stable)
%       https://github.com/sumiya11/Groebner.jl/blob/master/src/interface.jl
%   All other files have exactly same names as in the Julia implementation.

% Julia:
%   In Reduce, function names are formed as follows:
%   the file name + the Julia name of the function
%
%   So that Julia function "insert_in_hash_table" from the file "hashtable"
%   has name "hashtable_insert_in_hash_table" in Reduce and is 
%   defined in the file "hashtable.red".

% . f4putin - this file, contains package AM entry point with a single function `f4`;
% . groebner - high-level Groebner basis computation; calls to f4 as its backend;
% . f4 - low-level f4 implementation, the heart of the package;
% . hashtable - MonomialHashtable struct and related functions implementation;
% . basis - Pairset and Basis structs implementations, as well as some related functions;
% . matrix - MacaulayMatrix struct implementation; contains linear-algebra-style reduction;
% . internaltypes - type declarations (have no practical effect in Reduce, 
%                   but are left for consistency with Julia implementation)
% . io - input-output conversions of polynomial representations 
create!-package('(f4putin groebner f4 hashtable basis matrix internaltypes io), nil);

%--------------------------------------------------------------------------------------------------

asserted procedure f4putin_debug();
    t;

% Julia: Reduce-specific things;
%   As we want to mirror the behavior of Julia with integer division masks as close as possible,
%   we use generic Integers (not machine integers) and operate on the lowest `f4putin_sizeofInt32` bits.
fluid '(f4putin_sizeofInt32!*);
f4putin_sizeofInt32!* := 32;   

% Julia: analogue of x[end] -- access the last element of vector x
asserted procedure f4putin_getvlast(x);
    getv(x, length(x));

% Julia: convert list `x` to vector and return it
asserted procedure f4putin_list2vector(x: List): Vector;
    begin scalar len, v;
        len := length(x);
        v := dv_undef(len);
        for i := 1:len do
            putv(v, i, pop(x));
        return v
    end;

% Julia: convert vector `x` to list and return it
asserted procedure f4putin_vector2list(x: Vector): List;
    begin scalar v;
        for i := 1:upbv(x) do
            push(getv(x, i), v);
        return reversip(v)
    end;

%--------------------------------------------------------------------------------------------------

put('f4, 'psopfn, 'f4putin_groebner);

% AM entry point
asserted procedure f4putin_groebner(u: List): List;
    begin scalar polynomials, ring, exps, coeffs, 
                    bexps, bcoeffs;
        if null u or not (listp u) then
            f4putin_argumentError();
        
        % Extract the list of polynomials
        polynomials := reval pop u;
        if not (listp polynomials) or not (pop polynomials eq 'list) or null polynomials then
            f4_argumentError();
        if null u then
            f4_argumentError();
        
        % Extract the list of variables and the sort mode
        % to initialize the polynomial ring
        % variables and sort mode are specified in f4 call
        vars := reval pop u;
        if not (listp vars) or not (pop vars eq 'list) then
            f4_argumentError();
        for each var in vars do
            if not sfto_kernelp(var) then
                f4_argumentError();
        ord := pop u;
        % initialize global variables in io module 
        % for follow-up polynomial parsing
        io_init_globals(vars, ord);

        {ring, exps, coeffs} := io_convert_to_internal(polynomials);
        
        bexps . bcoeffs := groebner_groebner(ring, exps, coeffs);
        
        return io_convert_to_output(ring, bexps, bcoeffs)
    end;

% f4 argument error
asserted procedure f4putin_argumentError();
   rederr "usage: f4(polynomials: List, vars: List, order: Any). For example,

          > f4({x*y + 1, y*z + 1}, {x, y, z}, lex);";

%--------------------------------------------------------------------------------------------------

endmodule; % end of f4putin module

end;