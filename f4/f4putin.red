module f4putin;
% f4putin. Main module. The F4 algorithm for computing Groebner bases.
% Contain the AM entry point.

% DEVELOPER NOTE:
% Anywhere in this package sources, comments which start with the word "Julia"
% explain some of our design decisions during the migration of the package.

% This package is a migrated Julia package Groebner.jl maintained by A. Demin et al.
%   https://github.com/sumiya11/Groebner.jl/

% The f4putin module provides the implementation of the Faugere's F4 algorithm
%   https://doi.org/10.1016/S0022-4049(99)00005-5
%
% The interface contains the operator `f4` with the following signature
%     `f4(polynomials: List, vars: List, sortmode)`
% Where `polynomials` is a list of expressions, the ideal generators.
% The function returns a list of expressions, the ideal Groebner basis, each as a Lisp Prefix.

% Julia:
%   The file f4putin.red mirrors the file interface.jl from the Julia implementation.
%   *All other files have exactly same names as in the Julia implementation.*

% Julia:
%   Function names are formed as follows:
%   the file name + the Julia name of the function
%
%   So that julia function "insert_in_hash_table" from the file "hashtable"
%   has a name "hashtable_insert_in_hash_table" in our implementation and is 
%   defined in the file "hashtable.red".

% Julia: 
%   
%
%
%

create!-package('(f4putin groebner f4 hashtable basis matrix internaltypes io), nil);

%--------------------------------------------------------------------------------------------------

% Julia: Reduce-specific things;
%   As we want to mirror the behavior of Julia with integer division masks as close as possible,
%   we use generic Integers (not machine integers) and operate on the lowest `f4putin_sizeofInt32` bits.
fluid '(f4putin_sizeofInt32!*);
f4putin_sizeofInt32!* := 32;

asserted procedure f4putin_cdrn(x, n);
    if n = 0 then
        x
    else
        cdrn(cdr x, n - 1);    

% Julia: analogue of x[end] -- access to the last element of a vector
asserted procedure f4putin_getvlast(x);
    getv(x, length(x));

% Julia:
asserted procedure f4putin_list2vector(x: List): Vector;
    begin scalar len, v;
        len := length(x);
        v := dv_undef(len);
        for i := 1:len do
            putv(v, i, pop(x));
        return v
    end;

% Julia:
asserted procedure f4putin_vector2list(x: Vector): List;
    begin scalar v;
        for i := 1:upbv(x) do
            push(getv(x, i), v);
        return reversip(v)
    end;

%--------------------------------------------------------------------------------------------------

put('f4, 'psopfn, 'f4putin_groebner);

asserted procedure f4putin_groebner(u: List): List;
    begin scalar polynomials;
        if null u or not (listp u) then
            f4putin_argumentError();
        
        % Extract the list of polynomials
        polynomials := reval pop u;        if null u or not (listp u) then
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