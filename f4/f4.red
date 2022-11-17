module f4;
% f4. Main module. The F4 algorithm for computing Groebner bases.
% Contains the AM entryss point.

% This package is migrated from the Julia package Groebner.jl.
% Groebner.jl is maintained by Gowda S. and Demin A
%   https://github.com/sumiya11/Groebner.jl/

% The f4 module provides the implementation of the Faugere's F4 algorithm
%   https://doi.org/10.1016/S0022-4049(99)00005-5
%
% The interface contains the operator `f4` with the following signature
%     `f4(polynomials: List, vars: List, sortmode: Any)`
% Where 
%   . `polynomials` is a list of expressions, the ideal generators,
%   . `vars` is a list of Reduce kernels, the polynomial ring variables,
%   . `sortmode` is a sorting mode for monomial comparison ('lex, 'degrevlex).
% The `f4` function returns a list of expressions, which is the Groebner basis 
% of the ideal; each list element is a Lisp Prefix.

% The file f4.red mirrors the file interface.jl from the Julia implementation.
% (links may not be stable) https://github.com/sumiya11/Groebner.jl/blob/master/src/interface.jl
% All other files have exactly same names as in the Julia implementation
% with the prefix `f4`.

% In Reduce, function names are formed as follows:
% the Julia file name + the Julia name of the function
%
% So that Julia function "insert_in_hash_table" from the file "hashtable"
% has name "hashtable_insert_in_hash_table" in Reduce and is 
% defined in the file "f4hashtable.red".

% Files in this package:
% . f4.red - this file; contains package AM entry point with a single function `f4`;
% . f4basis.red - Pairset and Basis structs implementations, as well as some related functions;
% . f4coeffs.red - manipulations with basis coefficients (reduction modulo / reconstruction);
% . f4constants.red - some hacks to detect the largest possible modulo;
% . f4correctness.red - checking correctness in rational number groebner computation;
% . f4dv.red - dynamically-resized vector implementation;
% . f4f4.red - low-level f4 implementation, the heart of the package;
% . f4groebner.red - high-level Groebner bases over QQ and Zp; calls to f4 over Zp as its backend;
% . f4hashtable.red - MonomialHashtable struct and related functions implementation;
% . f4io.red - input-output conversions of polynomial representations;
% . f4isgroebner.red - checking that something is a groebner basis;
% . f4lucky.red - lucky prime numbers;
% . f4matrix.red - MacaulayMatrix struct implementation; implements linear-algebra-style reduction;
% . f4modular.red - rational number and crt reconstructions;
% . f4normalform.red - normal form computation;
% . f4poly.red - distributive polynomial implementation;
% . f4sorting.red - exponent vector comparators and sorting methods.
create!-package('(f4 f4groebner f4f4 f4hashtable f4basis f4matrix f4sorting f4io f4poly f4dv f4lucky f4coeffs f4correctness f4modular f4normalform f4isgroebner f4constants), nil);

loadtime load!-package 'vector88;
compiletime load!-package 'vector88;

load!-package 'rltools;

% Assertions should be OFF in production.
load!-package 'assert;
off1 'assert;
off1 'assert_procedures;
off1 'assert_inline_procedures;
off1 'assertinstall;
off1 'evalassert;

% from f4io.red
struct PolyRing checked by vectorp;

% from f4basis.red
struct SPair checked by vectorp;
struct Basis checked by vectorp;
struct Pairset checked by vectorp;

% from f4hashtable.red
struct Hashvalue checked by vectorp;
struct MonomialHashtable checked by vectorp;

% from f4matrix.red
struct MacaulayMatrix checked by vectorp;

% All supported coefficient types in F4
procedure f4_isCoeff(x); sqp(x) or integerp(x);
struct Coeff checked by f4_isCoeff;

% Polynomial monomial exponent vector type
struct ExponentVector checked by vectorp;
% The type of entry of an exponent vector
struct Degree checked by integerp;

struct ExponentIdx checked by integerp;

struct DivisionMask checked by integerp;

% Hash of exponent vector
struct ExponentHash checked by integerp;

% MonomsVector of a zero polynomial is an empty `MonomsVector` object.
% CoeffsVector of a zero polynomial is an empty `CoeffsVector` object.
% Vector of polynomial monomials.
struct MonomsVector checked by vectorp; 
% Vector of polynomial coefficients.
struct CoeffsVector checked by vectorp; 

% Column index of a matrix
struct ColumnIdx checked by integerp;

struct PrimeTracker checked by vectorp;

struct CoeffAccum checked by vectorp;

fluid '(!*backtrace);
fluid '(f4_largest!-small!-prime!* 
        f4_largest!-small!-modulus!*);

% Debug mode in f4; 
% prints debug info during execution of f4
switch f4debug=off;

% Uses modular lifting with f4
switch f4modular=on;

on1 'roundbf;

% Julia: Reduce-specific things;
%   As we want to mirror the behavior of Julia with integer division masks as close as possible,
%   we use generic Integers (not machine integers) and operate on the lowest `f4_sizeofInt32` bits.
fluid '(f4_sizeofInt32!*);
f4_sizeofInt32!* := 32;

% Julia: analogue of x[end] -- access the last element of vector x
asserted procedure f4_getvlast(x);
    getv(x, dv_length(x));

% Julia: convert list `x` to a dynamic vector and return it
asserted procedure f4_list2vector(x: List): Vector;
    begin scalar len, v;
        len := length(x);
        v := dv_undef(len);
        for i := 1:len do
            putv(v, i, pop(x));
        return v
    end;

put('f4, 'psopfn, 'f4_groebner);

% AM entry point
asserted procedure f4_groebner(u: List): List;
    begin scalar polynomials, ring, exps, coeffs, ans, 
                    bexps, bcoeffs, vars, ord, w, saveTorder,
                    varsNum, fsq, varsDen;
        if null u or not (listp u) then
            f4_argumentError();
        
        % Extract the list of polynomials
        polynomials := reval pop u;
        if not (listp polynomials) or not (pop polynomials eq 'list) or null polynomials then
            f4_argumentError();

        % Extract the list of variables and the sort mode
        % to initialize the polynomial ring
        % variables and sort mode are specified in f4 call
        saveTorder := if not null u then <<
            % variables and sort mode are specified in f4 call
            vars := reval pop u;
            if not (listp vars) or not (pop vars eq 'list) then
                f4_argumentError();
            for each var in vars do
            if not sfto_kernelp(var) then
                f4_argumentError();
            ord := pop u;
            poly_initRing({vars, ord})
        >> else if not null cdr poly_getGlobalVars() then <<
            % both variables and sort mode have been specified using torder
            poly_initRing(nil)
        >> else <<
            % sort mode has been specified using torder,
            % variables are taken from the inputBasis
            for each f in polynomials do <<
                fsq := simp f;
                varsNum := union(varsNum, kernels numr fsq);
                varsDen := union(varsDen, kernels denr fsq)
            >>;
            if varsDen then
                lprim {varsDen, "implicitly declared as parameters"};
            vars := lto_setminus(varsNum, varsDen);
            vars := sort(vars, 'ordp);
            poly_initRing({vars})
        >>;
        
        % get the current variables and order 
        vars . ord := poly_getVarsAndOrd();
        {ring, exps, coeffs} := io_convert_to_internal(polynomials, vars, ord);

        % run under error catch to recover the order in case of error 
        if !*f4modular then
            w := errorset({'groebner_groebner2, mkquote ring, mkquote exps, mkquote coeffs}, t, !*backtrace)
        else
            w := errorset({'groebner_groebner1, mkquote ring, mkquote exps, mkquote coeffs}, t, !*backtrace);
        
        if errorp w then
            return nil;
        bexps . bcoeffs := car w;
        ans := 'list . io_convert_to_output(ring, bexps, bcoeffs);

        torder cdr saveTorder;
        
        return ans
    end;

% f4 argument error
asserted procedure f4_argumentError();
   rederr "usage: f4(polynomials: List, vars: List, order: Any). For example,

          > f4({x*y + 1, y*z + 1}, {x, y, z}, lex);";

endmodule; % end of f4 module

end;