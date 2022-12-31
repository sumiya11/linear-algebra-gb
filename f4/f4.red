module f4;
% f4. Main module. The F4 algorithm for computing Groebner bases.
% Contains the AM and SM entry points.

revision('f4, "$Id$");

copyright('f4, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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

% This package is migrated from the Julia package Groebner.jl.
% Groebner.jl is maintained by Gowda S. and Demin A (email: asdemin_2@edu.hse.ru)
%   https://github.com/sumiya11/Groebner.jl/
% DOI: TBD

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

% Currently, there is one switch available, it is described below
% . f4modular (default is ON)

% Comments that explain the specific path taken by developers in adapting Julia 
% code to Reduce are usually preceded by the word "Julia:"

% Julia: In Reduce, function names are formed as follows:
% the Julia file name + the Julia name of the function.
% So that Julia function "insert_in_hash_table" from the file "hashtable"
% has name "hashtable_insert_in_hash_table" in Reduce and is 
% defined in the file "f4hashtable.red".

% Source code files in this package:
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

% For neat square bracket array indexing
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

% interface implemented in f4poly.red
procedure f4_isPolynomial(x); eqcar(x, 'p);
procedure f4_isCoeff(x); integerp(x);
struct Polynomial checked by f4_isPolynomial;
struct Terms checked by listp;
struct Term checked by listp;
struct Coeffs checked by listp;
% Coeff is an Integer.
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

% Needed for compatibility with torder from the dipoly package.
fluid '(compiled!-orders!* dipsortmode!* dipsortevcomp!*);
fluid '(dipvars!* global!-dipvars!* vdplastvar!* vdpsortmode!*);

fluid '(!*backtrace);

% These are described and set in f4constants.red
fluid '(f4_largest!-small!-prime!* 
        f4_largest!-small!-modulus!*);

% f4debug.
% Debug mode in f4; 
% Prints debug info during execution of f4.
% Is OFF by default and is not documented (used only for development).
switch f4debug=off;

% f4modular.
% Uses modular lifting with f4.
% In F4 algorithm, it is usually crucial to have this ON,
% since numbers in the matrix reduction get very large.
% Is ON by default. 
switch f4modular=on;

on1 'roundbf;

% Julia: Reduce-specific thing:
% As we want to mirror the behavior of Julia with integer division masks as close as possible,
% we use generic Integers (not machine integers) and operate on the lowest `f4_sizeofInt32` bits.
fluid '(f4_sizeofInt32!*);
f4_sizeofInt32!* := 32;

% Analogue of x[end] in Julia -- access the last element of vector x
asserted procedure f4_getvlast(x);
    getv(x, dv_length(x));

% Convert list `x` to a DynamicVector and return it
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
        >> else if not null cdr global!-dipvars!* then <<
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
            % if there are parameter variables
            if varsDen then
                rederr {"Variables in coefficient denominators are not supported, sorry", varsDen};
            % at this point, varsDen is nil
            vars := varsNum;
            vars := sort(vars, 'ordp);
            poly_initRing({vars})
        >>;
        
        % run under error catch to recover the term order in case of error 
        w := errorset({'f4_groebnerp, mkquote polynomials}, t, !*backtrace);
        torder cdr saveTorder;
        if errorp w then
            return nil;
        return 'list . car w
    end;

% Lisp Prefix --> Lisp Prefix
% Given a list of lisp prefix forms, converts them to polynomials,
% computes the Groebner basis of the corresponding ideal,
% and returns the basis elements as a list of lisp standard prefix.
% Important invariant: 
%    the term order is preserved during the out conversion.   
asserted procedure f4_groebnerp(inputBasisLp: List): List;
    begin scalar polys, gbpolys;
        polys := io_convert_lp_to_poly(inputBasisLp);
        gbpolys := f4_groebnerpoly(polys);
        return io_convert_poly_to_lp(gbpolys)
    end;

% Standard Form --> Standard Form
% Given a list of standard forms, converts them to SQ, then uses `f4_groebnerq`
% to compute the Groebner basis of the corresponding ideal,
% and returns the basis elements as a list of standard forms.
asserted procedure f4_groebnerf(inputBasisSf: List): List;
    begin scalar basis;
        basis := f4_groebnerq(for each f in inputBasisSf collect (f ./ 1));
        return (for each f in basis collect numr f)
    end;

% Standard Quotient --> Standard Quotient
% Given a list of standard quotients, converts them to polynmomials,
% computes the Groebner basis of the corresponding ideal,
% and returns the basis elements as a list of standard quotients.
asserted procedure f4_groebnerq(inputBasisSq: List): List;
    begin scalar basis;
        basis := f4_groebnerpoly(for each f in inputBasisSq collect 
                        poly_sq2poly f);
        return (for each f in basis collect poly_poly2sq f)
    end;

% Polynomial --> Polynomial
% Computes a Groebner basis of the ideal generated 
% by the given polynomials.
asserted procedure f4_groebnerpoly(inputBasis: List): List;
    begin scalar properIdeal, w, p, ring, exps, coeffs, gbexps, gbcoeffs;
        w := inputBasis; inputBasis := nil;
        properIdeal := t; while properIdeal and w do <<
            p := pop w;
            % get rid of zeros in input 
            if not poly_iszero!?(p) then <<
                push(p, inputBasis); 
                % if constant present in the input,
                % the ideal is improper
                if poly_isConst!?(p) then
                properIdeal := nil
            >>
        >>;
        if not properIdeal then
            return {poly_one()};
        if null inputBasis then
            % This is a bit unclear mathematically, 
            % but we go with the design decisions of the GROEBNER package
            return {poly_zero()};
        {ring, exps, coeffs} := io_convert_poly_to_internal(inputBasis);
        if !*f4modular then
            gbexps . gbcoeffs := groebner_groebner2(ring, exps, coeffs)
        else
            gbexps . gbcoeffs := groebner_groebner1(ring, exps, coeffs);
        return io_convert_internal_to_poly(ring, gbexps, gbcoeffs)
    end;

% f4 argument error
asserted procedure f4_argumentError();
   rederr "usage: f4(polynomials: List, vars: List, order: Any). For example,

          > f4({x*y + 1, y*z + 1}, {x, y, z}, lex);

          Or, using torder:

          > torder({x, y, z}, lex);
          > f4({x*y + 1, y*z + 1});
          ";

endmodule; % end of f4 module

end;
