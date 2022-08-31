module f4;
% f4. Main module. The F4 algorithm for computing Groebner bases.
% Contains the AM entry point.

% This package is migrated from Julia package Groebner.jl.
% Groebner.jl is maintained by Gowda S. and Demin A. 
%   https://github.com/sumiya11/Groebner.jl/

% DEVELOPER NOTE:
% Anywhere in this package sources, comments that start with the word "Julia"
% explain some of our design decisions during the migration of the package to Reduce.

% DEVELOPER NOTE [examples]:
% f4();
%

% The f4 module provides the implementation of the Faugere's F4 algorithm
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
%   The file f4.red mirrors the file interface.jl from the Julia implementation.
%       (links are not stable)
%       https://github.com/sumiya11/Groebner.jl/blob/master/src/interface.jl
%   All other files have exactly same names as in the Julia implementation.

% Julia:
%   In Reduce, function names are formed as follows:
%   the file name + the Julia name of the function
%
%   So that Julia function "insert_in_hash_table" from the file "hashtable"
%   has name "hashtable_insert_in_hash_table" in Reduce and is 
%   defined in the file "f4hashtable.red".

%--------------------------------------------------------------------------------------------------

% . f4 - this file, contains package AM entry point with a single function `f4`;
% . groebner - high-level Groebner basis computation; calls to f4 as its backend;
% . f4 - low-level f4 implementation, the heart of the package;
% . hashtable - MonomialHashtable struct and related functions implementation;
% . basis - Pairset and Basis structs implementations, as well as some related functions;
% . matrix - MacaulayMatrix struct implementation; contains linear-algebra-style reduction;
% . internaltypes - type declarations (have no practical effect in Reduce, 
%                   but are left for consistency with Julia implementation)
% . io - input-output conversions of polynomial representations
create!-package('(f4 f4groebner f4f4 f4hashtable f4basis f4matrix f4sorting f4io f4poly f4dv f4lucky f4coeffs f4correctness f4modular), nil);

loadtime load!-package 'vector88;
compiletime load!-package 'vector88;

load!-package 'rltools;

% Assertions should be OFF in production.
load!-package 'assert;
on1 'assert;
on1 'assert_procedures;
on1 'assert_inline_procedures;
on1 'assertinstall;
on1 'evalassert;

% from f5io.red
struct PolyRing;

% from basis.red
struct SPair;
struct Basis;
struct Pairset;

% from hashtable.red
struct Hashvalue;
struct MonomialHashtable;


% from matrix.red
struct MacaulayMatrix;

% All supported coefficient types in F4
struct Coeff; % checked by listp;

% Polynomial monomial exponent vector type
struct ExponentVector; % = Vector{UInt16}
% The type of entry of an exponent vector
struct Degree; % = eltype(ExponentVector)

struct ExponentIdx; % = Int32


struct DivisionMask; % = UInt32

% Hash of exponent vector
struct ExponentHash; %  = UInt32

% MonomsVector of a zero polynomial is an empty `MonomsVector` object.
% CoeffsVector of a zero polynomial is an empty `CoeffsVector` object.

% Vector of polynomial monomials
struct MonomsVector; % = Vector{ExponentIdx}
% Vector of polynomial coefficients
struct CoeffsVector; % = Vector{Coeff}

% Column index of a matrix
struct ColumnIdx; % = Int32

struct PrimeTracker;

struct CoeffAccum;

fluid '(coeff_plus!* coeff_times!*);

procedure coeff_modular(flag);
    if flag then <<

    >> else <<
        coeff_plus!* := function addsq;
        coeff_times!* := function multsq
    >>;

inline procedure coeff_plus(x, y);
    apply(coeff_plus!*, {x, y});

fluid '(!*backtrace);

#if (errorp (errorset (quote (itimes2 (expt 2 55) 2)) nil nil))
   procedure wuwu();
      "2^56 is bad";
#else
   procedure wuwu();
      "2^56 is good";
#endif

% turn to switch
asserted procedure f4_debug();
    nil;

% Julia: Reduce-specific things;
%   As we want to mirror the behavior of Julia with integer division masks as close as possible,
%   we use generic Integers (not machine integers) and operate on the lowest `f4_sizeofInt32` bits.
fluid '(f4_sizeofInt32!*);
f4_sizeofInt32!* := 32;

% Julia: analogue of x[end] -- access the last element of vector x
asserted procedure f4_getvlast(x);
    getv(x, dv_length(x));

% Julia: convert list `x` to vector and return it
asserted procedure f4_list2vector(x: List): Vector;
    begin scalar len, v;
        len := length(x);
        v := dv_undef(len);
        for i := 1:len do
            putv(v, i, pop(x));
        return v
    end;

%--------------------------------------------------------------------------------------------------

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
        w := errorset({'groebner_groebner2, mkquote ring, mkquote exps, mkquote coeffs}, t, !*backtrace);
        
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

%--------------------------------------------------------------------------------------------------

endmodule; % end of f4 module

end;

