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
%     `f5(polynomials: List)`
% Where `polynomials` is a list of expressions, the ideal generators.
% The function returns a list of expressions, the ideal Groebner basis, each as a Lisp Prefix.

% Julia:
%   The file f4putin.red mirrors the file interface.jl from Julia implementation.
%   *All other files have exactly same names as in Julia implementation.*

% Julia:
%   Function names are formed as follows:
%   the file name + the Julia name of the function
%
% So that, julia function "insert_in_hash_table" from the file "hashtable"
% has a name "hashtable_insert_in_hash_table" in this implementation.

create!-package('(f4putin groebner f4 hashtable basis matrix internaltypes io), nil);

%--------------------------------------------------------------------------------------------------

% Reduce-specific:
sizeofInt32 = 32;

asserted procedure cdrn(x, n);
    if n = 0 then
        x
    else
        cdrn(cdr x, n - 1);    

asserted procedure getvlast(x);
    getv(x, length(x));

%--------------------------------------------------------------------------------------------------

put('f4, 'psopfn, 'f4putin_groebner);

asserted procedure f4putin_groebner(u: List): List;
begin scalar x;
    polynomials := reval pop u;
    ring . internalpolys := io_convert_to_internal(polynomials);
    coeffs . exps := internalpolys;
    basis := groebner_groebner(ring, coeffs, exps);
    return io_convert_to_output(basis)
end;

endmodule; % end of f4putin module

end;