module f4internaltypes;
% 

%--------------------------------------------------------------------------------------------------

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

%--------------------------------------------------------------------------------------------------

endmodule; % end of internaltypes module

end;