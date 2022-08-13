module dv;
% Dynamic Vector.
% Implements DynamicVector interface.

struct DynamicVector;

% Julia:
%   In julia array indexing is 1-based vs. Reduce 0-based.
%   We emulate 1-based indexing by creating arrays of size 1 larger than needed.

% Vector of size n
asserted procedure dv_undef(n: Integer): DynamicVector;
    mkvect(n);

% Vector of size n filled with zeros
asserted procedure dv_zeros(n): DynamicVector;
    begin scalar v;
        v := dv_undef(n);
        for i := 1:n do
            putv(v, i, 0);
        return v
    end;

% Vector of size equal to size of `x`
asserted procedure dv_similar(x): DynamicVector;
    dv_undef(length(x));

% Resize vector x to be of length n;
% 
% . If length(x) == n, nothing happens
% . If length(x) < n, new array of size n is reallocated 
% and elements from x are copied to the beginning of new array
% . If length(x) > n, new array of size n is reallocated
% and the first n elements from x are copied to the beginning of new array
asserted procedure dv_resize(x, n): DynamicVector;
    begin scalar v;
        if length(x) = n then
            return x;
        v := dv_undef(n);
        for i := 1:min(n, length(x)) do
            putv(v, i, getv(x, i));
        return v
    end;

% Same as `dv_resize`, but fills new elements with zeros
asserted procedure dv_resizezeros(x, n): DynamicVector;
    begin scalar v;
        v := dv_resize(x, n);
        if n > length(x) then
            for i := length(x)+1:n do
                putv(v, i, 0);
        return v 
    end;

end;