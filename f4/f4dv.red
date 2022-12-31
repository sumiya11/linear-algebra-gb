module f4dv;
% Dynamic Vector. Implements interface of dynamic vector, which allows resizing.
% This has no analogue in Groebner.jl.

revision('f4dv, "$Id$");

copyright('f4dv, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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

struct DynamicVector checked by vectorp;

% In julia array indexing is 1-based vs. Reduce 0-based.
% We emulate 1-based indexing by creating arrays of size 1 larger than needed.

% Vector of size n
asserted procedure dv_undef(n: Integer): DynamicVector;
    mkvect(n);

% Length of vector x
asserted procedure dv_length(v: DynamicVector): Integer;
    upbv(v);

asserted procedure dv_isempty(v: DynamicVector): Boolean;
    dv_length(v) = 0;

% Vector of size n filled with zeros
asserted procedure dv_zeros(n): DynamicVector;
    begin scalar v;
        v := dv_undef(n);
        for i := 1:n do
            putv(v, i, 0);
        return v
    end;

% Vector of size n filled with zeros
asserted procedure dv_zeros_sq(n): DynamicVector;
    begin scalar v;
        v := dv_undef(n);
        for i := 1:n do
            putv(v, i, nil ./ 1);
        return v
    end;

% Vector of size equal to size of `x`
asserted procedure dv_similar(x): DynamicVector;
    dv_undef(dv_length(x));

% Deepcopy of x.
% Elements are *not* copied recursively.
asserted procedure dv_copy(x): DynamicVector;
    begin scalar v;
        v := dv_similar(x);
        for i := 1:dv_length(x) do
            putv(v, i, getv(x, i));
        return v
    end;

% Resize vector x to be of length n;
% 
% . If length(x) == n, nothing happens
% . If length(x) < n, new array of size n is reallocated 
% and elements from x are copied to the beginning of new array
% . If length(x) > n, new array of size n is reallocated
% and the first n elements from x are copied to the beginning of new array
asserted procedure dv_resize(x, n): DynamicVector;
    begin scalar v;
        if dv_length(x) = n then
            return x;
        v := dv_undef(n);
        for i := 1:min(n, dv_length(x)) do
            putv(v, i, getv(x, i));
        return v
    end;

% Same as `dv_resize`, but fills new elements with zeros
asserted procedure dv_resizezeros(x, n): DynamicVector;
    begin scalar v;
        v := dv_resize(x, n);
        if n > dv_length(x) then
            for i := dv_length(x)+1:n do
                putv(v, i, 0);
        return v 
    end;

endmodule; % end of dv module

end;