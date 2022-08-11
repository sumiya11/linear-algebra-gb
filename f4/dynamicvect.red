
struct DynamicVector;

asserted procedure dv_undef(n: Integer): DynamicVector;
    mkvect(n);

asserted procedure dv_zeros(n): DynamicVector;
    begin scalar v;
        v := dv_undef(n);
        for i := 1:n do
            putv(v, i, 0);
        return v
    end;

asserted procedure dv_similar(x): DynamicVector;
    dv_undef(length(x));

asserted procedure dv_resize(x, n): DynamicVector;
    begin scalar v;
        ASSERT(length(x) < n);
        v := dv_undef(n);
        for i := 1:length(x) do
            putv(v, i, getv(x, i));
        return v
    end;

asserted procedure dv_resizezeros(x, n): DynamicVector;
    begin scalar v;
        start := length(x) + 1;
        v := resize(x, n);
        for i := start:length(v) do
            putv(v, i, 0);
        return v
    end;

end;