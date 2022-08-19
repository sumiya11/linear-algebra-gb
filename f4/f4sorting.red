module f4sorting;
% Sorting routines used in f4.

%--------------------------------------------------------------------------------------------------

% degrevlex exponent vector comparison
asserted procedure sorting_exponent_isless_drl(ea: ExponentVector, eb: ExponentVector): Bool;
    begin integer i;
        if f4_getvlast(ea) < f4_getvlast(eb) then
            return t
        else if not (f4_getvlast(ea) = f4_getvlast(eb)) then
            return nil;
        
        i := dv_length(ea) - 1;
        while i > 1 and getv(ea, i) = getv(eb, i) do
            i := i - 1;
        
        if getv(ea, i) <= getv(eb, i) then
            return nil
        else
            return t
    end;

% lex exponent vector comparison
asserted procedure sorting_exponent_isless_lex(ea: ExponentVector, eb: ExponentVector): Bool;
    begin integer i;
        i := 1;
        while (i < dv_length(ea) - 1) and getv(ea, i) = getv(eb, i) do
            i := i + 1;
        
        if getv(ea, i) <= getv(eb, i) then
            return t
        else
            return nil
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_comparator_sort_gens_by_lead_increasing(x, y);
    sorting_exponent_isless_drl(car x, car y);

% sorts generators and corresponding coefficients from basis
% by their leading monomial in increasing ordering
%
% Used only once to sort input generators
asserted procedure sorting_sort_gens_by_lead_increasing(basis: Basis, ht: MonomialHashtable);
    begin scalar gens, exps, inds, coeffs, x;
        gens := basis_bget_gens(basis);
        coeffs := basis_bget_coeffs(basis);
        exps := hashtable_htget_exponents(ht);

        if f4_debug() then <<
            prin2t {"sort_gens_by_lead_increasing: input"};
            prin2t gens;
            prin2t exps;
            prin2t basis_bget_ntotal(basis)
        >>;

        inds := for i := 1:basis_bget_ntotal(basis) collect 
            {getv(exps, getv(getv(gens, i), 1)), getv(gens, i), getv(coeffs, i)};

        if f4_debug() then
            prin2t {"sort_gens_by_lead_increasing: before", inds};

        inds := sort(inds, 'sorting_comparator_sort_gens_by_lead_increasing);
        
        if f4_debug() then
            prin2t {"sort_gens_by_lead_increasing: after", inds};

        for i := 1:basis_bget_ntotal(basis) do <<
            x . inds := inds;
            putv(gens, i, cadr x);
            putv(coeffs, i, caddr x)
        >>
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_comparator_sort_gens_by_lead_increasing_in_standardize(x, y);
    sorting_exponent_isless_drl(car x, car y);

asserted procedure sorting_sort_gens_by_lead_increasing_in_standardize(basis: Basis, ht: MonomialHashtable, ord);
    begin scalar gens, exps, nnr, lead, genscopy, coeffs, coeffscopy, leadcopy, inds2;
        gens := basis_bget_gens(basis);
        exps := hashtable_htget_exponents(ht);
        lead := basis_bget_lead(basis);
        nnr := basis_bget_nonred(basis);
        coeffs := basis_bget_coeffs(basis);

        inds := for i := 1:basis_bget_nlead(basis) collect 
            {getv(exps, getv(getv(gens, i), 1)), getv(gens, i), getv(coeffs, i), getv(lead, i)};

        inds := sort(inds, 'sorting_comparator_sort_gens_by_lead_increasing_in_standardize);
        
        for i := 1:basis_bget_nlead(basis) do <<
            x . inds := inds;
            putv(gens, i, cadr x);
            putv(coeffs, i, caddr x);
            putv(lead, i, cadddr x)
        >>
    end;

%--------------------------------------------------------------------------------------------------

% returns degree(x) < degree(y)
asserted procedure sorting_comparator_sort_pairset_by_degree(x: SPair, y: SPair);
    basis_spget_deg(x) < basis_spget_deg(y);

% Sorts pairs from pairset in range [from, from+size]
% by lcm total degrees in increasing order
%
% Used in update once per one f4 iteration to sort pairs in pairset
% Also used in normal selection strategy also once per iteration
asserted procedure sorting_sort_pairset_by_degree(ps: Pairset, from: Integer, sz: Integer);
    begin scalar pairs, part;
        pairs := basis_psget_pairs(ps);
        part := for i := from:from+sz collect getv(pairs, i);

        part := sort(part, 'sorting_comparator_sort_pairset_by_degree);

        for i := from:from+sz do
            putv(pairs, i, pop(part))
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_comparator_sort_matrix_rows_decreasing(a, b);
    begin scalar va, vb;

        a := cadr a;
        b := cadr b;

        va := getv(a, 1);
        vb := getv(b, 1);

        if va > vb then
            return nil;
        if va < vb then
            return t;

        % same column index => compare density of rows

        va := dv_length(a);
        vb := dv_length(b);

        if va > vb then
            return nil;
        if va < vb then
            return t;

        return nil
    end;

asserted procedure sorting_comparator_sort_matrix_rows_increasing(a, b);
    begin scalar va, vb;

        a := cadr a;
        b := cadr b;

        va := getv(a, 1);
        vb := getv(b, 1);

        if va > vb then
            return nil;
        if va < vb then
            return t;

        % same column index => compare density of rows

        va := dv_length(a);
        vb := dv_length(b);

        if va > vb then
            return nil;
        if va < vb then
            return t;

        return nil
    end;

% by column index and density
asserted procedure sorting_sort_matrix_rows_decreasing(matrixj: MacaulayMatrix);
    begin scalar inds, uprows, up2coef, x;
        % smaller means  pivot being more left  %
        % and density being smaller             %
        uprows := matrix_mget_uprows(matrixj);
        up2coef := matrix_mget_up2coef(matrixj);

        inds := for i := 1 : matrix_mget_nup(matrixj) collect 
            {i, getv(uprows, i), getv(up2coef, i)};

        inds := sort(inds, 'sorting_comparator_sort_matrix_rows_decreasing);
        
        for i := 1:matrix_mget_nup(matrixj) do <<
            x . inds := inds;
            putv(uprows, i, cadr x);
            putv(up2coef, i, caddr x)
        >>
    end;

% by column index and density
asserted procedure sorting_sort_matrix_rows_increasing(matrixj: MacaulayMatrix);
    begin scalar inds, lowrows, low2coef, x;
        % smaller means  pivot being more left  %
        % and density being smaller             %
        lowrows := matrix_mget_lowrows(matrixj);
        low2coef := matrix_mget_low2coef(matrixj);

        inds := for i := 1 : matrix_mget_nlow(matrixj) collect 
            {i, getv(lowrows, i), getv(low2coef, i)};

        inds := sort(inds, 'sorting_comparator_sort_matrix_rows_increasing);
        
        for i := 1:matrix_mget_nlow(matrixj) do <<
            x . inds := inds;
            putv(lowrows, i, cadr x);
            putv(low2coef, i, caddr x)
        >>
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_comparator_sort_pairset_by_lcm_drl(x, y);
    sorting_exponent_isless_drl(cadr x, cadr y);

% sorts first npairs pairs from pairset by increasing of
% lcm exponent wrt the given monomial ordering
%
% Used only in normal selection strategy once per f4 iteration
asserted procedure sorting_sort_pairset_by_lcm(pairset: Pairset, npairs: Integer, ht: MonomialHashtable);
    begin scalar part, ps, exps;
        ps := basis_psget_pairs(pairset);
        exps := hashtable_htget_exponents(ht);

        if f4_debug() then
            prin2t {"sort_pairset_by_lcm:", ps, ", npairs:", npairs};

        part := for i := 1:npairs collect {getv(ps, i), getv(exps, basis_spget_lcmj(getv(ps, i)))};

        if f4_debug() then
            prin2t {"sort_pairset_by_lcm: before", part};

        part := sort(part, 'sorting_comparator_sort_pairset_by_lcm_drl);
        
        if f4_debug() then
            prin2t {"sort_pairset_by_lcm: after", part};

        for i := 1:npairs do
            putv(ps, i, car pop(part))
    end;

%--------------------------------------------------------------------------------------------------

% sorts generators selected in normal strategy function
% by their ordering in the current basis (identity sort)
asserted procedure sorting_sort_generators_by_position(gens: Vector, loadj: Int);
    begin scalar part;
        part := for i := 1:loadj collect getv(gens, i);

        if f4_debug() then
            prin2t {"sort_generators_by_position: before sort: ", part, loadj};

        part := sort(part, '<);

        if f4_debug() then
            prin2t {"sort_generators_by_position: after sort: ", part, loadj};

        for i := 1:loadj do
            putv(gens, i, pop(part))
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_comparator_sorting_sort_columns_by_hash_drl(a, b);
    begin scalar ha, hb, ea, eb;
        ha := cadr a;
        hb := cadr b;
        if not (hashtable_hvget_idx(ha) = hashtable_hvget_idx(hb)) then
            return hashtable_hvget_idx(ha) > hashtable_hvget_idx(hb);
        ea := caddr a;
        eb := caddr b;
        return not sorting_exponent_isless_drl(ea, eb)
    end;

asserted procedure sorting_sort_columns_by_hash(col2hash: Vector, symbol_ht: MonomialHashtable);
    begin scalar hd, es, len, col2hash_list;
        hd := hashtable_htget_hashdata(symbol_ht);
        es := hashtable_htget_exponents(symbol_ht);
        len := hashtable_htget_explen(symbol_ht);

        listtosort := for i := 1:dv_length(col2hash) collect 
            {getv(col2hash, i), getv(hd, getv(col2hash, i)), getv(es, getv(col2hash, i))};
        
        listtosort := sort(listtosort, 'sorting_comparator_sorting_sort_columns_by_hash_drl);

        for i := 1:dv_length(col2hash) do
            putv(col2hash, i, car pop(listtosort))
    end;

%--------------------------------------------------------------------------------------------------

endmodule; % end of sorting module

end; % end of file