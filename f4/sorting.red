module sorting;

% sorts generators and corresponding coefficients from basis
% by their leading monomial in increasing ordering
%
% Used only once to sort input generators
asserted procedure sorting_sort_gens_by_lead_increasing(basis: Basis, ht: MonomialHashtable);
    begin scalar gens;

    end;

procedure comparator(x, y);
    car x < car y;

procedure sort_two_at_the_same_time(list_1: List, list_2: List, values: Vector);
    begin scalar n, list_with_values;
        n := length(list_1);
        for i := 1:n do
            push({getv(values, i), pop(list_1), pop(list_2)}, list_with_values);
        prin2t list_with_values;
        list_with_values := sort(list_with_values, 'comparator);
        list_1 := for each x in list_with_values collect cadr x;
        list_2 := for each x in list_with_values collect caddr x;
        return list_1 . list_2
    end;

%--------------------------------------------------------------------------------------------------

% degrevlex exponent vector comparison
asserted procedure sorting_exponent_isless_drl(ea: ExponentVector, eb: ExponentVector): Bool;
    begin integer i;
        if f4putin_getvlast(ea) < f4putin_getvlast(eb) then
            return t
        else if not (f4putin_getvlast(ea) = f4putin_getvlast(eb)) then
            return nil;
        
        i := length(ea) - 1;
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
        while (i < length(ea) - 1) and getv(ea, i) = getv(eb, i) do
            i := i + 1;
        
        if getv(ea, i) <= getv(eb, i) then
            return t
        else
            return nil
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_sort_gens_by_lead_increasing_in_standardize(basis: Basis. ht: MonomialHashtable);
    begin scalar gens;

    end;

%--------------------------------------------------------------------------------------------------

% returns degree(x) < degree(y)
asserted procedure sorting_internal_pairset_degree(x: SPair, y: SPair);
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

        sort(part, 'sorting_internal_pairset_degree);

        for i := from:from+sz do
            putv(pairs, i, pop(part))
    end;

%--------------------------------------------------------------------------------------------------

%  by column index and density
asserted procedure sorting_sort_matrix_rows_decreasing(matrix: MacaulayMatrix);
    begin scalar inds;
        % smaller means  pivot being more left  %
        % and density being smaller             %

        
    end;

%  by column index and density
asserted procedure sorting_sort_matrix_rows_decreasing(matrix: MacaulayMatrix);
    begin scalar inds;
        % smaller means  pivot being more right =#
        % and density being larger =#
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_sort_pairset_by_lcm(pairset: Pairset, npairs: Integer, ht: MonomialHashtable);
    begin scalar part;

    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_sort_generators_by_position(gens: Vector, load: Int);
    begin scalar part;
        part := for i := 1:load collect getv(gens, i);

        sort(part, '<);

        for i := 1:load do
            putv(gens, i, pop(part))
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure sorting_sort_columns_by_hash(col2hash: Vector, symbol_ht: MonomialHashtable);
    being scalar x;

    end;

endmodule; % end of sorting module

end; % end of file