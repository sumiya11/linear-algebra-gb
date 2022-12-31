module f4sorting;
% Sorting routines used in f4.
% This file corresponds to file f4/sorting.jl in Groebner.jl 

revision('f4sorting, "$Id$");

copyright('f4sorting, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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

% Here, we want to use fast term orders on vectors for simple orders like
% lex or revgradlex, and allow all other orders available through torder.
% Thus, `sorting_exponent_isless` should be used as a generic comparator. 
% It calls either our implementation, or the implementation in the dipoly package.
% `sorting_exponent_isless` dispatches on the value of vdpsortmode!*.

% Degrevlex exponent vector comparator
asserted procedure sorting_exponent_isless_drl(ea: ExponentVector, eb: ExponentVector): Boolean;
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

% Deglex exponent vector comparator
asserted procedure sorting_exponent_isless_dl(ea: ExponentVector, eb: ExponentVector): Boolean;
    begin
        if f4_getvlast(ea) < f4_getvlast(eb) then
            return t
        else if not (f4_getvlast(ea) = f4_getvlast(eb)) then
            return nil;
        return sorting_exponent_isless_lex(cdr ea, cdr eb)
    end;

% Lex exponent vector comparator
asserted procedure sorting_exponent_isless_lex(ea: ExponentVector, eb: ExponentVector): Boolean;
    begin integer i;
        i := 1;
        while (i < dv_length(ea) - 1) and getv(ea, i) = getv(eb, i) do
            i := i + 1;
        
        if getv(ea, i) <= getv(eb, i) then
            return t
        else
            return nil
    end;

% Generic exponent vector comparator
asserted procedure sorting_exponent_isless(e1: ExponentVector, e2: ExponentVector): Boolean;
    begin scalar el1, el2, i, n;
        return if vdpsortmode!* eq 'lex then
            sorting_exponent_isless_lex(e1, e2)
        else if vdpsortmode!* eq 'gradlex then
            sorting_exponent_isless_dl(e1, e2)
        else if vdpsortmode!* eq 'revgradlex then
            sorting_exponent_isless_drl(e1, e2)
        else <<
            n := dv_length(e1);
            el1 := for i := 1:n-1 collect getv(e1, i);
            el2 := for i := 1:n-1 collect getv(e2, i);
            evcomp(el1, el2) = -1
        >>
   end;

asserted procedure sorting_comparator_sort_gens_by_lead_increasing(x, y);
    sorting_exponent_isless(car x, car y);

% Sorts terms and corresponding coefficients from the basis
% by their leading term in non-decreasing order.
%
% Used only once to sort input generators
asserted procedure sorting_sort_gens_by_lead_increasing(basis: Basis, ht: MonomialHashtable);
    begin scalar gens, exps, inds, coeffs, x;
        gens := basis_bget_gens(basis);
        coeffs := basis_bget_coeffs(basis);
        exps := hashtable_htget_exponents(ht);

        inds := for i := 1:basis_bget_ntotal(basis) collect 
            {getv(exps, getv(getv(gens, i), 1)), getv(gens, i), getv(coeffs, i)};

        inds := sort(inds, 'sorting_comparator_sort_gens_by_lead_increasing);

        for i := 1:basis_bget_ntotal(basis) do <<
            x := pop(inds);
            putv(gens, i, cadr x);
            putv(coeffs, i, caddr x)
        >>
    end;

% -//-
asserted procedure sorting_sort_gens_by_lead_increasing2(basis: Basis, ht: MonomialHashtable, coeffs_zz: Vector);
    begin scalar gens, exps, inds, coeffs, x;
        gens := basis_bget_gens(basis);
        coeffs := basis_bget_coeffs(basis);
        exps := hashtable_htget_exponents(ht);

        inds := for i := 1:basis_bget_ntotal(basis) collect 
            {getv(exps, getv(getv(gens, i), 1)), getv(gens, i), getv(coeffs, i), getv(coeffs_zz, i)};

        inds := sort(inds, 'sorting_comparator_sort_gens_by_lead_increasing);

        for i := 1:basis_bget_ntotal(basis) do <<
            x := pop(inds);
            putv(gens, i, cadr x);
            putv(coeffs, i, caddr x);
            putv(coeffs_zz, i, cadddr x)
        >>
    end;

asserted procedure sorting_comparator_sort_input_to_change_ordering(x, y);
    not sorting_exponent_isless(car x, car y);

% Sort polynomials to change the order of terms if needed
asserted procedure sorting_sort_input_to_change_ordering(exps, coeffs, ord);
    begin scalar inds, coeffs, x;

        inds := for i := 1:dv_length(exps) collect 
            {getv(exps, i), getv(coeffs, i)};

        inds := sort(inds, 'sorting_comparator_sort_input_to_change_ordering);

        for i := 1:dv_length(exps) do <<
            x := pop(inds);
            putv(exps, i, car x);
            putv(coeffs, i, cadr x)
        >>
    end;

asserted procedure sorting_comparator_sort_gens_by_lead_increasing_in_standardize(x, y);
    sorting_exponent_isless(car x, car y);

% Defined the order of polynomials in the output of f4
asserted procedure sorting_sort_gens_by_lead_increasing_in_standardize(basis: Basis, ht: MonomialHashtable, ord);
    begin scalar gens, exps, nnr, lead, coeffs, inds, x;
        gens := basis_bget_gens(basis);
        exps := hashtable_htget_exponents(ht);
        lead := basis_bget_lead(basis);
        nnr := basis_bget_nonred(basis);
        coeffs := basis_bget_coeffs(basis);

        inds := for i := 1:basis_bget_nlead(basis) collect 
            {getv(exps, getv(getv(gens, i), 1)), getv(gens, i), getv(coeffs, i), getv(lead, i)};

        inds := sort(inds, 'sorting_comparator_sort_gens_by_lead_increasing_in_standardize);
        
        for i := 1:basis_bget_nlead(basis) do <<
            x := pop(inds);
            putv(gens, i, cadr x);
            putv(coeffs, i, caddr x);
            putv(lead, i, cadddr x)
        >>
    end;

% returns degree(x) < degree(y), where x and y are critical pairs
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

% For sorting the rows in the upper part of the matrix.
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
            return t;
        if va < vb then
            return nil;

        return nil
    end;

% For sorting the rows in the lower part of the matrix.
asserted procedure sorting_comparator_sort_matrix_rows_increasing(a, b);
    begin scalar va, vb;

        a := cadr a;
        b := cadr b;

        va := getv(a, 1);
        vb := getv(b, 1);

        if va > vb then
            return t;
        if va < vb then
            return nil;

        % same column index => compare density of rows

        va := dv_length(a);
        vb := dv_length(b);

        if va > vb then
            return nil;
        if va < vb then
            return t;

        return nil
    end;

% Sort matrix rows by column index and density
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
            x := pop(inds);
            putv(uprows, i, cadr x);
            putv(up2coef, i, caddr x)
        >>
    end;

% Sort matrix rows by column index and density
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

asserted procedure sorting_comparator_sort_pairset_by_lcm(x, y);
    sorting_exponent_isless(cadr x, cadr y);

% sorts first npairs pairs from pairset by increasing of
% lcm exponent wrt the given monomial ordering
%
% Used only in normal selection strategy once per f4 iteration
asserted procedure sorting_sort_pairset_by_lcm(pairset: Pairset, npairs: Integer, ht: MonomialHashtable);
    begin scalar part, ps, exps;
        ps := basis_psget_pairs(pairset);
        exps := hashtable_htget_exponents(ht);

        part := for i := 1:npairs collect {getv(ps, i), getv(exps, basis_spget_lcmj(getv(ps, i)))};
        part := sort(part, 'sorting_comparator_sort_pairset_by_lcm);

        for i := 1:npairs do
            putv(ps, i, car pop(part))
    end;

% sorts generators selected in normal strategy function
% by their ordering in the current basis (identity sort)
asserted procedure sorting_sort_generators_by_position(gens: Vector, loadj: Integer);
    begin scalar part;
        part := for i := 1:loadj collect getv(gens, i);

        part := sort(part, '<);

        for i := 1:loadj do
            putv(gens, i, pop(part))
    end;

asserted procedure sorting_comparator_sorting_sort_columns_by_hash(a, b);
    begin scalar ha, hb, ea, eb;
        ha := cadr a;
        hb := cadr b;
        if not (hashtable_hvget_idx(ha) = hashtable_hvget_idx(hb)) then
            return hashtable_hvget_idx(ha) > hashtable_hvget_idx(hb);
        ea := caddr a;
        eb := caddr b;
        return not sorting_exponent_isless(ea, eb)
    end;

asserted procedure sorting_sort_columns_by_hash(col2hash: Vector, symbol_ht: MonomialHashtable);
    begin scalar hd, es, len, listtosort;
        hd := hashtable_htget_hashdata(symbol_ht);
        es := hashtable_htget_exponents(symbol_ht);
        len := hashtable_htget_explen(symbol_ht);

        listtosort := for i := 1:dv_length(col2hash) collect 
            {getv(col2hash, i), getv(hd, getv(col2hash, i)), getv(es, getv(col2hash, i))};
        
        listtosort := sort(listtosort, 'sorting_comparator_sorting_sort_columns_by_hash);

        for i := 1:dv_length(col2hash) do
            putv(col2hash, i, car pop(listtosort))
    end;

endmodule; % end of sorting module

end; % end of file
