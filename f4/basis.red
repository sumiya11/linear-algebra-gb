

% SPair is
%   poly1: Integer
%   poly2: Integer
%   lcm:   ExponentIdx
%   deg:   Degree

% Immutable
asserted procedure basis_SPair(poly1: Integer, poly2: Integer, 
                        lcmj: ExponentIdx, deg: Degree): SPair;
    {'sp, poly1, poly2, lcmj, deg};

asserted procedure basis_spget_poly1(sp: SPair): Integer;
    car cdrn(sp, 1);

asserted procedure basis_spget_poly2(sp: SPair): Integer;
    car cdrn(sp, 2);

asserted procedure basis_spget_lcmj(sp: SPair): ExponentIdx;
    car cdrn(sp, 3);

asserted procedure basis_spget_deg(sp: SPair): Degree;
    car cdrn(sp, 4);

%--------------------------------------------------------------------------------------------------

% Pairset is
%   pairs: Vector{SPair}
%   load:  Integer

% load -> ld
asserted procedure basis_Pairset(pairs: Vector, ld: Integer): Pairset;
    {'ps, pairs, ld};

% Initialize and return a pairset with capacity for `initial_size` pairs
asserted procedure basis_initialize_pairset();
    begin scalar initial_size;
        initial_size := 8;
        pairs := dv_undef(initial_size);
        return basis_Pairset(pairs, 0)
    end;

asserted procedure basis_psget_pairs(ps: Pairset): Vector;
    car cdrn(ps, 1);

asserted procedure basis_psget_load(ps: Pairset): Integer;
    car cdrn(ps, 2);

asserted procedure basis_psset_pairs(ps: Pairset, x);
    car cdrn(ps, 1) := x;

asserted procedure basis_psset_load(ps: Pairset, x);
    car cdrn(ps, 2) := x;

% checks if it's possible to add `added` number of pairs to the pairset,
% and extend the pairset if not
asserted procedure basis_check_enlarge_pairset(ps: Pairset, added: Integer);
    begin scalar sz;
        sz := length(basis_psget_pairs(ps));
        if basis_psget_load(ps) + added > sz then <<
            newsz := max(2*sz, basis_psget_load(ps) + added);
            basis_psset_pairs(ps, dv_resize(basis_psget_pairs(ps), newsz))
        >>;
    end;

%--------------------------------------------------------------------------------------------------

%
% 
%
asserted procedure basis_Basis(gens, coeffs, 
                            size, ndone, ntotal, 
                            isred, nonred, lead, nlead): Basis;
    {'b, gens, coeffs, size, ndone, ntotal, isred, nonred, lead, nlead};

asserted procedure basis_bget_gens(b: Basis);
    car cdrn(b, 1);

asserted procedure basis_bget_coeffs(b: Basis);
    car cdrn(b, 2);

asserted procedure basis_bget_size(b: Basis);
    car cdrn(b, 3);

asserted procedure basis_bget_ndone(b: Basis);
    car cdrn(b, 4);

asserted procedure basis_bget_ntotal(b: Basis);
    car cdrn(b, 5);

asserted procedure basis_bget_isred(b: Basis);
    car cdrn(b, 6);

asserted procedure basis_bget_nonred(b: Basis);
    car cdrn(b, 7);

asserted procedure basis_bget_lead(b: Basis);
    car cdrn(b, 8);

asserted procedure basis_bget_nlead(b: Basis);
    car cdrn(b, 9);

asserted procedure basis_bset_gens(b: Basis, x);
    car cdrn(b, 1) := x;

asserted procedure basis_bset_coeffs(b: Basis, x);
    car cdrn(b, 2) := x;

asserted procedure basis_bset_size(b: Basis, x);
    car cdrn(b, 3) := x;

asserted procedure basis_bset_ndone(b: Basis, x);
    car cdrn(b, 4) := x;

asserted procedure basis_bset_ntotal(b: Basis, x);
    car cdrn(b, 5) := x;

asserted procedure basis_bset_isred(b: Basis, x);
    car cdrn(b, 6) := x;

asserted procedure basis_bset_nonred(b: Basis, x);
    car cdrn(b, 7) := x;

asserted procedure basis_bset_lead(b: Basis, x);
    car cdrn(b, 8) := x;

asserted procedure basis_bset_nlead(b: Basis, x);
    car cdrn(b, 9) := x;

% Core structure that stores basis generators and some additional info
asserted procedure basis_initializeBasis(ring: PolyRing, ngens: Integer): Basis;
    begin scalar sz;
        sz := ngens * 2;
        ndone := 0;
        ntotal := 0;
        nlead := 0;
        gens := dv_undef(sz);
        coeffs := dv_undef(sz);
        isred := dv_zeros(sz);
        nonred := dv_undef(sz);;
        lead := dv_undef(sz);
        return basis_Basis(gens, coeffs, sz, ndone, 
            ntotal, isred, nonred, lead, nlead)
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure basis_copy_basis_thorough(basis: Basis): Basis;
    begin scalar gens;
        sz := hashtable_htget_size(basis);
        gens := dv_undef(sz);
        coeffs := dv_undef(sz);
        bgens := basis_bget_gens(basis);
        bcoeffs := basis_bget_coeffs(basis);
        for i := 1:basis_bget_ntotal(basis) do <<
            putv(gens, i, dv_undef(length(getv(bcoeffs, i))));
            putv(coeffs, i, dv_undef(length(getv(bcoeffs, i))));
            gensi := getv(gens, i);
            coeffsi := getv(coeffs, i);
            for j := 1:length(getv(bcoeffs, i)) do <<
                putv(gensi, j, getv(getv(bgens, i), j));
                putv(coeffsi, j, getv(getv(bcoeffs, i), j))
            >>
        >>;
        isred := copy(basis_bget_isred(basis));
        isred := copy(basis_bget_nonred(basis));
        isred := copy(basis_bget_lead(basis));
        return basis_Basis(gens, coeffs, hashtable_htget_size(basis), basis_bget_ndone(basis),
                        basis_bget_ntotal(basis), isred, nonred, lead, basis_bget_nlead(basis))
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure basis_check_enlarge_basis(basis: Basis, added: Integer);
    begin scalar x;
        if basis_bget_ndone(basis) + added >= hashtable_htget_size(basis) then <<
            sz := basis_bget_size(basis);
            sz := max(sz*2, basis_bget_ndone(basis) + added);
            hashtable_htset_size(basis, sz);
            basis_bset_gens(basis, dv_resize(basis_bget_gens(basis), sz));
            basis_bset_coeffs(basis, dv_resize(basis_bget_coeffs(basis), sz));
            basis_bset_isred(basis, dv_resizezeros(basis_bget_isred(basis), sz));
            basis_bset_nonred(basis, dv_resize(basis_bget_nonred(basis), sz));
            basis_bset_lead(basis, dv_resize(basis_bget_lead(basis), sz));
        >>
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure basis_normalize_basis(ring: PolyRing, basis: Basis): Basis;
    begin scalar cfs;
        cfs := basis_bget_coeffs(basis);
        for i := 1:basis_bget_ntotal(basis) do <<
            if not (null getv(cfs, i)) then <<
                cfsi := getv(cfs, i);
                mul := denr(getv(cfsi, 1)) ./ numr(getv(cfsi, 1));
                for j := 2:length(cfsi) do
                    putv(cfsi, j, multsq(getv(cfsi, j), mul));
                putv(cfsi, 1) := 1 ./ 1
            >>
        >>;
        return basis
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure basis_update_pairset(
                        pairset: Pairset, basis: Basis, 
                        ht: MonomialHashtable, update_ht: MonomialHashtable, 
                        idx: Integer, plcm: Vector);
    begin scalar pl;
        pl := basis_psget_load(pairset);
        bl := idx;
        nl := pl + bl;
        ps := basis_psget_pairs(pairset);

        isred := basis_bget_isred(basis);
        new_lead := getv(getv(basis_bget_gens(basis), idx), 1);
        for i := 1:bl-1 do <<
            putv(plcm, i, hashtable_get_lcm( getv(getv(basis_bget_gens(basis), i), 1), new_lead, ht, update_ht ));
            deg := basis_spget_deg(getv(hashtable_htget_hashdata(update_ht), getv(plcm, i)));
            newidx := pl + i;
            if getv(isred, i) = 0 then
                putv(ps, newidx, basis_SPair(i, idx, getv(plcm, i), deg))
            else
                putv(ps, newidx, basis_SPair(i, idx, 0, deg))
        >>;

        % traverse existing pairs
        for i := 1:pl do <<
            psi := getv(ps, i);
            j := basis_spget_poly1(psi);
            l := basis_spget_poly2(psi);
            m := max(basis_spget_deg(getv(ps, pl + l)), basis_spget_deg(getv(ps, pl + j)));
            
            % if an existing pair is divisible by lead of new poly
            % and has a greater degree than newly generated one
            if hashtable_is_monom_divisible(basis_spget_lcmj(psi), new_lead, ht) and basis_spget_deg(psi) > m then
                putv(ps, i, basis_SPair(basis_spget_poly1(psi), basis_spget_poly2(psi), 0, basis_spget_deg(psi)))
        >>;
        
        isred := basis_bget_isred(basis);
        j := 1;
        for i := 1:bl-1 do <<
            if getv(isred, i) = 0 then <<
                putv(ps, pl + j, getv(ps, plt + i));
                j := j + 1
            >>
        >>;

        sorting_sort_pairset_by_degree(pairset, pl + 1, j - 2);

        for i := 1:j-1 do
            putv(plcm, i, basis_spget_lcmj(getv(ps, pl + i)));
        putv(plcm, j, 0);
        pc := j;
        pc := pc - 1;

        for j := 1:pc do
            if getv(plcm, j) > 0 then
                hashtable_check_monomial_division_in_update(plcm, j + 1, pc, getv(plcm, j), update_ht)
    
        j := 1;
        for i := 1:basis_psget_load(pairset) do <<
            if not (basis_spget_lcmj(getv(ps, i)) = 0) then <<
                putv(ps, j, getv(ps, i));
                j := j + 1
            >>
        >>;

        if hashtable_htget_size(ht) - hashtable_htget_load(ht) <= pc then
            hashtable_enlarge_hash_table(ht);
        
        insert_plcms_in_basis_hash_table(pairset, pl, ht, update_ht, basis, plcm, j, pc + 1);

        nonred := basis_bget_nonred(basis);
        lcm := basis_bget_nlead(basis);
        isred := basis_bget_isred(basis);
        gens := basis_bget_gens(basis);
        for i := 1:lml do
            if getv(isred, getv(nonred, i)) then
                if hashtable_is_monom_divisible(getv(getv(gens, getv(nonred, i)), 1), new_leadm ht) then
                    putv(isred, getv(nonred, i), 1)
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure basis_update_basis(basis: Basis, ht: MonomialHashtable, update_ht: MonomialHashtable);
    begin scalar k;
        k := 1;
        lead := basis_bget_lead(basis);
        nonred := basis_bget_lead(basis);
        isred := basis_bget_isred(basis);
        for i := basis_bget_nlead(basis) do
            if getv(isred, getv(nonred, i)) = 0 then <<
                putv(lead, k, getv(lead, i));
                putv(nonred, k, getv(nonred, i));
                k := k + 1
            >>;
        basis_bset_nlead(basis, k - 1);

        htdata := hashtable_htget_hashdata(ht);
        for i := basis_bget_ndone(basis)+1:basis_bget_ntotal(basis) do
            if getv(isred, i) = 0 then <<
                putv(lead, k, hashtable_hvget_divmask(getv(ht, getv(getv(gens, i), 1))));
                putv(nonred, k, i);
                k := k + 1
            >>;

        basis_bset_nlead(basis, k - 1);
        basis_bset_ndone(basis, basis_bget_ntotal(basis))
    end;

% checks if element of basis at position idx is redundant

asserted procedure basis_is_redundant(pairset: Pairset, basis: Basis, 
                        ht: MonomialHashtable, update_ht: MonomialHashtable, 
                        idx: Integer);
    begin scalar ps;
        if 2 * hashtable_htget_load(update_ht) > hashtable_htget_size(update_ht) then
            hashtable_enlarge_hash_table(update_ht);
        
        gens := basis_bget_gens(basis);
        htdata := hashtable_htget_hashdata(ht);
        isred := basis_bget_isred(basis);

        ps := basis_psget_pairs(pairset);
        % lead of new polynomial
        lead_new := getv(getv(gens), 1);
        % degree of lead
        lead_deg := hashtable_hvget_deg(getv(htdata, lead_new));

        toreturn := nil;

        for i := idx+1:basis_bget_ntotal(basis) do <<
            if not (i = (idx + 1)) and not (getv(isred, i) = 1) then <<
                lead_i := getv(getv(gens, i), 1);
                if hashtable_is_monom_divisible(lead_new, lead_i, ht) then <<
                    lcm_new := hashtable_get_lcm(lead_i, lead_new, ht, ht);
                    psidx := basis_psget_load(pairset) + 1;
                    putv(ps, psidx, basis_SPair(i, idx, lcm_new, basis_spget_deg(getv(htdata, lcm_new))));
                    putv(isred, idx, 1);
                    basis_psset_load(pairset, basis_psget_load(pairset) + 1);

                    toreturn := t
                >>
            >>
        >>;

        return toreturn
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure basis_update(pairset: Pairset, basis: Basis, 
                        ht: MonomialHashtable, update_ht: MonomialHashtable, 
                        plcm: Vector): Integer;
    begin scalar npivs;
        % number of added elements
        npivs := basis_bget_ntotal(basis);
        % number of potential critical pairs to add
        npairs := basis_bget_ndone(basis) * npivs;
        for i := 1:npivs do
            npairs := npairs + i;
        % make sure pairset and update hashtable have enough
        % space to store new pairs
        basis_check_enlarge_pairset(pairset, npairs);
        pairset_size := length(basis_psget_pairs(pairset));
        if basis_bget_ndone(basis) + 1 <= basis_bget_ntotal(basis) then <<
            % for each new element in basis
            for i := basis_bget_ndone(basis)+1:basis_bget_ntotal(basis) do <<
                if not basis_is_redundant(pairset, basis, ht, update_ht, i) then <<
                    if length(plcm) < basis_bget_ntotal(basis) + 1 then
                        plcm := dv_resize(plcm, floor(basis_bget_ntotal(basis)*1.1) + 1);
                    basis_update_pairset(pairset, basis, ht, update_ht, i, plcm)
                >>
            >>
        >>;

        basis_update_basis(basis, ht, update_ht);

        return pairset_size
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure basis_fill_data(basis: Basis, ht: MonomialHashtable, 
                                    exponents: Vector, coeffs: Vector);
    begin scalar ngens;
        ngens := length(exponents);
        htexps := hashtable_htget_exponents(ht);
        bcoeffs := basis_bget_coeffs(basis);
        bgens := basis_bget_gens(basis);
        for i := 1:ngens do <<
            while length(getv(exponents, i)) >= hashtable_htget_size(ht) - hashtable_htget_load(ht) do
                hashtable_enlarge_hash_table(ht);
            % ht.exponents one can be reallocated together with ht,
            % so we need to reset it on each iteration
            etmp := getv(htexps, 1);

            nterms := length(getv(coeffs, i));
            putv(bcoeffs, i, getv(coeffs, i));
            putv(bgens, i, dv_undef(nterms));
            poly := getv(bgens, i);
            for j := 1:nterms do
                putv(poly, j, hashtable_insert_in_hash_table(ht, getv(getv(exponents, i), j))) 
        >>;

        basis_bset_ntotal(basis, ngens)
    end;

asserted procedure basis_filter_redundant(basis: Basis);
    begin scalar j;
        j := 1;
        isred := basis_bget_isred(basis);
        nonred := basis_bget_nonred(basis);
        lead := basis_bget_lead(basis);
        for i := 1:basis_bget_nlead(basis) do <<
            if getv(isred, getv(nonred, i)) = 0 then <<
                putv(lead, j, getv(lead, i));
                putv(nonred, j, getv(nonred, i));
                j := j + 1
            >>;
        >>;
        basis_bset_nlead(j - 1);
        ASSERT(basis_bget_ndone(basis) = basis_bget_ntotal(basis));
        return basis
    end;

asserted procedure basis_standardize_basis(ring: PolyRing, basis: Basis, 
                                            ht: MonomialHashtable, ord);
    begin scalar sz;
        nonred := basis_bget_nonred(basis);
        isred := basis_bget_isred(basis);
        lead := basis_bget_lead(basis);
        bcoeffs := basis_bget_coeffs(basis);
        bgens := basis_bget_gens(basis);
        for i := 1:basis_bget_nlead(basis) do <<
            idx := getv(nonred, i);
            putv(nonred, i, i);
            putv(isred, i, 0);
            putv(bcoeffs, getv(coeffs, idx));
            putv(bgens, i, getv(bgens, idx))
        >>;
        nlead := basis_bget_nlead(basis);
        basis_bset_ntotal(basis, nlead);
        basis_bset_ndone(basis, nlead);
        basis_bset_size(basis, nlead);
        
        basis_bset_coeffs(basis, dv_resize(bcoeffs, nlead));
        basis_bset_gens(basis, dv_resize(bgens, nlead));
        basis_bset_lead(basis, dv_resize(lead, nlead));
        basis_bset_nonred(basis, dv_resize(nonred, nlead));
        basis_bset_isred(basis, dv_resize(isred, nlead));

        sorting_sort_gens_by_lead_increasing_in_standardize(basis, ht, ord);
        return basis_normalize_basis(ring, basis)
    end;

asserted procedure basis_export_basis_data(basis: Basis, ht: MonomialHashtable);
    begin scalar exps;
        nlead := basis_bget_nlead(basis);
        exps := dv_undef(nlead);
        coeffs := dv_undef(nlead);

        htexps := hashtable_htget_exponents(ht);
        nonred := basis_bget_nonred(basis);
        gens := basis_bget_gens(basis);
        coeffs := basis_bget_coeffs(basis);
        for i := 1:nlead do <<
            idx := getv(nonred, i);
            poly := getv(gens, idx);
            putv(exps, i, dv_undef(length(poly)));
            expsi := getv(exps, i);
            for j := 1:length(poly) do
                putv(expsi, j, getv(htexps, getv(poly, j)));
            putv(coeffs, i, getv(coeffs, idx))
        >>;

        return exps . coeffs
    end;

asserted procedure basis_hash_to_exponents(basis: Basis, ht: MonomialHashtable);
    begin scalar exps;
        nlead := basis_bget_nlead(basis);
        exps := dv_undef(nlead);

        htexps := hashtable_htget_exponents(ht);
        nonred := basis_bget_nonred(basis);
        gens := basis_bget_gens(basis);

        for i := 1:nlead do <<
            idx := getv(nonred, i);
            poly := getv(gens, idx);
            putv(exps, i, dv_undef(length(poly)));
            expsi := getv(exps, i);
            for j := 1:length(poly) do
                putv(expsi, j, getv(htexps, getv(poly, j)));
        >>;

        return exps
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure basis_export_basis_data(pairset: Pairset, off: Integer, 
                    ht: MonomialHashtable, update_ht: MonomialHashtable, 
                    basis: Basis, plcm: Vector, ifirst: Integer, ilast: Integer);
    begin scalar x;
        % including ifirst and not including ilast

        explen := hashtable_htget_explen(ht);

        gens := basis_bget_gens(basis);
        mod := hashtable_htget_size(ht) - 1;
        ps := basis_psget_pairs(pairset);

        upddata := hashtable_htget_hashdata(update_ht);
        updexps := hashtable_htget_exponents(update_ht);
        htexps := hashtable_htget_exponents(ht);
        httable := hashtable_htget_hashtable(ht);
        htdata := hashtable_htget_hashdata(ht);

        m := ifirst;
        l := 1;
    Letsgo:
        while l < ilast do <<
            genspoly1 := getv(getv(gens, basis_spget_poly1(getv(ps, off+l))), 1);
            genspoly2 := getv(getv(gens, basis_spget_poly2(getv(ps, off+1))), 1);
            if getv(plcm, l) = 0 or hashtable_is_gcd_const(genspoly1, genspoly2, ht) then
                l := l + 1
            else <<
                putv(ps, m, getv(ps, off + l));
                h := hashtable_hvget_hash(getv(upddata, getv(plcm, l)));
                putv(htexps, hashtable_htget_load(ht) + 1, copy(getv(updexps, getv(plcm, l))));
                n := getv(htexps, hashtable_htget_load(ht) + 1);

                k := h;
                i := 1;
            Restart:
                while i <= hashtable_htget_size(ht) <<
                    k := hashtable_hashnextindex(h, i, mod);
                    hm := getv(httable, k);
                    if hm = 0 then
                        i := hashtable_htget_size(ht) + 1
                    else <<
                        if hashtable_hvget_hash(getv(htdata, hm)) = h then
                            i := i + 1
                        else <<
                            ehm := getv(htexps, hm);
                            
                            for j := 1:explen do
                                if not (getv(ehm, j) = getv(n, j)) then <<
                                    i := i + 1;
                                    go to Restart
                                >>;
                            
                            putv(ps, m, basis_SPair, basis_spget_poly1(getv(ps, m)), basis_spget_poly2(getv(ps, m)), hm, hashtable_hvget_deg(getv(ps, m)));
                            m := m + 1;
                            l := l + 1;
                            go to Letsgo
                        >>
                    >>
                >>
            >>;

            pos := hashtable_htget_load(ht) + 1;
            putv(httable, k, pos);
            ll := getv(plcm, l);
            putv(htdata, pos, hashtable_Hashvalue(h, hashtable_hvget_divmask(getv(upddata, ll), 0, hashtable_hvget_deg(getv(upddata, ll)))));

            hashtable_htset_load(hashtable_htget_load(ht) + 1);
            putv(ps, m, basis_SPair(basis_spget_poly1(getv(ps, m)), basis_spget_poly2(getv(ps, m)), pos, basis_spget_deg(getv(ps, m))));
            m := m + 1;
            l := l + 1
        >>;
        
        hashtable_htset_load(m - 1)
    end;

%--------------------------------------------------------------------------------------------------

end;
