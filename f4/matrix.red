
%--------------------------------------------------------------------------------------------------

asserted procedure matrix_MacaulayMatrix(uprows, lowrows, col2hash, coeffs, size, npivots, 
                            nrows, ncols, nup, nlow, nleft, nright, up2coef, low2coef);
    {'matrix, uprows, lowrows, col2hash, coeffs, size, npivots, 
        nrows, ncols, nup, nlow, nleft, nright, up2coef, low2coef};

asserted procedure matrix_mget_uprows(m: MacaulayMatrix);
    car cdrn(m, 1);

asserted procedure matrix_mget_lowrows(m: MacaulayMatrix);
    car cdrn(m, 2);

asserted procedure matrix_mget_col2hash(m: MacaulayMatrix);
    car cdrn(m, 3);

asserted procedure matrix_mget_coeffs(m: MacaulayMatrix);
    car cdrn(m, 4);

asserted procedure matrix_mget_size(m: MacaulayMatrix);
    car cdrn(m, 5);

asserted procedure matrix_mget_npivots(m: MacaulayMatrix);
    car cdrn(m, 6);

asserted procedure matrix_mget_nrows(m: MacaulayMatrix);
    car cdrn(m, 7);

asserted procedure matrix_mget_ncols(m: MacaulayMatrix);
    car cdrn(m, 8);

asserted procedure matrix_mget_nup(m: MacaulayMatrix);
    car cdrn(m, 9);

asserted procedure matrix_mget_nlow(m: MacaulayMatrix);
    car cdrn(m, 10);

asserted procedure matrix_mget_nleft(m: MacaulayMatrix);
    car cdrn(m, 11);

asserted procedure matrix_mget_nright(m: MacaulayMatrix);
    car cdrn(m, 11);

asserted procedure matrix_mget_up2coef(m: MacaulayMatrix);
    car cdrn(m, 12);
    
asserted procedure matrix_mget_low2coef(m: MacaulayMatrix);
    car cdrn(m, 13);

asserted procedure matrix_mset_uprows(m: MacaulayMatrix, x);
    car cdrn(m, 1) := x;

asserted procedure matrix_mset_lowrows(m: MacaulayMatrix, x);
    car cdrn(m, 2) := x;

asserted procedure matrix_mset_col2hash(m: MacaulayMatrix, x);
    car cdrn(m, 3) := x;

asserted procedure matrix_mset_coeffs(m: MacaulayMatrix, x);
    car cdrn(m, 4) := x;

asserted procedure matrix_mset_size(m: MacaulayMatrix, x);
    car cdrn(m, 5) := x;

asserted procedure matrix_mset_npivots(m: MacaulayMatrix, x);
    car cdrn(m, 6) := x;

asserted procedure matrix_mset_nrows(m: MacaulayMatrix, x);
    car cdrn(m, 7) := x;

asserted procedure matrix_mset_ncols(m: MacaulayMatrix, x);
    car cdrn(m, 8) := x;

asserted procedure matrix_mset_nup(m: MacaulayMatrix, x);
    car cdrn(m, 9) := x;

asserted procedure matrix_mset_nlow(m: MacaulayMatrix, x);
    car cdrn(m, 10) := x;

asserted procedure matrix_mset_nleft(m: MacaulayMatrix, x);
    car cdrn(m, 11) := x;

asserted procedure matrix_mset_nright(m: MacaulayMatrix, x);
    car cdrn(m, 12) := x;

asserted procedure matrix_mset_up2coef(m: MacaulayMatrix, x);
    car cdrn(m, 13) := x;
    
asserted procedure matrix_mset_low2coef(m: MacaulayMatrix, x);
    car cdrn(m, 14) := x;

asserted procedure matrix_initialize_matrix(ring: PolyRing): MacaulayMatrix;
    being scalar uprows;
        uprows := dv_undef(0);
        lowrows := dv_undef(0);
        col2hash := dv_undef(0);
        coeffs := dv_undef(0);

        size := 0;
        npivots := 0;
        nrows := 0;
        ncols := 0;

        nup := 0;
        nlow := 0;
        nleft := 0;
        nright := 0;

        up2coef := dv_undef(0);
        low2coef := dv_undef(0);

        return matrix_MacaulayMatrix(uprows, lowrows, col2hash, coeffs, size, npivots, 
        nrows, ncols, nup, nlow, nleft, nright, up2coef, low2coef)
    end;

asserted procedure matrix_reinitialize_matrix(matrix: MacaulayMatrix, npairs: Integer): MacaulayMatrix;
    begin scalar x;
        matrix_mset_uprows(matrix, dv_resize(matrix_mget_uprows(matrix), npairs*2));
        matrix_mset_lowrows(matrix, dv_resize(matrix_mget_lowrows(matrix), npairs*2));
        matrix_mset_up2coef(matrix, dv_resize(matrix_mget_up2coef(matrix), npairs*2));
        matrix_mset_low2coef(matrix, dv_resize(matrix_mget_low2coef(matrix), npairs*2));

        matrix_mset_size(matrix, 2*npairs);
        matrix_mset_ncols(matrix, 0);
        matrix_mset_nleft(matrix, 0);
        matrix_mset_nright(matrix, 0);
        matrix_mset_nup(matrix, 0);
        matrix_mset_nlow(matrix, 0);
        return matrix
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_normalize_sparse_row(row: Vector, ch): Vector;
    begin scalar pinv;
        pinv := denr(getv(row, 1)) ./ numr(getv(row, 1));
        for i := 2:length(row) do
            putv(row, i, multsq(getv(row, i), pinv));
        putv(row, 1, 1 ./ 1);
        return row
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_reduce_by_pivot(row: Vector, indices: Vector, cfs: Vector, magic);
    begin scalar mul;
        mul := negsq(getv(row, getv(indices, 1)));
        
        for j := 1:length(indices) do <<
            idx := getv(indices, j);
            putv(row, idx, addsq(getv(row, idx), multsq(mul, getv(cfs, j))))
        >>;
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_load_indexed_coefficients(densecoeffs: Vector, rowexps: Vector, cfsref: Vector);
    begin;
        for i := 1:length(densecoeffs) do
            putv(densecoeffs, nil . 1);
        for j := 1:length(rowexps) do
            putv(densecoeffs, getv(rowexps, j), getv(cfsref, j))
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure reduce_dense_row_by_known_pivots_sparse(densecoeffs: Vector, 
                        matrix: Matrix, basis: Basis, pivs: Vector, startcol: ColumnIdx,
                        tmp_pos: ColumnIdx, magic, exact_colmap);
    begin scalar ncols;
        
        ncols := matrix_mget_ncols(matrix);
        nleft := matrix_mget_nleft(matrix);

        mcoeffs := matrix_mget_coeffs(matrix);
        up2coef := matrix_mget_up2coef(matrix);
        low2coef := matrix_mset_low2coef(matrix);

        bcoeffs := basis_bget_coeffs(basis);

        % new row nonzero elements count
        k := 0;
        uzero := nil . 1;

        % new pivot index
        np := -1;

        for i := startcol:ncols do <<

            % if row element zero - no reduction
            if not (getv(densecoeffs, i) = uzero) then <<
                if null getv(pivs, i) or (not (tmp_pos = -1) and tmp_pos = i) then <<
                    if np = 1 then
                        np := i;
                    k := k + 1
                >> else <<
                    reducerexps := getv(pivs, i);
                    if exact_colmap
                        cfs := getv(mcoeffs, tmp_pos)
                    else if i <= nleft
                        cfs := getv(bcoeffs, getv(up2coef, i))
                    else
                        cfs := getv(mcoeffs, getv(low2coef, i))
                    >>
                    matrix_reduce_by_pivot(densecoeffs, reducerexps, cfs, magic)
                >>
            >>
        >>;

        newrow := dv_undef(k);
        newcfs := dv_undef(k);

        % all reduced!
        if k = 0 then
            return {t, newrow, newcfs};
        
        j := 1
        for i := np:ncols do <<
            if getv(densecoeffs, i) = uzero then <<
                putv(newrow, j, i);
                putv(newcfs, j, getv(densecoeffs, i));
                j := j + 1
            >>
        >>

        return {nil, newrow, newcfs}
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_exact_sparse_rref(ring: PolyRing, matrix: MacaulayMatrix, basis: Basis);
    begin scalar ncols;
        ncols := matrix_mget_ncols(matrix);
        nlow := matrix_mget_nlow(matrix);
        nright := matrix_mget_nright(matrix);
        nleft := matrix_mget_nleft(matrix);

        uprows := matrix_mget_uprows(matrix);
        lowrows := matrix_mget_low(matrix);
        low2coef := matrix_mget_low2coef(matrix);
        up2coef := matrix_mget_up2coef(matrix);

        % Julia
        % magic :=

        pivs := dv_undef(ncols);
        for i := 1:matrix_mget_nup(matrix) do
            putv(pivs, i, getv(uprows, i));

        l2c_tmp := dv_undef(max(ncols, nlow));
        for i := 1:nlow do
            putv(l2c_tmp, getv(getv(lowrows, i), 1), getv(low2coef, i));

        rowidx2coef := low2coef;
        matrix_mset_low2coef(l2c_tmp);

        upivs := matrix_mget_lowrows(matrix);

        densecoeffs := dv_zeros(ncols);

        bcoeffs := basis_bget_coeffs(basis);

        for i := 1:nlow do <<
            % select next row to be reduced
            rowexps := getv(upivs, i);
            
            % corresponding coefficients from basis
            % (no need to copy here)
            cfsref := getv(bcoeffs, getv(rowidx2coef, i));
            
            % we load coefficients into dense array
            % into rowexps indices
            matrix_load_indexed_coefficients(densecoeffs, rowexps, cfsref);

            % reduce it with known pivots from matrix.uprows
            % first nonzero in densecoeffs is at startcol position
            startcol = getv(rowexps, 1);

            result = matrix_reduce_dense_row_by_known_pivots_sparse(densecoeffs, matrix, basis, pivs, startcol, -1, magic, nil);
            zeroed := pop(result); 
            newrow := pop(result); 
            newcfs := pop(result);

            % if not fully reduced
            if not zeroed then <<
                mcoeffs := matrix_mget_coeffs(matrix);
                low2coef := matrix_mget_low2coef(matrix);

                putv(mcoeffs, i, newcfs);
                
                putv(pivs, getv(newrow, 1)) := newrow;

                putv(low2coef, getv(newrow, 1)) := i;

                if not (getv(getv(mcoeffs, i), 1) = 1 ./ 1) then
                    matrix_normalize_sparse_row(getv(mcoeffs, i), magic)
            >>
        >>;

        newpivs := 0;



    end;

%--------------------------------------------------------------------------------------------------

end; % of file