module f4matrix;
% The matrixj used in f4-style reduction.

%--------------------------------------------------------------------------------------------------

% struct MacaulayMatrix
%
% Provides the `MacaulayMatrix` interface, which implements Macaulay matrixj
% of the folliwing structure, where A,B,C,D are sparse matrices itself.
%        | A  B |
%        | C  D |
% Part A contains known pivots of reducing rows,
%        and part CD are rows to be reduced by AB.
%
%
% The matrixj stores the following elements:
%   . uprows - rows from upper, AB part of the matrixj
%   . lowrows - rows from lower, CD part of the matrixj,
%   . col2hash - maps column idx {1 ... ncols} to monomial hash idx {2 ... load}
%   . coeffs - rows coefficients
%   . size - total number of allocated tows
%   . npivots - the number of new basis elements discovered during matrixj reduction
%   . nrows - number of filled rows
%   . ncols - number of filled columns
%   . nup - number of upper rows (rows in AB section)
%   . nlow - number of lower rows (rows in CD section)
%   . nleft - number of right rows (rows in AC section)
%   . nright - number of right rows (rows in BD section)
asserted procedure matrix_MacaulayMatrix(uprows: Vector, lowrows: Vector, col2hash: Vector, coeffs: Vector, size: Integer, npivots: Integer, 
                            nrows: Integer, ncols: Integer, nup: Integer, nlow: Integer, nleft: Integer, nright: Integer, up2coef: Vector, low2coef: Vector): MacaulayMatrix;
    begin scalar v;
        v := dv_undef(14);
        putv(v, 1, uprows);
        putv(v, 2, lowrows);
        putv(v, 3, col2hash);
        putv(v, 4, coeffs);
        putv(v, 5, size);
        putv(v, 6, npivots);
        putv(v, 7, nrows);
        putv(v, 8, ncols);
        putv(v, 9, nup);
        putv(v, 10, nlow);
        putv(v, 11, nleft);
        putv(v, 12, nright);
        putv(v, 13, up2coef);
        putv(v, 14, low2coef);
        return v
    end;

asserted procedure matrix_mget_uprows(m: MacaulayMatrix): Vector;
    getv(m, 1);

asserted procedure matrix_mget_lowrows(m: MacaulayMatrix): Vector;
    getv(m, 2);

asserted procedure matrix_mget_col2hash(m: MacaulayMatrix): Vector;
    getv(m, 3);

asserted procedure matrix_mget_coeffs(m: MacaulayMatrix): Vector;
    getv(m, 4);

asserted procedure matrix_mget_size(m: MacaulayMatrix): Integer;
    getv(m, 5);

asserted procedure matrix_mget_npivots(m: MacaulayMatrix): Integer;
    getv(m, 6);

asserted procedure matrix_mget_nrows(m: MacaulayMatrix): Integer;
    getv(m, 7);

asserted procedure matrix_mget_ncols(m: MacaulayMatrix): Integer;
    getv(m, 8);

asserted procedure matrix_mget_nup(m: MacaulayMatrix): Integer;
    getv(m, 9);

asserted procedure matrix_mget_nlow(m: MacaulayMatrix): Integer;
    getv(m, 10);

asserted procedure matrix_mget_nleft(m: MacaulayMatrix): Integer;
    getv(m, 11);

asserted procedure matrix_mget_nright(m: MacaulayMatrix): Integer;
    getv(m, 12);

asserted procedure matrix_mget_up2coef(m: MacaulayMatrix): Vector;
    getv(m, 13);
    
asserted procedure matrix_mget_low2coef(m: MacaulayMatrix): Vector;
    getv(m, 14);

asserted procedure matrix_mset_uprows(m: MacaulayMatrix, x);
    putv(m, 1, x);

asserted procedure matrix_mset_lowrows(m: MacaulayMatrix, x);
    putv(m, 2, x);

asserted procedure matrix_mset_col2hash(m: MacaulayMatrix, x);
    putv(m, 3, x);

asserted procedure matrix_mset_coeffs(m: MacaulayMatrix, x);
    putv(m, 4, x);

asserted procedure matrix_mset_size(m: MacaulayMatrix, x);
    putv(m, 5, x);

asserted procedure matrix_mset_npivots(m: MacaulayMatrix, x);
    putv(m, 6, x);

asserted procedure matrix_mset_nrows(m: MacaulayMatrix, x);
    putv(m, 7, x);

asserted procedure matrix_mset_ncols(m: MacaulayMatrix, x);
    putv(m, 8, x);

asserted procedure matrix_mset_nup(m: MacaulayMatrix, x);
    putv(m, 9, x);

asserted procedure matrix_mset_nlow(m: MacaulayMatrix, x);
    putv(m, 10, x);

asserted procedure matrix_mset_nleft(m: MacaulayMatrix, x);
    putv(m, 11, x);

asserted procedure matrix_mset_nright(m: MacaulayMatrix, x);
    putv(m, 12, x);

asserted procedure matrix_mset_up2coef(m: MacaulayMatrix, x);
    putv(m, 13, x);
    
asserted procedure matrix_mset_low2coef(m: MacaulayMatrix, x);
    putv(m, 14, x);

asserted procedure matrix_initialize_matrix(ring: PolyRing): MacaulayMatrix;
    begin scalar uprows, lowrows, col2hash, coeffs, size, npivots,
                    nrows, ncols, nup, nlow, nleft, nright, up2coef, low2coef;
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

asserted procedure matrix_reinitialize_matrix(matrixj: MacaulayMatrix, npairs: Integer): MacaulayMatrix;
    begin
        matrix_mset_uprows(matrixj, dv_resize(matrix_mget_uprows(matrixj), npairs*2));
        matrix_mset_lowrows(matrixj, dv_resize(matrix_mget_lowrows(matrixj), npairs*2));
        matrix_mset_up2coef(matrixj, dv_resize(matrix_mget_up2coef(matrixj), npairs*2));
        matrix_mset_low2coef(matrixj, dv_resize(matrix_mget_low2coef(matrixj), npairs*2));

        matrix_mset_size(matrixj, 2*npairs);
        matrix_mset_ncols(matrixj, 0);
        matrix_mset_nleft(matrixj, 0);
        matrix_mset_nright(matrixj, 0);
        matrix_mset_nup(matrixj, 0);
        matrix_mset_nlow(matrixj, 0);
        return matrixj
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_zero_coeff_vector(ring, ncols);
    begin scalar v, ch;
        v := dv_undef(ncols);
        ch := io_prget_ch(ring);
        for i := 1:dv_length(v) do
            putv(v, i, if ch = 0 then (nil ./ 1) else 0);
        return v
    end;

%--------------------------------------------------------------------------------------------------

% Normalize `row` by first coefficient
asserted procedure matrix_normalize_sparse_row(ring: PolyRing, row: Vector): Vector;
    if io_prget_ch(ring) = 0 then
        matrix_normalize_sparse_row_qq(ring, row)
    else
        matrix_normalize_sparse_row_ff(ring, row);

% Normalize `row` by first coefficient
asserted procedure matrix_normalize_sparse_row_qq(ring: PolyRing, row: Vector): Vector;
    begin scalar pinv;
        pinv := denr(getv(row, 1)) ./ numr(getv(row, 1));
        for i := 2:dv_length(row) do
            putv(row, i, multsq(getv(row, i), pinv));
        putv(row, 1, 1 ./ 1);
        return row
    end;

% Normalize `row` by first coefficient
asserted procedure matrix_normalize_sparse_row_ff(ring: PolyRing, row: Vector): Vector;
    begin scalar ch, pinv;
        ch := io_prget_ch(ring);
        pinv := remainder(modular_invmod(getv(row, 1), ch), ch);
        ASSERT(pinv > 0);

        for i := 2:dv_length(row) do
            putv(row, i, remainder(getv(row, i) #* pinv, ch));
        putv(row, 1, 1);
        return row
    end;

%--------------------------------------------------------------------------------------------------

% reduces row by mul*cfs at indices positions
asserted procedure matrix_reduce_by_pivot(ring: PolyRing, row: Vector, indices: Vector, cfs: Vector);
    if io_prget_ch(ring) = 0 then
        matrix_reduce_by_pivot_qq(ring, row, indices, cfs)
    else
        matrix_reduce_by_pivot_ff(ring, row, indices, cfs);

% reduces row by mul*cfs at indices positions
asserted procedure matrix_reduce_by_pivot_qq(ring: PolyRing, row: Vector, indices: Vector, cfs: Vector);
    begin scalar mul, idx;
        mul := negsq(getv(row, getv(indices, 1)));
        
        for j := 1:dv_length(indices) do <<
            idx := getv(indices, j);
            putv(row, idx, addsq(getv(row, idx), multsq(mul, getv(cfs, j))))
        >>;
    end;

% reduces row by mul*cfs at indices positions
asserted procedure matrix_reduce_by_pivot_ff(ring: PolyRing, row: Vector, indices: Vector, cfs: Vector);
    begin scalar mul, idx, ch;
        ch := io_prget_ch(ring);
        mul := ch #- getv(row, getv(indices, 1));
        ASSERT(mul > 0);

        for j := 1:dv_length(indices) do <<
            idx := getv(indices, j);
            putv(row, idx, remainder(getv(row, idx) #+ mul #* getv(cfs, j), ch))
        >>;
    end;

%--------------------------------------------------------------------------------------------------

% zero entries of densecoeffs and load coefficients cfsref to indices rowexps
asserted procedure matrix_load_indexed_coefficients(ring: PolyRing, densecoeffs: Vector, rowexps: Vector, cfsref: Vector);
    if io_prget_ch(ring) = 0 then
        matrix_load_indexed_coefficients_qq(ring, densecoeffs, rowexps, cfsref)
    else
        matrix_load_indexed_coefficients_ff(ring, densecoeffs, rowexps, cfsref);

% zero entries of densecoeffs and load coefficients cfsref to indices rowexps
asserted procedure matrix_load_indexed_coefficients_qq(ring: PolyRing, densecoeffs: Vector, rowexps: Vector, cfsref: Vector);
    <<
        for i := 1:dv_length(densecoeffs) do
            putv(densecoeffs, i, nil . 1);
        for j := 1:dv_length(rowexps) do
            putv(densecoeffs, getv(rowexps, j), getv(cfsref, j))
    >>;

% zero entries of densecoeffs and load coefficients cfsref to indices rowexps
asserted procedure matrix_load_indexed_coefficients_ff(ring: PolyRing, densecoeffs: Vector, rowexps: Vector, cfsref: Vector);
    <<
        for i := 1:dv_length(densecoeffs) do
            putv(densecoeffs, i, 0);
        for j := 1:dv_length(rowexps) do
            putv(densecoeffs, getv(rowexps, j), getv(cfsref, j))
    >>;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_reduce_dense_row_by_known_pivots_sparse(ring: PolyRing, densecoeffs: Vector, 
                        matrixj: MacaulayMatrix, basis: Basis, pivs: Vector, startcol: ColumnIdx,
                        tmp_pos: ColumnIdx, exact_colmap: Boolean);
    begin scalar ncols, nleft, mcoeffs, up2coef, low2coef, bcoeffs, k, uzero,
                    np, cfs, j, np, newrow, newcfs, reducerexps;
        
        ncols := matrix_mget_ncols(matrixj);
        nleft := matrix_mget_nleft(matrixj);

        mcoeffs := matrix_mget_coeffs(matrixj);
        up2coef := matrix_mget_up2coef(matrixj);
        low2coef := matrix_mget_low2coef(matrixj);

        bcoeffs := basis_bget_coeffs(basis);

        % new row nonzero elements count
        k := 0;
        uzero := if io_prget_ch(ring) = 0 then (nil ./ 1) else 0;

        % new pivot index
        np := -1;

        for i := startcol:ncols do <<

            % if row element zero - no reduction
            if getv(densecoeffs, i) neq uzero then <<
                if null getv(pivs, i) or ((tmp_pos neq -1) and (tmp_pos = i)) then <<
                    if np = -1 then
                        np := i;
                    k := k + 1
                >> else <<

                    % exponents of reducer row at column i
                    reducerexps := getv(pivs, i);
                    
                    if exact_colmap then     % if exact mapping in matrixj.low2coef
                        cfs := getv(mcoeffs, tmp_pos)
                    else if i <= nleft then   % if reducer is from upper
                        cfs := getv(bcoeffs, getv(up2coef, i))
                    else                 % -//- of lower part of matrixj
                        cfs := getv(mcoeffs, getv(low2coef, i));

                    matrix_reduce_by_pivot(ring, densecoeffs, reducerexps, cfs)
                >>
            >>
        >>;

        newrow := dv_undef(k);
        newcfs := dv_undef(k);

        % all reduced!
        if k = 0 then
            return {t, newrow, newcfs};
        
        % store new row in sparse format
        % where k - number of structural nonzeros in new reduced row, k > 0
        j := 1;
        for i := np:ncols do <<  % from new pivot
            if getv(densecoeffs, i) neq uzero then <<
                putv(newrow, j, i);
                putv(newcfs, j, getv(densecoeffs, i));
                j := j + 1
            >>
        >>;

        return {nil, newrow, newcfs}
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_exact_sparse_rref(ring: PolyRing, matrixj: MacaulayMatrix, basis: Basis);
    begin scalar ncols, nlow, nright, nleft, uprows, lowrows, low2coef, up2coef, pivs, l2c_tmp,
                    rowidx2coef, upivs, densecoeffs, bcoeffs, rowexps, cfsref, startcol, 
                    result, zeroed, newrow, newcfs, mcoeffs, newpivs, k;
        ncols := matrix_mget_ncols(matrixj);
        nlow := matrix_mget_nlow(matrixj);
        nright := matrix_mget_nright(matrixj);
        nleft := matrix_mget_nleft(matrixj);

        uprows := matrix_mget_uprows(matrixj);
        lowrows := matrix_mget_lowrows(matrixj);
        low2coef := matrix_mget_low2coef(matrixj);
        up2coef := matrix_mget_up2coef(matrixj);

        % Julia
        % magic :=

        % known pivots
        % no_copy
        pivs := dv_undef(ncols);
        for i := 1:matrix_mget_nup(matrixj) do
            putv(pivs, i, getv(uprows, i));

        l2c_tmp := dv_undef(max(ncols, nlow));
        for i := 1:nlow do
            putv(l2c_tmp, getv(getv(lowrows, i), 1), getv(low2coef, i));

        % no_copy
        rowidx2coef := low2coef;
        matrix_mset_low2coef(matrixj, l2c_tmp);

        % unknown pivots
        % (not discovered yet)
        % we will modify them inplace when reducing by pivs
        upivs := matrix_mget_lowrows(matrixj);

        densecoeffs := matrix_zero_coeff_vector(ring, ncols);

        bcoeffs := basis_bget_coeffs(basis);

        for i := 1:nlow do <<
            % select next row to be reduced
            rowexps := getv(upivs, i);
            
            % corresponding coefficients from basis
            % (no need to copy here)
            cfsref := getv(bcoeffs, getv(rowidx2coef, i));
            
            % we load coefficients into dense array
            % into rowexps indices
            matrix_load_indexed_coefficients(ring, densecoeffs, rowexps, cfsref);

            % reduce it with known pivots from matrixj.uprows
            % first nonzero in densecoeffs is at startcol position
            startcol := getv(rowexps, 1);

            {zeroed, newrow, newcfs} := matrix_reduce_dense_row_by_known_pivots_sparse(ring, densecoeffs, matrixj, basis, pivs, startcol, -1, nil);

            % if not fully reduced
            if not zeroed then <<
                mcoeffs := matrix_mget_coeffs(matrixj);
                low2coef := matrix_mget_low2coef(matrixj);

                % matrixj coeffs sparsely stores coefficients of new row
                putv(mcoeffs, i, newcfs);
                
                % add new pivot at column index newrow[1]
                % (which is the first nnz column of newrow)
                putv(pivs, getv(newrow, 1), newrow);

                % set ref to coefficient to matrixj
                % guaranteed to be from lower part
                putv(low2coef, getv(newrow, 1), i);

                matrix_normalize_sparse_row(ring, getv(mcoeffs, i))
            >>
        >>;
        
        if f4_debug() then <<
            prin2t {"exact_sparse_rref: lower part reduced:", matrixj};
            prin2t {"exact_sparse_rref: pivs:", pivs}
        >>;

        % number of new pivots
        newpivs := 0;

        % a row to be reduced for each column
        lowrows := matrix_mget_lowrows(matrixj);
        lowrows := dv_resize(lowrows, matrix_mget_nright(matrixj));
        matrix_mset_lowrows(matrixj, lowrows);

        bcoeffs := basis_bget_coeffs(basis);
        mcoeffs := matrix_mget_coeffs(matrixj);

        % interreduce new pivots..
        % for each right (non-pivotal) column
        for i := 1:nright do <<
            k := ncols - i + 1;
            if not null getv(pivs, k) then <<

                if k <= nleft then
                    cfsref := getv(bcoeffs, getv(up2coef, k))
                else   % of lower part of matrixj
                    cfsref := getv(mcoeffs, getv(low2coef, k));
                
                startcol := getv(getv(pivs, k), 1);

                matrix_load_indexed_coefficients(ring, densecoeffs, getv(pivs, k), cfsref);
                
                newpivs := newpivs + 1;

                {zeroed, newrow, newcfs} := matrix_reduce_dense_row_by_known_pivots_sparse(ring, densecoeffs, matrixj, basis, pivs, startcol, startcol, nil);

                % update row and coeffs
                putv(lowrows, newpivs, newrow);
                putv(mcoeffs, getv(low2coef, k), newcfs);
                putv(low2coef, k, getv(low2coef, k));
                putv(pivs, k, getv(lowrows, newpivs))
            >>
        >>;

        % shrink matrixj
        matrix_mset_size(matrixj, newpivs);
        matrix_mset_nrows(matrixj, newpivs);
        matrix_mset_npivots(matrixj, newpivs);
        matrix_mset_lowrows(matrixj, dv_resize(matrix_mget_lowrows(matrixj), newpivs))
    end;

asserted procedure matrix_linear_algebra(ring: PolyRing, matrixj: MacaulayMatrix, basis: Basis);
    begin scalar coeffs;
        coeffs := matrix_mget_coeffs(matrixj);
        matrix_mset_coeffs(matrixj, dv_resize(coeffs, matrix_mget_nlow(matrixj)));
        matrix_exact_sparse_rref(ring, matrixj, basis)
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_interreduce_matrix_rows(ring: PolyRing, matrixj: MacaulayMatrix, basis: Basis);
    begin scalar pivs, densecfs, uprows, lowrows, low2coef, up2coef, coeffs, k,
                    l, cfs, reducexps, startcol, zeroes, newrow, newcfs;
        matrix_mset_lowrows(matrixj, dv_resize(matrix_mget_lowrows(matrixj), matrix_mget_ncols(matrixj)));
        matrix_mset_up2coef(matrixj, dv_resize(matrix_mget_up2coef(matrixj), matrix_mget_ncols(matrixj)));
        matrix_mset_low2coef(matrixj, dv_resize(matrix_mget_low2coef(matrixj), matrix_mget_ncols(matrixj)));
        matrix_mset_coeffs(matrixj, dv_resize(matrix_mget_coeffs(matrixj), matrix_mget_ncols(matrixj)));

        % Julia..
        % magic = 

        uprows := matrix_mget_uprows(matrixj);
        lowrows := matrix_mget_lowrows(matrixj);
        low2coef := matrix_mget_low2coef(matrixj);
        up2coef := matrix_mget_up2coef(matrixj);
        coeffs := matrix_mget_coeffs(matrixj);

        % prin2t {"In interreduce"};
        % prin2t {"Matrix", matrixj};
        % prin2t {"Basis", basis};

        % same pivs as for rref
        % pivs: column idx --> vector of present columns
        pivs := dv_undef(matrix_mget_ncols(matrixj));
        for i := 1:matrix_mget_nrows(matrixj) do <<
            putv(pivs, getv(getv(uprows, i), 1), getv(uprows, i));
            putv(low2coef, getv(getv(uprows, i), 1), i);
            putv(coeffs, i, dv_copy(getv(basis_bget_coeffs(basis), getv(up2coef, i))))
        >>;

        densecfs := matrix_zero_coeff_vector(ring, matrix_mget_ncols(matrixj));

        k := matrix_mget_nrows(matrixj);

        % prin2t "uwu";

        for i := 1:matrix_mget_ncols(matrixj) do <<
            l := matrix_mget_ncols(matrixj) - i + 1;
            if not (null getv(pivs, l)) then <<
                cfs := getv(coeffs, getv(low2coef, l));
                reducexps := getv(pivs, l);

                startcol := getv(reducexps, 1);

                for j := 1:dv_length(reducexps) do
                    putv(densecfs, getv(reducexps, j), getv(cfs, j));

                % prin2t {"Reducting, i =", i};
                % prin2t {"densecfs = ", densecfs};
                % prin2t {"pivs = ", pivs};
                % prin2t {"l = ", l};
                {zeroes, newrow, newcfs} := matrix_reduce_dense_row_by_known_pivots_sparse(ring, densecfs, matrixj, basis, pivs, l, l, nil);

                putv(lowrows, k, newrow);
                putv(coeffs, getv(low2coef, l), newcfs);
                putv(pivs, l, getv(lowrows, k));
                k := k - 1;

                for idx := 1:dv_length(newrow) do
                    putv(densecfs, getv(newrow, idx), if io_prget_ch(ring) = 0 then nil ./ 1 else 0)
            >>
        >>;

        matrix_mset_npivots(matrixj, matrix_mget_nrows(matrixj))
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_convert_hashes_to_columns(matrixj: MacaulayMatrix, 
                            symbol_ht: MonomialHashtable);
    begin scalar hdata, loadj, col2hash, j, k, nterms, uprows, row,
                    lowrows;

        hdata := hashtable_htget_hashdata(symbol_ht);
        loadj := hashtable_htget_load(symbol_ht);

        % monoms from symbolic table represent one column in the matrixj
        
        if f4_debug() then <<
            prin2t {"convert_hashes_to_columns: matrix", matrixj};
            prin2t {"convert_hashes_to_columns: symbol_ht", symbol_ht}
        >>;

        col2hash := dv_undef(loadj - 1);
        j := 1;
        % number of pivotal cols
        k := 0;
        for i := hashtable_htget_offset(symbol_ht) : loadj do <<
            % column to hash index
            putv(col2hash, j, i);
            j := j + 1;

            % meaning the column is pivoted
            if hashtable_hvget_idx(getv(hdata, i)) = 2 then
                k := k + 1
        >>;

        % sort columns
        sorting_sort_columns_by_hash(col2hash, symbol_ht);
        
        if f4_debug() then <<
            prin2t {"convert_hashes_to_columns: sorted by hash", col2hash}
        >>;

        matrix_mset_nleft(matrixj, k);
        % -1 as long as hashtable loadj is always 1 more than actual
        matrix_mset_nright(matrixj, loadj - k - 1);

        % store the other direction of mapping,
        % hash -> column
        for k := 1:dv_length(col2hash) do
            hashtable_hvset_idx(getv(hdata, getv(col2hash, k)), k);
        
        nterms := 0;
        uprows := matrix_mget_uprows(matrixj);
        for k := 1:matrix_mget_nup(matrixj) do <<
            row := getv(uprows, k);

            for j := 1:dv_length(row) do
                putv(row, j, hashtable_hvget_idx(getv(hdata, getv(row, j))));
            
            nterms := nterms + dv_length(row)
        >>;

        lowrows := matrix_mget_lowrows(matrixj);
        for k := 1:matrix_mget_nlow(matrixj) do <<
            row := getv(lowrows, k);

            for j := 1:dv_length(row) do
                putv(row, j, hashtable_hvget_idx(getv(hdata, getv(row, j))));
            
            nterms := nterms + dv_length(row)
        >>;

        matrix_mset_ncols(matrixj, matrix_mget_nleft(matrixj) + matrix_mget_nright(matrixj));

        % Julia: asserts

        matrix_mset_col2hash(matrixj, col2hash)
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure matrix_convert_matrix_rows_to_basis_elements(matrixj: MacaulayMatrix, 
                            basis: Basis, ht: MonomialHashtable, symbol_ht: MonomialHashtable);
    begin scalar rows, crs, bcoeffs, bgens, mcoeffs, low2coef, colidx;
        
        basis_check_enlarge_basis(basis, matrix_mget_npivots(matrixj));

        rows := matrix_mget_lowrows(matrixj);
        crs := basis_bget_ndone(basis);

        bcoeffs := basis_bget_coeffs(basis);
        bgens := basis_bget_gens(basis);

        mcoeffs := matrix_mget_coeffs(matrixj);
        low2coef := matrix_mget_low2coef(matrixj);

        for i := 1:matrix_mget_npivots(matrixj) do <<
            colidx := getv(getv(rows, i), 1);
            % Julia: 
            hashtable_insert_in_basis_hash_table_pivots(getv(rows, i), ht, symbol_ht, matrix_mget_col2hash(matrixj));

            putv(bcoeffs, crs + i, getv(mcoeffs, getv(low2coef, colidx)));
            putv(bgens, crs + i, getv(rows, i))
        >>;

        basis_bset_ntotal(basis, basis_bget_ntotal(basis) + matrix_mget_npivots(matrixj))
    end;

%--------------------------------------------------------------------------------------------------

endmodule; % end of matrixj module

end; % of file