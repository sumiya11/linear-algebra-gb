module f4isgroebner;
% Checking if a basis is a groebner basis.
% This file corresponds to file gb/isgroebner.jl in Groebner.jl

revision('f4isgroebner, "$Id$");

copyright('f4isgroebner, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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

% This functionality is used in checking the correctness of 
% reconstructed Groebner basis.
%
% Checks that the given `basis` stores a Groebner basis.
asserted procedure isgroebner_isgroebner_f4(ring: PolyRing, basis: Basis, ht: MonomialHashtable): Boolean;
    begin scalar matrixj, symbol_ht, update_ht, pairset, plcm;
        matrixj := matrix_initialize_matrix(ring);
        symbol_ht := hashtable_initialize_secondary_hash_table(ht);
        update_ht := hashtable_initialize_secondary_hash_table(ht);

        pairset := basis_initialize_pairset();

        plcm := dv_undef(0);
        basis_update(pairset, basis, ht, update_ht, plcm);

        if basis_psget_load(pairset) = 0 then
            return t;

        f4_select_isgroebner(pairset, basis, matrixj, symbol_ht);

        f4_symbolic_preprocessing(basis, matrixj, ht, symbol_ht);

        matrix_convert_hashes_to_columns(matrixj, symbol_ht);

        sorting_sort_matrix_rows_increasing(matrixj);
        sorting_sort_matrix_rows_decreasing(matrixj);

        return matrix_linear_algebra_isgroebner(ring, matrixj, basis)
    end;

endmodule; % end of f4isgroebner

end; % eof