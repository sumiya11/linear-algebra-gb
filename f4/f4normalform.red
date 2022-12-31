module f4normalform;
% Computing normal form of several polynomials at once.
% This file corresponds to file gb/normalform.jl in Groebner.jl

revision('f4normalform, "$Id$");

copyright('f4normalform, "(c) 2023 A. Demin, T. Sturm, MPI Informatics, Germany");

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
% Compute the normal of `tobereduced` with respect to `basis`.
asserted procedure normalform_normal_form_f4(ring: PolyRing, basis: Basis, ht: MonomialHashtable, tobereduced: Basis): Basis;
    begin scalar matrixj, symbol_ht;
        matrixj := matrix_initialize_matrix(ring);
        symbol_ht := hashtable_initialize_secondary_hash_table(ht);

        f4_select_tobereduced(basis, tobereduced, matrixj, symbol_ht, ht);

        f4_symbolic_preprocessing(basis, matrixj, ht, symbol_ht);

        matrix_convert_hashes_to_columns(matrixj, symbol_ht);

        sorting_sort_matrix_rows_decreasing(matrixj);

        matrix_linear_algebra_nf(ring, matrixj, tobereduced, basis);

        matrix_convert_nf_rows_to_basis_elements(matrixj, tobereduced, ht, symbol_ht);

        return tobereduced
    end;

endmodule; % end of f4normalform

end; % eof