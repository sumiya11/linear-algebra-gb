module f4coeffs;
% Polynomial coefficients manipulations

% Julia: no CoeffBuffer

asserted procedure CoeffAccum(gb_coeffs_zz: Vector, prev_gb_coeffs_zz: Vector, gb_coeffs_qq: Vector): CoeffAccum;
    begin scalar v;
        v := dv_undef(3);
        putv(v, 1, gb_coeffs_zz);
        putv(v, 2, prev_gb_coeffs_zz);
        putv(v, 3, gb_coeffs_qq);
        return v;
    end;



endmodule; % end of f4coeffs

end; % eof