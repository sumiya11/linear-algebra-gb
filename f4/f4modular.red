module f4modular;

asserted procedure modular_gcdext(a: Integer, m: Integer): List;
    begin integer buf, u1, u2, u3, v1, v2, v3,
                    buf1, buf2, buf3;
        
        ASSERT(a < m);
        ASSERT(gcdf(a, m) = 1);

        u1 := 1;
        u2 := 0;
        u3 := m;
        v1 := 0;
        v2 := 1;
        v3 := a;

        while v3 neq 0 do <<
            buf := u3 / v3;

            buf1 := buf * v1;
            buf2 := buf * v2;
            buf3 := buf * v3;

            buf1 := u1 - buf1;
            buf2 := u2 - buf2;
            buf3 := u3 - buf3;

            u1 := v1;
            u2 := v2;
            u3 := v3;

            v1 := buf1;
            v2 := buf2;
            v3 := buf3
        >>;

        return {u3, u2, u1} 
    end;

asserted procedure modular_invmod(a: Integer, m: Integer): Integer;
    begin scalar x;
        x := cadr modular_gcdext(a, m);
        if x < 0 then
            x := x + m;
        return x
    end;

%--------------------------------------------------------------------------------------------------

asserted procedure modular_rational_reconstruction_bound(modulo: Integer): Integer;
    floor(sqrt(modulo / 2) + 1);

asserted procedure modular_rational_reconstruction(bnd: Integer, a: Integer, m: Integer): DottedPair;
    begin scalar buf, u1, u2, u3, v1, v2, v3, result,
                    buf1, buf2, buf3, den, num;
        if a = 0 then
            return t . (nil ./ 1);
        
        % assumes input is nonnegative
        ASSERT(a > 0);

        if a = 1 then
            return t . (1 ./ 1);
        
        u1 := 1;
        u2 := 0;
        u3 := m;
        v1 := 0;
        v2 := 1;
        v3 := a;

        result := t;
        while t do <<
            if v2 > bnd then <<
                result := false;
                go to Return_
            >>;

            buf := v3;
            if buf < 0 then
                buf := -buf;
            
            if buf < bnd then
                go to Return_;
            
            buf := u3 / v3;

            buf1 := buf * v1;
            buf2 := buf * v2;
            buf3 := buf * v3;

            buf1 := u1 - buf1;
            buf2 := u2 - buf2;
            buf3 := u3 - buf3;

            u1 := v1;
            u2 := v2;
            u3 := v3;

            v1 := buf1;
            v2 := buf2;
            v3 := buf3
        >>;
    Return_:
    
        den := v2;
        num := v3;

        if den < 0 then <<
            den := -den;
            num := -num
        >>;

        return result . (num ./ den)
    end;

%--------------------------------------------------------------------------------------------------

% Julia: in Julia, the function modifies its argument;
% Here we return the answer.
asserted procedure modular_CRT(M: Integer, a1: Integer, minv1: Integer, 
                                a2: Integer, minv2: Integer, m1: Integer, m2: Integer): Integer;
    begin scalar buf, n1, n2;
        buf := m1 * minv1;
        n1 := buf * a2;

        buf := m2 * minv2;
        n2 := buf * a1;

        buf := n1 + n2;
        return remainder(buf, M)
    end;

endmodule; % end of modular 

end; % eof