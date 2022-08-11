
asserted procedure io_PolyRing(nvars: Integer, explen: Integer, ord): PolyRing;
    {'pr, nvars, explen, ord};

asserted procedure io_pr_nvars(pr: PolyRing): Integer;
    cadr pr;

asserted procedure io_pr_explen(pr: PolyRing): Integer;
    caddr pr;

asserted procedure io_pr_ord(pr: PolyRing);
    cadddr pr;

%--------------------------------------------------------------------------------------------------

% list of lp --> list of Poly1
asserted procedure io_tointernal(polys)
    begin

    end;

% list of Poly1 --> list of lp
asserted procedure io_output(polys)
    begin
        
    end;

end;