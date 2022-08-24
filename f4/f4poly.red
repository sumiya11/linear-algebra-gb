module f4poly;

fluid '(global!-dipvars!*);
fluid '(vdpsortmode!*);

procedure f5_isPolynomial(x); eqcar(x, 'p);
procedure f5_isCoeff(x); sqp(x) or integerp(x);
struct Polynomial checked by f5_isPolynomial;
struct Terms checked by listp;
struct Term checked by listp;
struct Coeffs checked by listp;
% Coeff can be either an Integer or a Standard Quotient.
struct Coeff checked by f5_isCoeff;

asserted procedure poly_initRing(vars, ord);
    <<
        global!-dipvars!* := vars;
        vdpsortmode!* := ord
    >>;

% Invariant: the first entry in the exponent list is the sum of subsequent entries

% Returns the first entry in the exponent list - the total degree
asserted inline procedure poly_totalDegExp(e1: List): Integer;
   car e1;

% Returns exponent list of zeros of the appropriate length
asserted inline procedure poly_zeroExp(): List;
   for x := 1:length(global!-dipvars!*) collect 0;

asserted procedure poly_sumExp(e1: List, e2: List): List;
   <<
      ASSERT(not null e1 and not null e2);
      if car e1 #= 0 then
         e2
      else if car e2 #= 0 then
         e1
      else
         for each x in e1 collect x #+ pop e2
   >>;

% Return the elementwise subtraction of exponent lists e1, e2
asserted procedure poly_subExp(e1: List, e2: List): List;
   <<
      ASSERT(not null e1 and not null e2);
      for each x in e1 collect x #- pop e2
   >>;

% Returns the elementwise maximum of exponent lists e1, e2
asserted procedure poly_elmaxExp(e1: List, e2: List): List;
   begin scalar w;
      ASSERT(not null e1 and not null e2);
      e1 := cdr e1;
      e2 := cdr e2;
      w := for each x in e1 collect max(x, pop e2);
      return (for each x in w sum x) . w
   end;

% Checks if e1 is elementwise less or equal than e2
asserted procedure poly_divExp!?(e1: List, e2: List): Boolean;
   begin scalar brk;
      % We deliberately also check the total degree.
      while not brk and e1 do
         brk := pop e1 #> pop e2;
      return not brk
   end;

% Checks that at least one of e1[i] or e2[i] is zero for each i
asserted procedure poly_disjExp!?(e1: List, e2: List): Boolean;
   begin scalar ok;
      e1 := cdr e1;
      e2 := cdr e2;
      ok := t;
      while ok and e1 do
         ok := pop e1 #= 0 or pop e2 #= 0;
      return ok
   end;

% Insert the "dg" into the ev in the place of variable v
asserted procedure poly_insertExp(ev: List, v: Any, dg: Integer, vars: List): List;
   (car ev #+ dg) . poly_insertExp1(cdr ev, v, dg, vars);

asserted procedure poly_insertExp1(ev: List, v: Any, dg: Integer, vars: List): List;
   if null ev or null vars then
      nil
   else if car vars eq v then
      dg . cdr ev
   else
      car ev . poly_insertExp1(cdr ev, v, dg, cdr vars);

asserted inline procedure poly_2aExp(e);
   % Returns list of prefix equivalents of exponent vector e.
   ev_2aExp1(cdr e,  cdr global!-dipvars!*);

procedure ev_2aExp1(u,v);
   if null u then
      nil
   else if car u #= 0 then
      ev_2aExp1(cdr u,cdr v)
   else if car u #= 1 then
      car v . ev_2aExp1(cdr u,cdr v)
   else
      {'expt,car v,car u} . ev_2aExp1(cdr u,cdr v);

%--------------------------------------------------------------------------------------------------

% Returns the identity Term (just one).
asserted inline procedure poly_identityTerm(): Term;
   poly_zeroExp();

asserted inline procedure poly_isIdentityTerm!?(tm: Term): Boolean;
   poly_totalDegTerm(tm) #= 0;

% Returns the total degree of the Term
asserted inline procedure poly_totalDegTerm(a: Term): Integer;
   poly_totalDegExp(a);

% Returns a * b
asserted inline procedure poly_mulTerm(a: Term, b: Term): Term;
   poly_sumExp(a, b);

% Returns a / b,
% Assuming a >= b.
asserted inline procedure poly_divTerm(a: Term, b: Term): Term;
   poly_subExp(a, b);

% Checks that a | b
asserted inline procedure poly_dividesTerm!?(a: Term, b: Term): Boolean;
   poly_divExp!?(a, b);

% Returns lcm(a, b)
asserted inline procedure poly_lcmTerm(a: Term, b: Term): Term;
   poly_elmaxExp(a, b);

% Returns a < b in the current term order term order
asserted inline procedure poly_cmpTerm(a: Term, b: Term): Boolean;
   poly_cmpExp(a, b);

% Checks if gcd(a, b) is one
asserted inline procedure poly_disjTerm!?(a: Term, b: Term): Boolean;
   poly_disjExp!?(a, b);

% Checks if a = b
asserted inline procedure poly_eqTerm!?(a: Term, b: Term): Boolean;
   poly_eqExp!?(a, b);

%--------------------------------------------------------------------------------------------------

asserted procedure poly_cmpExpLex(e1: List, e2: List): Boolean;
   <<
      e1 := cdr e1;
      e2 := cdr e2;
      while e1 and (car e1 #= car e2) do <<
         e1 := cdr e1;
         e2 := cdr e2
      >>;
      e1 and (car e1 #< car e2)
   >>;

% Compares exponent lists e1, e2 w.r.t. graded lex term order,
% and returns e1 < e2
asserted inline procedure poly_cmpExpGradLex(e1: List, e2: List): Boolean;
   (car e1 #< car e2) or (car e1 #= car e2 and poly_cmpExpLex(e1, e2));

% Compares exponent lists e1, e2 w.r.t. graded reversed lex term order,
% and returns e1 < e2
asserted inline procedure poly_cmpExpRevGradLex(e1: List, e2: List): Boolean;
   (car e1 #< car e2) or (car e1 #= car e2 and not poly_cmpExpLex(cdr e1, cdr e2));

% Compares exponent lists e1, e2 w.r.t. the current term order
asserted procedure poly_cmpExp(e1: List, e2: List): Boolean;
   <<
   if vdpsortmode!* eq 'lex then
      poly_cmpExpLex(e1, e2)
   else if vdpsortmode!* eq 'gradlex then
      poly_cmpExpGradLex(e1, e2)
   else if vdpsortmode!* eq 'revgradlex then
      poly_cmpExpRevGradLex(e1, e2)
   >>;

% Compares exponent lists w.r.t. the total degree
asserted inline procedure poly_tdegCmpExp(e1: List, e2: List): Boolean;
   car e1 #< car e2;

% Checks that e1 = e2 elementwise
asserted inline procedure poly_eqExp!?(e1: List, e2: List): Boolean;
   e1 = e2;

%--------------------------------------------------------------------------------------------------

% Constructor of Polynomial, forms a Polynomial
% from a list of `Term`s, a list of `Coeff`s, and a sugar degree.
asserted procedure poly_PolynomialWithSugar(ts: Terms, cfs: Coeffs, sugar: Integer): Polynomial;
   {'p, ts, cfs, sugar};

% Same as above, but doesn't care about sugar
asserted procedure poly_Polynomial(ts: Terms, cfs: Coeffs): Polynomial;
   poly_PolynomialWithSugar(ts, cfs, 0);

% Returns zero polynomial
asserted procedure poly_zero(): Polynomial;
   poly_Polynomial(nil, nil);

asserted procedure poly_iszero!?(p: Polynomial): Boolean;
   null poly_getTerms(p);

%--------------------------------------------------------------------------------------------------

% Returns the tail of the polynomial `poly`.
% That is, `poly - lead(poly)`,
% and preserves the sugar degree.
asserted inline procedure poly_tail(poly: Polynomial): Polynomial;
   poly_PolynomialWithSugar(poly_tailTerms(poly), poly_tailCoeffs(poly), poly_getSugar(poly));

% Returns the leading term of `poly`
asserted inline procedure poly_leadTerm(poly: Polynomial): Term;
   car poly_getTerms(poly);

% Returns the leading coefficient of `poly`
asserted inline procedure poly_leadCoeff(poly: Polynomial): Coeff;
   car poly_getCoeffs(poly);

% Returns the tail terms of `poly`
asserted inline procedure poly_tailTerms(poly: Polynomial): Terms;
   cdr poly_getTerms(poly);

% Returns the tail coefficients of `poly`
asserted inline procedure poly_tailCoeffs(poly: Polynomial): Coeffs;
   cdr poly_getCoeffs(poly);

% Returns the length of `poly`, i.e., the number of terms
asserted inline procedure poly_length(poly: Polynomial): Integer;
   length(poly_getTerms(poly));

% Returns true if the polynomial is constant
asserted inline procedure poly_isConst!?(poly: Polynomial): Boolean;
   poly_iszero!?(poly) or poly_eqExp!?(poly_leadTerm(poly), poly_zeroExp());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ADAPTIVE COEFFICIENT ARITHMETIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted inline procedure poly_iszeroCoeff!?(a: Coeff): Boolean;
   not numr(a);

asserted inline procedure poly_isoneCoeff!?(a: Coeff): Boolean;
   eqn(numr(a), 1) and eqn(denr(a), 1);

asserted inline procedure poly_zeroCoeff(): Coeff;
   poly_2Coeff(0);

asserted inline procedure poly_oneCoeff(): Coeff;
   poly_2Coeff(1);

asserted inline procedure poly_2Coeff(c: Any): Coeff;
   if sqp c then
      c
   else
      c ./ 1;

asserted inline procedure poly_2aCoeff(c: Coeff);
   prepsq c;

asserted inline procedure poly_addCoeff(a: Coeff, b: Coeff): Coeff;
   addsq(a, b);

asserted inline procedure poly_subCoeff(a: Coeff, b: Coeff): Coeff;
   subtrsq(a, b);

asserted inline procedure poly_mulCoeff(a: Coeff, b: Coeff): Coeff;
   multsq(a, b);

asserted inline procedure poly_negCoeff(a: Coeff): Coeff;
   negsq(a);

asserted inline procedure poly_isNegCoeff!?(a: Coeff): Boolean;
   minusf numr a;

asserted inline procedure poly_divCoeff(a: Coeff, b: Coeff): Coeff;
   quotsq(a, b);

asserted inline procedure poly_invCoeff(a: Coeff): Coeff;
   <<
      % ASSERT(not null numr a);
      denr(a) ./ numr(a)
   >>;

%--------------------------------------------------------------------------------------------------

% Returns s = fmult*f - fcoeff*gmult*g
%
% !! Assuming the leading terms of `fmult*f` and `fcoeff*gmult*g`
%    mutually cancel each other
asserted inline procedure poly_paircombTail(f: Polynomial, fmult: Term, fcoeff: Coeff, 
                                             g: Polynomial, gmult: Term, gcoeff: Coeff): Polynomial;
   % Assuming the leading monomials of `fmult*f` and `fcoeff*gmult*g`
   % are mutually canceled, we can use the general reduction applied to
   % the tails of input polynomials
   poly_paircomb(poly_tail(f), fmult, fcoeff, poly_tail(g), gmult, gcoeff);

% Returns s = gcoeff*fmult*f - fcoeff*gmult*g
asserted procedure poly_paircomb(f: Polynomial, fmult: Term, fcoeff: Coeff, 
                                 g: Polynomial, gmult: Term, gcoeff: Coeff): Polynomial;
   begin scalar fterms, fcoeffs, gterms, gcoeffs, gmultcoeff, isOneFmult,
                sterms, scoeffs, ft, gt, fc, gc, newc, sugar, isOneGmult,
                isOneFmultCf, isOneGmultCf, fmultcoeff;
      % We return s, a new polynomial,
      % constructed as s = fmultcoeff*fmult*f - gmultcoeff*gmult*g.
      % We form two lists, sterms and scoeffs, which would be the list of
      % terms and the list of coefficients of s.
      % The sterms list is formed by merging two sorted lists:
      % the list of terms of f each multiplied by fmult,
      % with the list of terms of g each multiplied by gmult.
      % In parallel (in the same loop), scoeffs list is formed
      % by merging the list of coefficients of f each multiplied by fmultcoeff,
      % with the list of coefficients of g each multiplied by gmultcoeff,
      % in the same merge order as for the terms.
      fterms  := poly_getTerms(f); fcoeffs := poly_getCoeffs(f);
      gterms  := poly_getTerms(g); gcoeffs := poly_getCoeffs(g);
      gmultcoeff := fcoeff;
      fmultcoeff := gcoeff;
      % identity check is 1 car and 1 comparison
      isOneFmult := poly_isIdentityTerm!?(fmult);
      isOneGmult := poly_isIdentityTerm!?(gmult);
      % identity check is 1 comparison
      isOneFmultCf := poly_isoneCoeff!?(fmultcoeff);
      isOneGmultCf := poly_isoneCoeff!?(gmultcoeff);
      % Merge two sorted lists: fterms and gterms, multiplied by fmult and gmult, respectively.
      % Merge in the same order two other lists:
      % fcoeffs and gcoeffs, multiplied by gmultcoeff and fmultcoeff, respectively.
      %
      % if n = length(fterms), m = length(gterms), then in the worst case
      %   2(n+m)*E  +    m*E     +     2m*C    +   (n+m)
      % comparison    term mult.    coeff op.     reverse
      %
      % where E is the cost of iterating exponent vector,
      % and C is the cost of one arithmetic operation on coefficients
      while fterms and gterms do <<
         if null ft then <<
            ft := car fterms;
            if not isOneFmult then
               ft := poly_mulTerm(ft, fmult)
         >>;
         if null gt then <<
            gt := car gterms;
            if not isOneGmult then
               gt := poly_mulTerm(gt, gmult)
         >>;
         % Optimization: return -1,0,1 just as C comparator;
         if poly_cmpTerm(gt, ft) then <<   % if term gt < term ft
            push(ft, sterms);
            if isOneFmultCf then
               push(car fcoeffs, scoeffs)
            else
               push(poly_mulCoeff(fmultcoeff, car fcoeffs), scoeffs);
            pop(fterms); pop(fcoeffs);
            ft := nil
         >> else if poly_eqTerm!?(gt, ft) then <<  % if term gt = term ft
            fc := if isOneFmultCf then
               car fcoeffs
            else
               poly_mulCoeff(car fcoeffs, fmultcoeff);
            gc := if isOneGmultCf then
               car gcoeffs
            else
               poly_mulCoeff(car gcoeffs, gmultcoeff);
            newc := poly_subCoeff(fc, gc);
            if not poly_iszeroCoeff!?(newc) then <<
               push(gt, sterms);
               push(newc, scoeffs)
            >>;
            pop(fterms); pop(fcoeffs);
            pop(gterms); pop(gcoeffs);
            gt := nil;
            ft := nil
         >> else <<   % if term gt > term ft
            push(gt, sterms);
            if isOneGmultCf then
               push(poly_negCoeff(car gcoeffs), scoeffs)
            else
               push(poly_negCoeff(poly_mulCoeff(car gcoeffs, gmultcoeff)), scoeffs);
            pop(gterms); pop(gcoeffs);
            gt := nil
         >>
      >>;
      if null gterms and null fterms then <<
         scoeffs := reversip(scoeffs);
         sterms  := reversip(sterms)
      >>;
      % Merge what is left from gterms and gcoeffs
      if not null gterms then <<
         if poly_isIdentityTerm!?(gmult) then
            sterms := nconc(reversip sterms, gterms)
         else <<
            while gterms do
               push(poly_mulTerm(pop(gterms), gmult), sterms);
            sterms := reversip(sterms)
         >>;
         while gcoeffs do
            push(poly_negCoeff(poly_mulCoeff(pop(gcoeffs), gmultcoeff)), scoeffs);
         scoeffs := reversip(scoeffs)
      >>;
      % Merge what is left from fterms and fcoeffs
      if not null fterms then <<
         if isOneFmult then
            sterms := nconc(reversip sterms, fterms)
         else <<
            while fterms do
               push(poly_mulTerm(pop(fterms), fmult), sterms);
            sterms := reversip(sterms)
         >>;
         while fcoeffs do
            push(poly_mulCoeff(pop(fcoeffs), fmultcoeff), scoeffs);
         scoeffs := reversip(scoeffs)
      >>;
      sugar := max(poly_getSugar(f) + poly_totalDegTerm(fmult),
                        poly_getSugar(g) + poly_totalDegTerm(gmult));
      return poly_PolynomialWithSugar(sterms, scoeffs, sugar)
   end;

% Returns poly*c
asserted procedure poly_multCoeffs(poly: Polynomial, cf: Coeff): Polynomial;
   if poly_isoneCoeff!?(cf) then
      poly
   else
      poly_PolynomialWithSugar(
         poly_getTerms(poly), 
         for each c in poly_getCoeffs(poly) collect poly_mulCoeff(c, cf), 
         poly_getSugar(poly)
      );

%--------------------------------------------------------------------------------------------------

asserted inline procedure poly_sumPoly(f: Polynomial, g: Polynomial): Polynomial;
   if poly_iszero!?(f) then
      g
   else if poly_iszero!?(g) then
      f
   else 
      poly_paircomb(f, poly_identityTerm(), poly_negCoeff(poly_oneCoeff()),
                     g, poly_identityTerm(), poly_oneCoeff());

%--------------------------------------------------------------------------------------------------

% Constructs a Polynomial from a SF `u`
asserted inline procedure poly_sf2poly(u: SF): Polynomial;
   poly_sq2poly(u ./ 1);

% Conversion to Polynomial: scan the standard form. ev and bc are the
% exponent and coefficient parts collected so far from higher parts.
asserted procedure poly_sf2poly1(u: SF, ev: List, bc: Coeff): Polynomial;
  if null u then
     poly_zero()
  else if domainp u then
     poly_PolynomialWithSugar({ev}, {poly_mulCoeff(bc, poly_2Coeff(u))}, poly_totalDegExp(ev))
  else
     poly_sumPoly(poly_sf2poly2(mvar u,ldeg u,lc u,ev,bc), poly_sf2poly1(red u,ev,bc));

% Returns true if `var` is the variable in the current polynomial ring,
% nil otherwise
asserted procedure poly_isPolyVar(var: Any): Boolean;
   begin scalar vars;
      vars := global!-dipvars!*;
      while vars and not (car vars eq var) do
         vars := cdr vars;
      return not null vars
   end;

% returns var^deg as a Standard Quotient
asserted inline procedure poly_groundCoeff(var, deg);
   ((((var . deg) . 1) . nil) ./ 1);

% Conversion to Polynomial: multiply leading power either into exponent vector
% or into coefficient
asserted procedure poly_sf2poly2(var,dg,c,ev,bc): Polynomial;
   if poly_isPolyVar(var) then 
      poly_sf2poly1(c,poly_insertExp(ev, var, dg, cdr global!-dipvars!*), bc)
   else
      poly_sf2poly1(c, ev, poly_mulCoeff(bc, poly_groundCoeff(var, dg)));

% Constructs a Polynomial from a SQ `u`.
% Assuming the Polynomial is correct, `poly_poly2sq` is the inverse function, 
% so that poly_poly2sq(poly_sq2poly(x)) = x.
asserted inline procedure poly_sq2poly(u: SQ): Polynomial;
   poly_sq2poly1(u, poly_zeroExp(), poly_oneCoeff());

% Recusrively converts a Standard Quotient X/Y to a Polynomial.
% Polynomial variables are encoded in exponent vector `ev`, 
% and parameters are encoded in coefficient `bc`.
% If polynomial variable is encountered in denominator Y, an error is raised.
asserted procedure poly_sq2poly1(u: SQ, ev: List, bc: Coeff): Polynomial;
   begin scalar numpoly;
      for each var in kernels denr u do <<
         if poly_isPolyVar(var) then
            rederr {"Polynomial variable in denominator in input:", var}
      >>;
      numpoly := poly_sf2poly1(numr u, ev, bc);
      return poly_multCoeffs(numpoly, 1 ./ denr u)
   end;

% SQ --> Poly

% Converts cf*tm where cf is a Coeff and tm is a Term
% to a Standard Quotient
asserted procedure poly_lead2sq(cf: Coeff, tm: Term): SQ;
   begin scalar vs;
      vs := cdr global!-dipvars!*;
      for each e in cdr tm do    % not very general
         cf := multsq(mksp!*(numr simp (pop vs), e) ./ 1, cf);
      return cf
   end;

% Converts a Polynomial to a Standard Quotient
asserted procedure poly_poly2sq(p: Polynomial): SQ;
   addsq(
      poly_lead2sq(poly_leadCoeff(p), poly_leadTerm(p)), 
      poly_poly2sq(poly_tail(p))
   );

% Poly --> Lisp Prefix 

% Returns prefix equivalent to the sum of elements of u
asserted procedure poly_replus(u: List): Any;
   if atom u then u else if null cdr u then car u else 'plus . u;

% Returns prefix equivalent to the product of elements of u.
% u is a list of prefix expressions the first of which is a number.
procedure poly_retimes(u: List): Any;
   if car u = 1 then
      if cdr u then poly_retimes cdr u else 1
   else if null cdr u then
      car u
   else
      'times . u;

% Returns prefix equivalent to the Polynomial u.
asserted procedure poly_poly2lp1(u: Polynomial): Any;
   begin scalar x,y;
      if poly_iszero!?(u) then
         return nil;
      x := poly_leadCoeff u;
      y := poly_2aExp poly_leadTerm u;
      if poly_isNegCoeff!?(x) then <<
         return {'minus,poly_retimes(poly_2aCoeff(poly_negCoeff(x)) . y)} . poly_poly2lp1 poly_tail(u)
      >>;
      return poly_retimes(poly_2aCoeff x . y) . poly_poly2lp1 poly_tail(u)
   end;

% Converts a Polynomial to an equivalent Lisp Prefix
asserted procedure poly_poly2lp(f: Polynomial): Any;
   if poly_iszero!?(f) then 0 else poly_replus poly_poly2lp1(f);

endmodule; % end of poly

end;