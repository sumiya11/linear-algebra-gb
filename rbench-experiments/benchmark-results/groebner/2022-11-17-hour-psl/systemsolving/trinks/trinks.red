% trinks system in lex
% characteristic 0

load_package groebner;

system := {35*p + 40*z + 25*t - 27*s,
          45*p + 35*s - 165*b - 36*h,
          -11*s*b + 3b^2 + 99*w*h,
          25*p*s - 165*b^2 + 15*w*h + 30*z*h - 18*t*h,
          15*p*t + 20*z*s - 9*w*h,
          -11*b^3 + w*p*h + 2*z*t*h}$

vars := {w,p,z,t,s,b,h}$

gb := groebner(system, vars, 'revgradlex)$
torder(vars, lex)$

gb := groebner(system)$

end;
