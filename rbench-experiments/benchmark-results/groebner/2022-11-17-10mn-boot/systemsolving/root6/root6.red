% root-6 system in lex
% characteristic 0
%
% from "Tropical Approach to the Cyclic n-Roots Problem"

load_package groebner;

system := {
x1 + x2 + x3 + x4 + x5 + x6, x1*x2 + x1*x3 + x1*x4 + x1*x5 + x1*x6 + x2*x3 + x2*x4 + x2*x5 + x2*x6 + x3*x4 + x3*x5 + x3*x6 + x4*x5 + x4*x6 + x5*x6, x1*x2*x3 + x1*x2*x4 + x1*x2*x5 + x1*x2*x6 + x1*x3*x4 + x1*x3*x5 + x1*x3*x6 + x1*x4*x5 + x1*x4*x6 + x1*x5*x6
+ x2*x3*x4 + x2*x3*x5 + x2*x3*x6 + x2*x4*x5 + x2*x4*x6 + x2*x5*x6 + x3*x4*x5 + x3*x4*x6 + x3*x5*x6 + x4*x5*x6, x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x3*x6 + x1*x2*x4*x5 + x1*x2*x4*x6 + x1*x2*x5*x6 + x1*x3*x4*x5 + x1*x3*x4*x6 + x1*x3*x5*x6 + x1*x4*x5*x6 + x2*x3*x4*x5 + x2*x3*x4*x6 + x2*x3*x5*x6 + x2*x4*x5*x6 + x3*x4*x5*x6, x1*x2*x3*x4*x5 + x1*x2*x3*x4*x6 + x1*x2*x3*x5*x6 + x1*x2*x4*x5*x6 + x1*x3*x4*x5*x6 + x2*x3*x4*x5*x6, x1*x2*x3*x4*x5*x6 + 1
}$

vars := {x1, x2, x3, x4, x5, x6}$
torder(vars, lex)$

gb := groebner(system)$

end;
