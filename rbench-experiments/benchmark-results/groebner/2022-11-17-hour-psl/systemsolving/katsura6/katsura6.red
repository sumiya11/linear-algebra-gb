% katsura-6 system in revgradlex
% characteristic 0
% 0 dim
%
% PoSSo test suite
% https://www-sop.inria.fr/saga/POL/BASE/2.multipol/

load_package groebner;

system := {x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 - 1,
    2*x1*x6 + 2*x2*x5 + 2*x2*x7 + 2*x3*x4 - x6,
    2*x1*x5 + 2*x2*x4 + 2*x2*x6 + x3^2 + 2*x3*x7 - x5,
    2*x1*x4 + 2*x2*x3 + 2*x2*x5 + 2*x3*x6 + 2*x4*x7 - x4,
    2*x1*x3 + x2^2 + 2*x2*x4 + 2*x3*x5 - x3 + 2*x4*x6 + 2*x5*x7,
    2*x1*x2 + 2*x2*x3 - x2 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7,
    x1^2 - x1 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2}$

vars := {x1, x2, x3, x4, x5, x6, x7}$
torder(vars, revgradlex)$

gb := groebner(system)$

end;
