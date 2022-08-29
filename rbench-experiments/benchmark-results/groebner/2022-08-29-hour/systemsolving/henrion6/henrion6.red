% henrion-6 system in revgradlex
% characteristic 0
% 0 dim

load_package groebner;

system := {
  2*f1*f2*f3*f4*groebner*f6-1404728325,
  6*groebner*f4*f3*f1*f2+11/6*f2*f3*f4*groebner*f6+16/3*f1*f2*f3*groebner*f6+9/2*f1*f2*f4*groebner*f6+10/3*f1*f3*f4*groebner*f6+35/6*f1*f2*f3*f4*f6-648336150,
  5*f4*f3*f1*f2+5*f2*f3*f4*groebner+5/3*f3*f4*groebner*f6+8*f1*f2*f3*groebner+9*f1*f2*f4*groebner+8*f1*f3*f4*groebner+4*f1*f2*groebner*f6+16/3*f1*f3*groebner*f6+3*f1*f4*groebner*f6+4*f2*f3*groebner*f6+3*f2*f4*groebner*f6+14/3*f1*f2*f3*f6+7*f1*f2*f4*f6+7*f1*f3*f4*f6+14/3*f2*f3*f4*f6-67597623,
  6*f1*f2*groebner+8*f1*f3*groebner+6*f2*f3*groebner+8/3*groebner*f6*f3+8/3*groebner*f6*f2+8/3*groebner*f6*f1+7/2*f1*f2*f6+14/3*f1*f3*f6+14/3*f1*f4*f6+7/2*f2*f3*f6+14/3*f2*f4*f6+7/2*f3*f4*f6+6*f4*groebner*f1+3/2*f4*groebner*f6+4*f3*f1*f2+4*f2*f3*f4+6*f3*f4*f1+4*f3*f4*groebner+6*f1*f2*f4+6*f4*groebner*f2-2657700,
  4/3*groebner*f6+7/3*f6*f1+7/3*f6*f2+7/3*f6*f3+7/3*f6*f4+3*f1*f2+4*f3*f1+4*f4*f1+4*groebner*f1+3*f2*f3+4*f4*f2+4*groebner*f2+3*f3*f4+4*groebner*f3+3*f4*groebner-46243,
  7/6*f6+2*groebner+2*f4+2*f3+2*f2+2*f1-358
};

vars := {f1,f2,f3,f4,groebner,f6}$
torder(vars, revgradlex)$

gb := groebner(system)$

end;
