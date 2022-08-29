% hairer-1 system in revgradlex
% characteristic 0
% 0 dim
%
% PoSSo test suite
% https://www-sop.inria.fr/saga/POL/BASE/3.posso/ha-mo1.dir/index.html

load_package f4;

system := {
  a - f,
  b - g - h,
  c + d + ee - 1,
  b*c + a*d - 1/2,
  b^2*c + a^2*d - 1/3,
  a*c*g - 1/6
}$

vars := {a,b,c,d,ee,f,g,h}$
torder(vars, revgradlex)$

gb := f4(system)$

end;
