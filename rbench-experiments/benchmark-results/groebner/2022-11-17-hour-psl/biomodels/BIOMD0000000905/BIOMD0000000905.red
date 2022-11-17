load_package groebner$
torder({}, revgradlex)$

k1 := 1/2;
k2 := 1/20;
k3 := 1/5;
k4 := 1/25;
k5 := 4/5;
k6 := 3/10;
k7 := 1/10;
k8 := 12/5;
k9 := 1/10;
k10 := 5;
k11 := 2/5;
k12 := 6/5;
k13 := 1/50;
k14 := 3/5;
k15 := 1/2;
k16 := 1;
operator diff$
odes := { diff(x1, t) = (1*k16*k1*x1 + (-1)*k16*(k2*x2*x1 + k3*x1*x1))/k16,
  diff(x2, t) = (1*k16*(k4 + k6*x2*x1) + (-1)*k16*(k5*x2 + k7*k2*x2*x1 + k15*k14*x5*x2))/k16,
  diff(x3, t) = (1*k16*k8*x1 + (-1)*k16*k9*x3)/k16,
  diff(x4, t) = (1*k16*k10 + (-1)*k16*k11*x4)/k16,
  diff(x5, t) = (1*k16*k12*x4 + (-1)*k16*(k13*x5 + k14*x5*x2))/k16 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file