load_package groebner$
torder({}, revgradlex)$

k1 := 419/500000000000;
k2 := 3/500;
k3 := 59413/10000;
k4 := 7/20;
k5 := 63/25000000000;
k6 := 2244;
k7 := 7/20;
k8 := 19200000000;
k9 := 8730000000;
k10 := 3110;
k11 := 139/10;
k12 := 1;
k13 := 1;
k14 := 1;
operator diff$
odes := { diff(x1, t) = (-1)*k1*x4*x1*x5*x3/k12,
  diff(x2, t) = (1*k12*k2*x1*x4 + (-1)*k12*(k4 + k3)*x2 + (-1)*k5*x4*x1*x2*x6)/k12,
  diff(x3, t) = (1*k14*k6 + (-1)*k14*k7*x3 + (-1)*k8*x3*x7)/k14,
  diff(x4, t) = (-1)*k1*x4*x1*x2*x6/k12,
  diff(x5, t) = (1*k13*k10 + (-1)*k13*(k3 + k4)*x5 + (-1)*k5*x4*x1*x5*x3)/k13,
  diff(x6, t) = (1*k13*k6 + (-1)*k13*k7*x6 + (-1)*k8*x6*x7)/k13,
  diff(x7, t) = ((-1)*k12*k11*x7 + 1*k3*(x5*x4 + x2*x1) + (-1)*k9*(x6*x4 + x3*x1)*x7)/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file