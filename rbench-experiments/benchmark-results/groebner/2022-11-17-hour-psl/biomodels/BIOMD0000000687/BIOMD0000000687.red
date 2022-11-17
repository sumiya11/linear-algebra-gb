load_package groebner$
torder({}, revgradlex)$

k1 := 175/11;
k2 := 1/10;
k3 := 1/5;
k4 := 1;
k5 := 1;
k6 := 10;
k7 := 1/10;
k8 := 1/10;
k9 := 1/2;
k10 := 1/5;
k11 := 1/10;
k12 := 1/100;
k13 := 1/1000000;
k14 := 31/2;
k15 := 1;
k16 := 1/10;
k17 := 1;
operator diff$
odes := { diff(x1, t) = (1*k17*k6 + (-1)*k17*k7*x1 + (-1)*k17*k8*x5*x1 + (-1)*k17*k9*x5*x1)/k17,
  diff(x2, t) = (1*k17*k8*x5*x1 + (-1)*k17*k2*x2 + (-1)*k17*k12*x2 + 1*k17*k11*x4 + (-1)*k17*k13*x6*x2)/k17,
  diff(x3, t) = (1*k17*k12*x2 + (-1)*k17*k3*x3 + (-1)*k17*k13*x6*x3)/k17,
  diff(x4, t) = (1*k17*k9*x5*x1 + (-1)*k17*k11*x4 + (-1)*k17*k7*x4)/k17,
  diff(x5, t) = (1*k17*k4*x3 + (-1)*k17*k5*x5)/k17,
  diff(x6, t) = (1*k17*k10*x15 + 1*k17*k14*(x2 + x3)*x6 + (-1)*k17*k16*x6)/k17,
  diff(x7, t) = (-1)*k17*k15*x7/k17,
  diff(x8, t) = (2*k17*k15*x7 + (-1)*k17*k15*x8)/k17,
  diff(x9, t) = (2*k17*k15*x8 + (-1)*k17*k15*x9)/k17,
  diff(x10, t) = (2*k17*k15*x9 + (-1)*k17*k15*x10)/k17,
  diff(x11, t) = (2*k17*k15*x10 + (-1)*k17*k15*x11)/k17,
  diff(x12, t) = (2*k17*k15*x11 + (-1)*k17*k15*x12)/k17,
  diff(x13, t) = (2*k17*k15*x12 + (-1)*k17*k15*x13)/k17,
  diff(x14, t) = (2*k17*k15*x13 + (-1)*k17*k15*x14)/k17,
  diff(x15, t) = (2*k17*k15*x14 + (-1)*k17*k15*x15)/k17 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file