load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 22649/50000000;
k3 := 130837/10000000000;
k4 := 499767/500000;
k5 := 1;
k6 := 1;
k7 := 1;
k8 := 1;
k9 := 1;
k10 := 1;
k11 := 49991/200000;
k12 := 241033/1000000;
k13 := 46949/250000;
k14 := 1;
k15 := 1;
k16 := 1;
k17 := 1;
k18 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*k2*x1*x2 + (-1)*k1*k3*x1*x3 + (-1)*k1*k4*x1*x2*x3 + 1*k1*k5*x7 + (-1)*k1*k8*x1 + 1*k1*k11*x4 + 1*k1*k12*x5 + 1*k1*k13*x6)/k1,
  diff(x2, t) = ((-1)*k1*k2*x1*x2 + (-1)*k1*k4*x1*x2*x3 + 1*k1*k6*x8 + (-1)*k1*k9*x2 + 1*k1*k11*x4 + 1*k1*k13*x6)/k1,
  diff(x3, t) = ((-1)*k1*k3*x1*x3 + (-1)*k1*k4*x1*x2*x3 + 1*k1*k7*x9 + (-1)*k1*k10*x3 + 1*k1*k12*x5 + 1*k1*k13*x6)/k1,
  diff(x4, t) = (1*k1*k2*x1*x2 + (-1)*k1*k11*x4 + (-1)*k1*k14*x4)/k1,
  diff(x5, t) = (1*k1*k3*x1*x3 + (-1)*k1*k12*x5 + (-1)*k1*k15*x5)/k1,
  diff(x6, t) = (1*k1*k4*x1*x2*x3 + (-1)*k1*k13*x6 + (-1)*k1*k16*x6)/k1,
  diff(x7, t) = 0/k1,
  diff(x8, t) = 0/k1,
  diff(x9, t) = 0/k1,
  diff(x10, t) = (1*k1*k17*x1 + (-1)*k1*k18*x10)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file