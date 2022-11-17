load_package f4$
torder({}, revgradlex)$

k1 := 1;
k2 := 5000000;
k3 := 10;
k4 := 100000000;
k5 := 1/10;
k6 := 5;
k7 := 100000;
k8 := 5000000;
k9 := 55;
k10 := 1;
k11 := 2;
k12 := 31/1000000;
k13 := 41/1000000;
k14 := 1/10000000000;
k15 := 100001/10000000000;
k16 := 1/1000000;
k17 := 11/1000000;
operator diff$
odes := { diff(x1, t) = (-1)*k1*(k2*x1*x10 - k3*x2)/k1,
  diff(x2, t) = (1*k1*(k2*x1*x10 - k3*x2) + (-1)*k1*(k4*x2*x4 - k5*x3) + 1*k1*k10*x7)/k1,
  diff(x3, t) = (1*k1*(k4*x2*x4 - k5*x3) + (-1)*k1*(k6*x3 - k7*x6*x5))/k1,
  diff(x4, t) = ((-1)*k1*(k4*x2*x4 - k5*x3) + 1*k1*k11*x9)/k1,
  diff(x5, t) = (1*k1*(k6*x3 - k7*x6*x5) + (-1)*k1*(k8*x5*x8 - k9*x7))/k1,
  diff(x6, t) = 1*k1*(k6*x3 - k7*x6*x5)/k1,
  diff(x7, t) = (1*k1*(k8*x5*x8 - k9*x7) + (-1)*k1*k10*x7)/k1,
  diff(x8, t) = (-1)*k1*(k8*x5*x8 - k9*x7)/k1,
  diff(x9, t) = (1*k1*k10*x7 + (-1)*k1*k11*x9)/k1,
  diff(x10, t) = (-1)*k1*(k2*x1*x10 - k3*x2)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file