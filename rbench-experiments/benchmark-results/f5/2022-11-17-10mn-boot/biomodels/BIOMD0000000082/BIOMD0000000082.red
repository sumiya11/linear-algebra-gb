load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 5000000;
k3 := 1/2;
k4 := 100000000;
k5 := 1/10;
k6 := 1/10;
k7 := 100000;
k8 := 10000000;
k9 := 1/10;
k10 := 1/20;
k11 := 1/10;
k12 := 1/100000000;
k13 := 1001/100000000;
k14 := 1/1000000000;
k15 := 10001/1000000000;
k16 := 1/1000000000;
k17 := 10001/1000000000;
operator diff$
odes := { diff(x1, t) = (-1)*k1*(k2*x1*x8 - k3*x2)/k1,
  diff(x2, t) = (1*k1*(k2*x1*x8 - k3*x2) + (-1)*k1*(k4*x2*x9 - k5*x3) + 1*k1*k10*x6)/k1,
  diff(x3, t) = (1*k1*(k4*x2*x9 - k5*x3) + (-1)*k1*(k6*x3 - k7*x4*x5))/k1,
  diff(x4, t) = (1*k1*(k6*x3 - k7*x4*x5) + (-1)*k1*(k8*x4*x7 - k9*x6))/k1,
  diff(x5, t) = 1*k1*(k6*x3 - k7*x4*x5)/k1,
  diff(x6, t) = (1*k1*(k8*x4*x7 - k9*x6) + (-1)*k1*k10*x6)/k1,
  diff(x7, t) = (-1)*k1*(k8*x4*x7 - k9*x6)/k1,
  diff(x8, t) = (-1)*k1*(k2*x1*x8 - k3*x2)/k1,
  diff(x9, t) = ((-1)*k1*(k4*x2*x9 - k5*x3) + 1*k1*k11*x10)/k1,
  diff(x10, t) = (1*k1*k10*x6 + (-1)*k1*k11*x10)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file