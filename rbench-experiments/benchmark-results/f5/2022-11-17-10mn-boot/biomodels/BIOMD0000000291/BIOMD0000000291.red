load_package f5$
torder({}, revgradlex)$

k1 := 121/250;
k2 := 3979/50000;
k3 := 1019/200000000;
k4 := 83/3125000;
k5 := 5489/1000000;
k6 := 1613/5000000000;
k7 := 301/100000;
k8 := 1011/10000000000;
k9 := 337/20000;
k10 := 53/400;
k11 := 95000;
k12 := 3000;
k13 := 1;
k14 := 1;
operator diff$
odes := { diff(x8, t) = 0/k14,
  diff(x9, t) = 0/k14,
  diff(x10, t) = 0/k14,
  diff(x1, t) = k11*k3*x5*x6 - k12*k4*x1*x6 - k3*x1 - k9*x1*x7^(k13 + 1) + k4*x4 + k10*x4*x7,
  diff(x2, t) = k6*x7*x6 - k8*x2 + k9*x1*x7^(k13 + 1) + k10*x4*x7,
  diff(x3, t) = k5*x7^k13*x5 - k7*x3 + k9*x1*x7^(k13 + 1),
  diff(x4, t) = k12*k4*x1*x6 - k4*x4 - k10*x4*x7 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file