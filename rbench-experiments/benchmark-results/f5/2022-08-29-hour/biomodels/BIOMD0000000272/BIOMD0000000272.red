load_package f5$
torder({}, revgradlex)$

k1 := 164683/5000000;
k2 := 76;
k3 := 114701/50000000000;
k4 := 339973/50000000;
k5 := 1101/100000;
k6 := 317871/100000000;
k7 := 82021/5000000;
k8 := 999293/1000;
k9 := 0;
k10 := 1;
k11 := 1;
k12 := 1;
k13 := 999293/1000;
k14 := 0;
k15 := 0;
operator diff$
odes := { diff(x1, t) = (1*k1*k2*k12 + (-1)*k1*x1*k12 + (-1)*k3*x2*x1*k12 + 1*k4*x3*k12)/k11,
  diff(x2, t) = ((-1)*k3*x2*x1*k12 + 1*k4*x3*k12 + 1*k5*x4*k12)/k10,
  diff(x3, t) = (1*k3*x2*x1*k12 + (-1)*k4*x3*k12 + (-1)*k1*x3*k12)/k11,
  diff(x4, t) = (1*k1*x3*k12 + (-1)*k5*x4*k12 + (-1)*k6*x4*k12 + (-1)*k7*x4*k12)/k12,
  diff(x5, t) = 1*k6*x4*k12/k12,
  diff(x6, t) = 1*k7*x4*k12/k10 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file