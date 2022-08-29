load_package f5$
torder({}, revgradlex)$

k1 := 1/125;
k2 := 1/1000;
k3 := 1/1000;
k4 := 1/1000;
k5 := 100000000;
k6 := 1000000000;
k7 := 9/100;
k8 := 9/100;
k9 := 1/50;
k10 := 1/50;
k11 := 9/100;
k12 := 1/20;
k13 := 3/100;
k14 := 3/100;
k15 := 1/50;
k16 := 1/50;
k17 := 1;
operator diff$
odes := { diff(x1, t) = (1*k17*k1*x3 + 1*k17*k7*x1*x3*(1 - (x1 + x2)/k5) + (-1)*k17*k13*x1)/k17,
  diff(x2, t) = (1*k17*k2*x4 + 1*k17*k8*x2*x4*(1 - (x2 + x1)/k5) + (-1)*k17*k14*x2)/k17,
  diff(x3, t) = (1*k17*k3*x1 + 1*k17*k9*x3*(1 - (x3 + x4)/k6) + (-1)*k17*k15*x3 + (-1)*k17*k12*x3 + 1*k17*k11*x4)/k17,
  diff(x4, t) = (1*k17*k12*x3 + (-1)*k17*k11*x4 + 1*k17*k4*x2 + 1*k17*k10*x4*x2*(1 - (x4 + x3)/k6) + (-1)*k17*k16*x4)/k17 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file