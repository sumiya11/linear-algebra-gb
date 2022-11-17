load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 41/100;
k3 := 5/4;
k4 := 57/200;
k5 := 11/10;
k6 := 17/200;
k7 := 3/1000;
k8 := 13/25;
k9 := 31/2000;
k10 := 3/25;
k11 := 19/10;
k12 := 41/100;
k13 := 1;
operator diff$
odes := { diff(x1, t) = (1*k13*k11 + (-1)*k13*k4*x4*x1 + (-1)*k13*(x1*k3*x2 + x1))/k13,
  diff(x2, t) = (1*k13*(x2*k6*x1 + k8*x3) + (-1)*k13*(x2*k2/k1 + x2*k7*x3))/k13,
  diff(x3, t) = (1*k13*k4*x4*x1 + (-1)*k13*k5*x2*x3)/k13,
  diff(x4, t) = (1*k13*x4*k10*(1 - k9*x4) + (-1)*k13*k4*x4*x1)/k13 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file