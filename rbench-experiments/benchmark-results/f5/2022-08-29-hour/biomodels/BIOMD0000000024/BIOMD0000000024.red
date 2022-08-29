load_package f5$
torder({}, revgradlex)$

k1 := 1/1000000000000000;
k2 := 1;
k3 := 1;
k4 := 2;
k5 := 1;
k6 := 3;
k7 := 4;
k8 := 21/100;
k9 := 21/100;
k10 := 0;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = (1*k1*k2/(1 + (x3/k3)^k4) + (-1)*k1*k8*x2)/k1,
  diff(x3, t) = (1*k1*k5*delay(x2, k7)^k6 + (-1)*k1*k9*x3)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file