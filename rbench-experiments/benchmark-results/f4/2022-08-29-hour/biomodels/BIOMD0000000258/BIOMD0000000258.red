load_package f4$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 11/10;
k4 := 1/100;
k5 := 1/100;
k6 := 1/100;
k7 := 1/100;
k8 := 1;
k9 := 1;
k10 := 1;
k11 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k8*x1/k4/(1 + x1/k4 + x2/k6) + 1*k2*k8/k3*x2/k5/(1 + x3/k7 + x2/k5))/k10,
  diff(x2, t) = (1*k8*x1/k4/(1 + x1/k4 + x2/k6) + (-1)*k2*k8/k3*x2/k5/(1 + x3/k7 + x2/k5) + (-1)*k1*k8*x2/k6/(1 + x1/k4 + x2/k6) + 1*k8/k3*x3/k7/(1 + x3/k7 + x2/k5))/k10,
  diff(x3, t) = (1*k1*k8*x2/k6/(1 + x1/k4 + x2/k6) + (-1)*k8/k3*x3/k7/(1 + x3/k7 + x2/k5))/k10 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file