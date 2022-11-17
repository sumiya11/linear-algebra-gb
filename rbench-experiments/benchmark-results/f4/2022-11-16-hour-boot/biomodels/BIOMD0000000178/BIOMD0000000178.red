load_package f4$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 1;
k4 := 1;
k5 := 141/1000;
k6 := 13/1000;
k7 := 29/500;
k8 := 3/20000;
k9 := 1;
k10 := 0;
k11 := 0;
operator diff$
odes := { diff(x2, t) = (-1)*k8*x2*k2/k2,
  diff(x3, t) = ((-1)*k7*x3*k2 + 1*k8*x2*k2)/k2,
  diff(x4, t) = ((-1)*k5*x4*k2 + 1*k7*x3*k2)/k2,
  diff(x5, t) = (1*k5*x4*k2 + (-1)*k6*x5*k3)/k3,
  diff(x6, t) = 1*k6*x5*k3/k4 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file