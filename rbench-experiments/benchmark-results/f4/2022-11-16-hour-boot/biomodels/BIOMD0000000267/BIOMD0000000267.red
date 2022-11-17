load_package f4$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 1;
k4 := 1;
k5 := 141/1000;
k6 := 13/1000;
k7 := 29/500;
k8 := 1;
k9 := 0;
k10 := 0;
operator diff$
odes := { diff(x1, t) = (-1)*k7*x1*k2/k2,
  diff(x2, t) = ((-1)*k5*x2*k2 + 1*k7*x1*k2)/k2,
  diff(x3, t) = (1*k5*x2*k2 + (-1)*k6*x3*k3)/k3,
  diff(x4, t) = 1*k6*x3*k3/k4 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file