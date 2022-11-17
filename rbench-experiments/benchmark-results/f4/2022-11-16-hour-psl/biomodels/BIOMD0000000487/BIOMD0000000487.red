load_package f4$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 1/10;
k4 := 1;
k5 := 1;
k6 := 1/10;
k7 := 1;
k8 := 0;
k9 := 0;
k10 := 0;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*x1*x2 + 1*k2*x3 + 1*k3*x3)/k7,
  diff(x2, t) = ((-1)*k1*x1*x2 + 1*k2*x3 + 1*k6*x6)/k7,
  diff(x3, t) = (1*k1*x1*x2 + (-1)*k2*x3 + (-1)*k3*x3)/k7,
  diff(x4, t) = ((-1)*k4*x4*x5 + 1*k5*x6 + 1*k6*x6)/k7,
  diff(x5, t) = (1*k3*x3 + (-1)*k4*x4*x5 + 1*k5*x6)/k7,
  diff(x6, t) = (1*k4*x4*x5 + (-1)*k5*x6 + (-1)*k6*x6)/k7 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file