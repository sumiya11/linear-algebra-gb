load_package f5$
torder({}, revgradlex)$

k1 := 6;
k2 := 18;
k3 := 1/25;
k4 := 1/50;
k5 := 7/50;
k6 := 6/125;
k7 := 1;
k8 := 11081000;
k9 := 11081000;
operator diff$
odes := { diff(x1, t) = (-1)*k7*(k1*k3*x3*x1 + k2*k4*x2*x1)/k8/k7,
  diff(x2, t) = (1*k7*(k1*k3*x3*x1 + k2*k4*x2*x1)/k8 + (-1)*k7*k5*x2)/k7,
  diff(x3, t) = (1*k7*k5*x2 + (-1)*k7*k6*x3)/k7,
  diff(x4, t) = 1*k7*k6*x3/k7,
  diff(x5, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file