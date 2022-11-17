load_package f4$
torder({}, revgradlex)$

k1 := 10;
k2 := 9/10;
k3 := 3/10;
k4 := 4/5;
k5 := 1/10;
k6 := 1/50;
k7 := 4/5;
k8 := 7/10;
k9 := 3/100;
k10 := 1;
operator diff$
odes := { diff(x1, t) = (1*k10*(k1 + k2*x1*(1 - x1/k4)) + (-1)*k10*k3*x1*x2)/k10,
  diff(x2, t) = (1*k10*k5*x2*x3 + (-1)*k10*k6*x2)/k10,
  diff(x3, t) = (1*k10*k7*x3*(1 - x3/k8) + (-1)*k10*(k5*x2*x3 + k9*x3))/k10 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file