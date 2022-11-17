load_package f5$
torder({}, revgradlex)$

k1 := 20;
k2 := 30;
k3 := 361/1000;
k4 := 1/20;
k5 := 1/5000;
k6 := 67/500;
k7 := 171/500;
k8 := 1/400;
k9 := 501/10000;
k10 := 1;
operator diff$
odes := { diff(x1, t) = 1*k10*(k1*k2 - x1*(k3 + k4))/k10,
  diff(x2, t) = 1*k10*(x1*k3 - x2*(k6 + k9) - k5*x2^2)/k10,
  diff(x3, t) = 1*k10*(x2*k6 - x3*(k7 + k8))/k10 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file