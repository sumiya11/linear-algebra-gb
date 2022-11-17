load_package f5$
torder({}, revgradlex)$

k1 := 0;
k2 := 1;
k3 := 1/5;
k4 := 1/20;
k5 := 19/50000;
k6 := 19/50000;
k7 := 1/200;
k8 := 3/12500;
k9 := 5/2;
k10 := 10;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = (1*k2*k3*k9 + (-1)*k2*k6*x2 + (-1)*k2*(k7*x2*x4 - k8*x5))/k2,
  diff(x3, t) = (1*k2*k4*x5 + (-1)*k2*k5*x3)/k2,
  diff(x4, t) = (1*k2*k4*x5 + (-1)*k2*(k7*x2*x4 - k8*x5))/k2,
  diff(x5, t) = ((-1)*k2*k4*x5 + 1*k2*(k7*x2*x4 - k8*x5))/k2 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file