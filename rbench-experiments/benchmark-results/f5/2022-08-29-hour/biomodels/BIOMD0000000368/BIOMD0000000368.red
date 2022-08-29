load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 1/10;
k3 := 1/10;
k4 := 1/10;
k5 := 1/10;
k6 := 1;
k7 := 1;
k8 := 5;
k9 := 5;
k10 := 1;
k11 := 1;
operator diff$
odes := { diff(x1, t) = (-(k1*x6 + k10*x8))*x1 + k6*x5,
  diff(x2, t) = (-k2)*x5*x2 + k7*x6,
  diff(x3, t) = (-(k3*x6 + k4*x8))*x3 + k8*x7,
  diff(x4, t) = (-k5)*x7*x4 + k9*x8,
  diff(x5, t) = (k1*x6 + k10*x8)*x1 - k6*x5,
  diff(x6, t) = k2*x5*x2 - k7*x6,
  diff(x7, t) = (k3*x6 + k4*x8)*x3 - k8*x7,
  diff(x8, t) = k5*x7*x4 - k9*x8 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file