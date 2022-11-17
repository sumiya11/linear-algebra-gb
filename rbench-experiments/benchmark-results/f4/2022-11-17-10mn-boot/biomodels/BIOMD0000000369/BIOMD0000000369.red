load_package f4$
torder({}, revgradlex)$

k1 := 1;
k2 := 1/10;
k3 := 1;
k4 := 1;
k5 := 1;
k6 := 1;
k7 := 1;
k8 := 1;
k9 := 0;
k10 := 1/1000;
k11 := 1;
operator diff$
odes := { diff(x1, t) = (-(k1*x5 + k9*x7))*x1,
  diff(x2, t) = (-k2)*(1 + k10)*x4*x2,
  diff(x3, t) = (-k4)*x6*x3,
  diff(x4, t) = (k1*x5 + k9*x7)*x1 - k5*x4,
  diff(x5, t) = k2*x4*x2 - k3*x7*x5 - k6*x5,
  diff(x6, t) = k2*k10*x4*x2 + k3*x7*x5 - k7*x6,
  diff(x7, t) = k4*x6*x3 - k8*x7 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file