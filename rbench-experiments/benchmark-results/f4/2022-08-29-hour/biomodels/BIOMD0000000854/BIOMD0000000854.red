load_package f4$
torder({}, revgradlex)$

k1 := 11/20;
k2 := 7/20;
k3 := 7/500;
k4 := 23/1000;
k5 := 11/5;
k6 := 0;
k7 := 1;
k8 := 0;
k9 := 0;
k10 := 1;
k11 := 1;
operator diff$
odes := { diff(x1, t) = (1*k10*k2*x2 + (-1)*k10*k3*k1*x1 + 1*k10*k1*x3)/k10,
  diff(x2, t) = ((-1)*k10*k2*x2 + 1*k10*k1*x4 + (-1)*k10*k3*k1*x2)/k10,
  diff(x3, t) = (1*k10*k3*k1*x1 + (-1)*k10*k1*x3 + 1*k10*k2*x4 + (-1)*k10*k5*k2*x3)/k10,
  diff(x4, t) = ((-1)*k10*k1*x4 + 1*k10*k3*k1*x2 + (-1)*k10*k2*x4 + 1*k10*k5*k2*x3)/k10 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file