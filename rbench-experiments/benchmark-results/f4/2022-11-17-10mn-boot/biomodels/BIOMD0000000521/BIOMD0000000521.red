load_package f4$
torder({}, revgradlex)$

k1 := 4833/100;
k2 := 713/100;
k3 := 206/5;
k4 := 121/1000;
k5 := 59/20000;
k6 := 31/10000;
k7 := 87/10000;
k8 := 729/1000;
k9 := 6/25;
k10 := 100;
k11 := 1;
k12 := 1;
operator diff$
odes := { diff(x1, t) = (-k9)*x1,
  diff(x2, t) = k4*x2*(1 - (x2 + x3 + x4)/k10) + k6*x4 - k5*x2 - k8*x1*k9*x2,
  diff(x3, t) = k5 - k8*x1*k9*x3,
  diff(x4, t) = k8*x1*k9*x3 - k6*x4 - k7*x4 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file