load_package f5$
torder({}, revgradlex)$

k1 := 167/500;
k2 := 411/500000;
k3 := 79/100;
k4 := 449/100;
k5 := 7/40;
k6 := 23/20;
k7 := 243/500;
k8 := 79/25000;
k9 := 37/2;
k10 := 61/2;
k11 := 1;
k12 := 923/50;
operator diff$
odes := { diff(x1, t) = ((-1)*k2*x1*x2 + 1*k1*x3 + 1*k10 + (-1)*k3*x1)/k11,
  diff(x2, t) = ((-1)*k2*x1*x2 + 1*k1*x3)/k11,
  diff(x3, t) = (1*k2*x1*x2 + (-1)*k1*x3 + (-1)*k6*x3*x5 + 1*k4*x4 + 1*k5*x4)/k11,
  diff(x4, t) = (1*k6*x3*x5 + (-1)*k4*x4 + (-1)*k5*x4)/k11,
  diff(x5, t) = ((-1)*k6*x3*x5 + 1*k4*x4 + 1*k7 + (-1)*k8*x5)/k11 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file