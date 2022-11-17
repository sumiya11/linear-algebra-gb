load_package f5$
torder({}, revgradlex)$

k1 := 14/5;
k2 := 153/20;
k3 := 29/50;
k4 := 1/1000;
k5 := 39/25;
k6 := 1/4;
k7 := 47/50;
k8 := 27/100;
k9 := 1/2;
k10 := 7/200;
k11 := 1;
k12 := 17/200;
k13 := 44000;
k14 := 6;
k15 := 1;
k16 := 44000;
operator diff$
odes := { diff(x1, t) = (-1)*k15*(k1*x3*x1/k13 + k5*k1*x6*x1/k13 + k2*x4*x1/k13)/k15,
  diff(x2, t) = (1*k15*(k1*x3*x1/k13 + k5*k1*x6*x1/k13 + k2*x4*x1/k13) + (-1)*k15*k6*(1 - k3 - k4)*x2 + (-1)*k15*k6*k4*x2 + (-1)*k15*k6*k3*x2)/k15,
  diff(x3, t) = (1*k15*k6*k3*x2 + (-1)*k15*k10*x3 + (-1)*k15*k7*x3 + (-1)*k15*k8*x3)/k15,
  diff(x4, t) = (1*k15*k6*k4*x2 + (-1)*k15*k11*x4 + (-1)*k15*k7*x4 + (-1)*k15*k8*x4)/k15,
  diff(x5, t) = 1*k15*k6*(1 - k3 - k4)*x2/k15,
  diff(x6, t) = (1*k15*k7*x4 + 1*k15*k7*x3 + (-1)*k15*k9*x6 + (-1)*k15*k12*x6)/k15,
  diff(x7, t) = (1*k15*k8*x4 + 1*k15*k8*x3 + 1*k15*k9*x6)/k15,
  diff(x8, t) = (1*k15*k11*x4 + 1*k15*k10*x3 + 1*k15*k12*x6)/k15 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file