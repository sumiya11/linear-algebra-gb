load_package f4$
torder({}, revgradlex)$

k1 := 11/250;
k2 := 1/5;
k3 := 199/1000;
k4 := 187/1000;
k5 := 361/500;
k6 := 101/1000;
k7 := 9/200;
k8 := 1723/1000;
k9 := 27/1000;
k10 := 911/1000;
k11 := 743/1000;
k12 := 71/250;
k13 := 1523/1000;
k14 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k14*k4*x1*x2 + 1*k14*k1*x1*(k13*x5 + 1))/k14,
  diff(x2, t) = (1*k14*(k3*x3 + k2*x4) + (-1)*k14*(k3*x3*x2 + k2*x4*x2 + k5*x2*x1))/k14,
  diff(x3, t) = (1*k14*(k6*x1 + k7*x2) + (-1)*k14*(k6*x1*x3 + k7*x2*x3 + k8*x3*x1))/k14,
  diff(x4, t) = (1*k14*(k9 + k9*x3*x4) + (-1)*k14*(k9*x3 + k9*x4 + k10*x4*x2))/k14,
  diff(x5, t) = (1*k14*(k11*k12*x5 + k11*x3*x5*x5) + (-1)*k14*(k11*k12*x5*x5 + k11*x3*x5))/k14 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file