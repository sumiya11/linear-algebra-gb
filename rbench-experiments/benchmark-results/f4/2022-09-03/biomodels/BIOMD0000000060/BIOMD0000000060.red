load_package f4$
torder({}, revgradlex)$

k1 := 0;
k2 := 1;
k3 := 144/5;
k4 := 1500;
k5 := 9/10;
k6 := 4;
k7 := 1500;
k8 := 9/10;
k9 := 3;
k10 := 3859/10;
k11 := 7/4;
k12 := 1/10;
k13 := 1;
operator diff$
odes := { diff(x1, t) = 1*(k3*x3 - k4*k5^k6*x1)/k2,
  diff(x2, t) = 1*(k7*k8^k9*x3 - k10*x2)/k2,
  diff(x3, t) = ((-1)*(k3*x3 - k4*k5^k6*x1) + (-1)*(k7*k8^k9*x3 - k10*x2) + (-1)*(k11*x3 - k12*x4))/k2,
  diff(x4, t) = 1*(k11*x3 - k12*x4)/k2 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file