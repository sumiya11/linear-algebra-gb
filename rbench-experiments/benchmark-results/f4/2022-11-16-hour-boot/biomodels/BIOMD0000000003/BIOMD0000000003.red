load_package f4$
torder({}, revgradlex)$

k1 := 1/17;
k2 := 1/100;
k3 := 3;
k4 := 1;
k5 := 1/2;
k6 := 1;
k7 := 1/40;
k8 := 1/100;
k9 := 1/4;
k10 := 1/50;
k11 := 1/200;
k12 := 3/2;
k13 := 1/200;
k14 := 1/200;
k15 := 1/200;
k16 := 1/2;
operator diff$
odes := { diff(x1, t) = (1*k6*k7 + (-1)*x1*k6*k8 + (-1)*x1*k6*k9*x3*(x1 + k10)^(-1))/k6,
  diff(x2, t) = (1*k6*(1 + (-1)*x2)*x1*k3*(x1 + k5)^(-1)*(k11 + (-1)*x2 + 1)^(-1) + (-1)*k6*x2*k12*(k13 + x2)^(-1))/k6,
  diff(x3, t) = (1*k6*x2*k4*(1 + (-1)*x3)*(k14 + (-1)*x3 + 1)^(-1) + (-1)*k6*k16*x3*(k15 + x3)^(-1))/k6 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file