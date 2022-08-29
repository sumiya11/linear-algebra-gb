load_package f4$
torder({}, revgradlex)$

k1 := 409/250;
k2 := 1/500;
k3 := 3743/10000;
k4 := 1/25;
k5 := 19/50;
k6 := 11/200;
k7 := 1/100;
k8 := 1/500;
k9 := 1;
operator diff$
odes := { diff(x1, t) = (1*k9*k1*x1*(1 - k2*x1) + (-1)*k9*x1*x2)/k9,
  diff(x2, t) = (1*k9*k4*x1*x2 + (-1)*k9*k3*x2 + 1*k9*k7*x2*x3)/k9,
  diff(x3, t) = (1*k9*k5 + 1*k9*k8*x1*x3 + (-1)*k9*k6*x3)/k9 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file