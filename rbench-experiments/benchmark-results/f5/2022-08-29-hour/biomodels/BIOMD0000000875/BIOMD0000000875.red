load_package f5$
torder({}, revgradlex)$

k1 := 343/10000000000;
k2 := 1/2;
k3 := 480;
k4 := 3;
k5 := 1/2;
k6 := 134000;
k7 := 10;
k8 := 3/100;
k9 := 1;
operator diff$
odes := { diff(x1, t) = (1*k9*k7 + (-1)*k9*k8*x1 + (-1)*k9*k1*x3*x1)/k9,
  diff(x2, t) = (1*k9*k1*x3*x1 + (-1)*k9*k2*x2)/k9,
  diff(x3, t) = (1*k9*(1 - k5)*k3*k2*x2 + (-1)*k9*k4*x3)/k9,
  diff(x4, t) = (1*k9*k5*k3*k2*x2 + (-1)*k9*k4*x4)/k9 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file