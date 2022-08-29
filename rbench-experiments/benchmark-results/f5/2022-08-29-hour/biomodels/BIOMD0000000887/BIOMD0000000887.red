load_package f5$
torder({}, revgradlex)$

k1 := 10;
k2 := 1/1000;
k3 := 11/1000;
k4 := 3/1000;
k5 := 29/1000;
k6 := 9/250;
k7 := 3/250;
k8 := 3/100;
k9 := 3/100;
k10 := 3/100;
k11 := 1;
operator diff$
odes := { diff(x1, t) = (1*k11*k1 + (-1)*k11*k2*x1*x3 + (-1)*k11*k7*x1)/k11,
  diff(x2, t) = (1*k11*k2*x1*x3 + 1*k11*k3*x3 + (-1)*k11*k4*x2 + (-1)*k11*k8*x2)/k11,
  diff(x3, t) = (1*k11*k4*x2 + (-1)*k11*k5*x3*x4 + (-1)*k11*k9*x3)/k11,
  diff(x4, t) = (1*k11*k6*x3 + (-1)*k11*k10*x4)/k11 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file