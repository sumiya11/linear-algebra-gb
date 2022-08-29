load_package groebner$
torder({}, revgradlex)$

k1 := 1/50;
k2 := 1/50;
k3 := 1/5;
k4 := 1/2;
k5 := 4/5;
k6 := 1/100;
k7 := 4/5;
k8 := 1/125;
k9 := 16/3;
k10 := 100/7;
k11 := 1;
operator diff$
odes := { diff(x1, t) = (1*k11*k1 + (-1)*k11*k2*x1 + (-1)*k11*k3*x1*x2 + (-1)*k11*k4*x1*x3)/k11,
  diff(x2, t) = (1*k11*k5*k3*x1*x2 + (-1)*k11*(k2 + k6)*x2)/k11,
  diff(x3, t) = (1*k11*k7*k4*x1*x3 + (-1)*k11*(k2 + k8)*x3)/k11 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file