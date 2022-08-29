load_package groebner$
torder({}, revgradlex)$

k1 := 10;
k2 := 3/100;
k3 := 1500;
k4 := 1/50;
k5 := 3/125000;
k6 := 3/1000;
k7 := 6/25;
k8 := 12/5;
k9 := 1000;
k10 := 1000;
k11 := 1;
operator diff$
odes := { diff(x1, t) = (1*k11*(k1 + k2*x1) + (-1)*k11*(k4*x1 + k5*x4*x1 + k2*x1*(x1 + x2 + x3)/k3))/k11,
  diff(x2, t) = (1*k11*k5*x4*x1 + (-1)*k11*(k4*x2 + k6*x2))/k11,
  diff(x3, t) = (1*k11*k6*x2 + (-1)*k11*k7*x3)/k11,
  diff(x4, t) = (1*k11*k9*k7*x3 + (-1)*k11*(k5*x4*x1 + k8*x4))/k11 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file