load_package groebner$
torder({}, revgradlex)$

k1 := 9/25;
k2 := 11/100;
k3 := 12/25;
k4 := 4/25;
k5 := 1/5;
k6 := 3/5;
k7 := 9/250;
k8 := 9;
k9 := 1;
operator diff$
odes := { diff(x1, t) = (1*k9*k1*x1 + (-1)*k9*k1*x1*x1 + (-1)*k9*k1*x1*x2 + (-1)*k9*k2*x1*x3)/k9,
  diff(x2, t) = ((-1)*k9*k3*x2*x4 + (-1)*k9*x2 + 1*k9*k2*x1*x3)/k9,
  diff(x3, t) = (1*k9*k8*x2 + (-1)*k9*k5*x3 + (-1)*k9*k4*x3*x4 + (-1)*k9*k2*x1*x3)/k9,
  diff(x4, t) = (1*k9*k6*x2*x4 + (-1)*k9*k7*x4)/k9 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file