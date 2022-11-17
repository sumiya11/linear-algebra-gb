load_package groebner$
torder({}, revgradlex)$

k1 := 9/10;
k2 := 1;
k3 := 0;
k4 := 11/10;
k5 := 4/5;
k6 := 4/5;
k7 := 1/10000;
k8 := 17/10;
k9 := 1;
operator diff$
odes := { diff(x1, t) = (1*k9*k1*k2 + (-1)*k9*k3*x1 + (-1)*k9*k8*x2*x1/(x1 + k7))/k9,
  diff(x2, t) = (1*k9*k6*x3 + (-1)*k9*k5*x2)/k9,
  diff(x3, t) = (1*k9*k4*x1*k2 + (-1)*k9*k6*x3)/k9 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file