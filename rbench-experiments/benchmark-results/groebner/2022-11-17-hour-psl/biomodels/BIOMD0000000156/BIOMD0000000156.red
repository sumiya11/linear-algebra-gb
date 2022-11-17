load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 0;
k3 := 37/10;
k4 := 3/2;
k5 := 9/10;
k6 := 11/10;
k7 := 2;
k8 := 1;
operator diff$
odes := { diff(x1, t) = (1*k8*k7*x1*k1 + (-1)*k8*k3*x2*x1)/k8,
  diff(x2, t) = (1*k8*k6*x3 + (-1)*k8*k5*x2)/k8,
  diff(x3, t) = (1*k8*k4*x1*k1 + (-1)*k8*k6*x3)/k8 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file