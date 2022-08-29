load_package groebner$
torder({}, revgradlex)$

k1 := 123/1000;
k2 := 5/8;
k3 := 3/10;
k4 := 307/500;
k5 := 6/5;
k6 := 5/2;
k7 := 1;
k8 := 1;
operator diff$
odes := { diff(x1, t) = (1*k7*k4*x6 + (-1)*k7*k2*x3)/k7,
  diff(x3, t) = (2*k7*k3*x2*x6 + (-1)*k7*k2*x3 + (-1)*k7*k1*x4)/k7,
  diff(x4, t) = (1*k7*k2*x3 + (-1)*k7*k1*x4)/k7,
  diff(x5, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file