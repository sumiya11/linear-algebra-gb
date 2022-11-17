load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 8;
k3 := 1;
k4 := 1;
k5 := 3/2;
k6 := 1;
k7 := 1;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = 0,
  diff(x3, t) = (2*k2*k6*x4 + (-1)*k3*x3^2 + (-1)*k4*x3*x4 + (-1)*k5*x3)/k1,
  diff(x4, t) = ((-1)*k2*k6*x4 + 1*k3*x3^2)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file