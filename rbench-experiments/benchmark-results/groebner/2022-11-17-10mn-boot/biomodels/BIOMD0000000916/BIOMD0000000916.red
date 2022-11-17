load_package groebner$
torder({}, revgradlex)$

k1 := 18/5;
k2 := 6/5;
k3 := 133/5;
k4 := 6/5;
k5 := 6/5;
k6 := 1;
k7 := 1;
operator diff$
odes := { diff(x1, t) = (-1)*k6*k3*x1/k6,
  diff(x2, t) = (1*k6*k3*x1 + (-1)*k6*(k4*x2 - k5*x4) + (-1)*k6*k1*x2)/k6,
  diff(x3, t) = 1*k6*k2*x4/k6,
  diff(x4, t) = (1*k6*(k4*x2 - k5*x4) + (-1)*k6*k2*x4)/k6,
  diff(x5, t) = 1*k6*k1*x2/k6 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file