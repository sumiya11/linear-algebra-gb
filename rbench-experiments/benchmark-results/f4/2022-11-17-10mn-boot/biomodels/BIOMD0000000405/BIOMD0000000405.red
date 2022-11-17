load_package f4$
torder({}, revgradlex)$

k1 := 500;
k2 := 500;
k3 := 10;
k4 := 1000;
k5 := 3465735902799/100000000000000;
k6 := 1;
k7 := 100;
operator diff$
odes := { diff(x1, t) = ((-1)*k6*k4*x1*x5 + 1*k6*k1 + (-1)*k6*k5*x1)/k6,
  diff(x2, t) = ((-1)*k6*k4*x2*x5 + 1*k6*k2 + (-1)*k6*k5*x2)/k6,
  diff(x3, t) = (1*k6*k4*x1*x5 + (-1)*k6*k3*x3)/k6,
  diff(x4, t) = (1*k6*k4*x2*x5 + (-1)*k6*k3*x4)/k6,
  diff(x5, t) = ((-1)*k6*k4*x1*x5 + (-1)*k6*k4*x2*x5 + 1*k6*k3*x3 + 1*k6*k3*x4)/k6,
  diff(x6, t) = 0/k6 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file