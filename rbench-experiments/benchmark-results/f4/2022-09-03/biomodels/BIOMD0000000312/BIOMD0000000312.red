load_package f4$
torder({}, revgradlex)$

k1 := 2;
k2 := 2;
k3 := 1;
k4 := 1;
k5 := 4;
k6 := 1;
operator diff$
odes := { diff(x1, t) = (1*k6*k1*x3 + (-1)*k6*k2*x1*x2)/k6,
  diff(x2, t) = (1*k6*k3*x3 + (-1)*k6*k4*x2)/k6 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file