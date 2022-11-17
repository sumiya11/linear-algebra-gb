load_package f4$
torder({}, revgradlex)$

k1 := 1;
k2 := 1/200;
k3 := 1/100;
k4 := 1/100000;
k5 := 1/40000;
k6 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*k2*x1 + (-1)*k1*k4*x1)/k1,
  diff(x2, t) = (1*k1*k2*x1 + (-1)*k1*k3*x2)/k1,
  diff(x3, t) = (1*k1*k3*x2 + 1*k1*k5*x4)/k1,
  diff(x4, t) = (1*k1*k4*x1 + (-1)*k1*k5*x4)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file