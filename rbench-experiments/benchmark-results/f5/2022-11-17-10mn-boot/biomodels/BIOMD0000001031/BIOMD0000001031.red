load_package f5$
torder({}, revgradlex)$

k1 := 0;
k2 := 3/20;
k3 := 1/50;
k4 := 1;
k5 := 1;
operator diff$
odes := { diff(x1, t) = (1*k4*k1*x1 + (-1)*k4*x1*x3)/k4,
  diff(x2, t) = (1*k4*x1*x3 + (-1)*k4*k5*x2)/k4,
  diff(x3, t) = (1*k4*k5*x2 + (-1)*k4*(k3*x1*x3 + k2*x3))/k4 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file