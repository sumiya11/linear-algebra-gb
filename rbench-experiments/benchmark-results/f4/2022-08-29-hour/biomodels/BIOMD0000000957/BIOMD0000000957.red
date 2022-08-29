load_package f4$
torder({}, revgradlex)$

k1 := 209/1000000000;
k2 := 909/1000;
k3 := 1/10;
k4 := 1;
k5 := 5999815;
operator diff$
odes := { diff(x1, t) = (-1)*k4*k1*x2*x1/k4,
  diff(x2, t) = (1*k4*k1*x2*x1 + (-1)*k4*k3*x2 + (-1)*k4*k2*x2)/k4,
  diff(x3, t) = 1*k4*k3*x2/k4,
  diff(x4, t) = 1*k4*k2*x2/k4 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file