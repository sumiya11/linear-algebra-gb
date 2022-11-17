load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 1/250;
k3 := 1000;
k4 := 21/100000;
k5 := 27/50000;
k6 := 3/125000;
k7 := 3/125000;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*k2*x1 + (-1)*k1*(k3*x2*x1 - k4*x4))/k1,
  diff(x2, t) = (1*k1*k2*x1 + (-1)*k1*(k3*x2*x1 - k4*x4) + 2*k1*k5*x4)/k1,
  diff(x3, t) = (1*k1*k2*x1 + 1*k1*k5*x4)/k1,
  diff(x4, t) = (1*k1*(k3*x2*x1 - k4*x4) + (-1)*k1*k5*x4)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file