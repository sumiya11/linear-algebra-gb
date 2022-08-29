load_package groebner$
torder({}, revgradlex)$

k1 := 7/1000;
k2 := 33/100;
k3 := 21/5000;
k4 := 1/100;
k5 := 1;
operator diff$
odes := { diff(x1, t) = (-1)*k5*k1*x1*x2/k5,
  diff(x2, t) = (1*k5*k2 + (-1)*k5*k3*x1 + (-1)*k5*k4*x2)/k5 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file