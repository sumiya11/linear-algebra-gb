load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 9/20;
k3 := 500000000;
k4 := 10;
k5 := 100;
k6 := 0;
k7 := 10;
k8 := 199/100;
operator diff$
odes := { diff(x1, t) = (1*k2 + (-1)*k4*x1*(1 + x1)*(1 + x2)^2/(k3 + (1 + x1)^2*(1 + x2)^2))/k1,
  diff(x2, t) = (50*k4*x1*(1 + x1)*(1 + x2)^2/(k3 + (1 + x1)^2*(1 + x2)^2) + (-1)*k7*x2*(1 + k6*x2)*(1 + x3)^2/(k5 + (1 + k6*x2)^2*(1 + x3)^2))/k1,
  diff(x3, t) = (1/50*k7*x2*(1 + k6*x2)*(1 + x3)^2/(k5 + (1 + k6*x2)^2*(1 + x3)^2) + (-1)*k8*x3)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file