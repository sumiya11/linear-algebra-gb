load_package groebner$
torder({}, revgradlex)$

k1 := 602213670000000000000000;
k2 := 1;
k3 := 1;
k4 := 1;
k5 := 1/10;
k6 := 200;
k7 := 1/2;
k8 := 1/10;
k9 := 1/10;
k10 := 10;
k11 := 3/100;
k12 := 1/20;
k13 := 200;
k14 := 0;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = (1*k3*k4/(1 + (x3*(1 - 2/(1 + (1 + 8*k6*x3)^(1/2)))/(2*k5))^2) + (-1)*k8*x2*k3)/k3,
  diff(x3, t) = (1*k7*x2*k3 + (-1)*k9*x3*k3 + (-1)*k3*(k10*x3*2/(1 + (1 + 8*k13*x3)^(1/2)) + k11*x3)/(k12 + x3))/k3 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file