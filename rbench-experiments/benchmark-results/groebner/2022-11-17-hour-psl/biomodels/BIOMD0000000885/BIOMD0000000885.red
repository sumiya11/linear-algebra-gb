load_package groebner$
torder({}, revgradlex)$

k1 := 1/10;
k2 := 1/100;
k3 := 1/100;
k4 := 1000000;
k5 := 1/10;
k6 := 100;
k7 := 1000;
k8 := 1721/250;
k9 := 431/1000;
k10 := 980000000;
k11 := 30218000;
k12 := 1;
operator diff$
odes := { diff(x1, t) = (1*k12*k1*k5*x1*(1 - x1/k4) + (-1)*k12*k2*(1 - k5)*x1)/k12,
  diff(x2, t) = (1*k12*k2*(1 - k5)*x1 + (-1)*k12*k3*x2)/k12,
  diff(x3, t) = (1*k12*(k6*x1 + k7*x2) + (-1)*k12*k8*x3)/k12,
  diff(x4, t) = (1*k12*k9*x4*(1 - x4/k10) + (-1)*k12*k11*x3*x4)/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file