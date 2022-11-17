load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 199/100;
k3 := 36000000;
k4 := 1/100;
k5 := 9/40000000;
k6 := 1;
k7 := 29/10;
k8 := 6;
k9 := 0;
k10 := 1/2;
k11 := 4/5;
k12 := 1;
operator diff$
odes := { diff(x1, t) = (1*k12*(k1 + k2*x1*(1 - (x1 + 1)/k3)) + (-1)*k12*(k4*x1 + (1 - k10*k11)*k5*x3*x1))/k12,
  diff(x2, t) = (1*k12*(1 - k10*k11)*k5*x3*x1 + (-1)*k12*k6*x2)/k12,
  diff(x3, t) = (1*k12*(1 - (k9 + k11)/2)*k7*x2 + (-1)*k12*k8*x3)/k12,
  diff(x4, t) = (1*k12*(k9 + k11)/2*k7*x2 + (-1)*k12*k8*x4)/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file