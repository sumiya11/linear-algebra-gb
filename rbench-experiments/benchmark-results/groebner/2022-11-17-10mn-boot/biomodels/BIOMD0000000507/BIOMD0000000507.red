load_package groebner$
torder({}, revgradlex)$

k1 := 625/4;
k2 := 78/5;
k3 := 5/2;
k4 := 1;
k5 := 14809/500000000;
k6 := 4003/2000;
k7 := 0;
k8 := 14809/500000000;
k9 := 4003/2000;
k10 := 1;
k11 := 1;
k12 := 1;
operator diff$
odes := { diff(x1, t) = (1*k10*k1/(1 + x2^k3) + (-1)*k10*k11*x1)/k10,
  diff(x2, t) = (1*k10*k2/(1 + (x1/(1 + x3/k8)^k9)^k4) + (-1)*k10*k12*x2)/k10,
  diff(x3, t) = 0/k10 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file