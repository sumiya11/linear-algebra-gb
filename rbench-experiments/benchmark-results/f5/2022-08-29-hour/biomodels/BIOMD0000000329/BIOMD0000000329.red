load_package f5$
torder({}, revgradlex)$

k1 := 1/10000000000000;
k2 := 53/250;
k3 := 29/10;
k4 := 38/25;
k5 := 19/100;
k6 := 122/25;
k7 := 59/50;
k8 := 31/25;
k9 := 2909/100;
k10 := 806/25;
k11 := 679/50;
k12 := 4/25;
k13 := 153;
operator diff$
odes := { diff(x1, t) = (1*k1*k2 + 1*k1*k3*x1 + (-1)*k1*k4*x2*x1/(k5 + x1) + (-1)*k1*k6*x3*x1/(k7 + x1))/k1,
  diff(x2, t) = (1*k1*k8*x1 + (-1)*k1*k10*x2/(k9 + x2))/k1,
  diff(x3, t) = (1*k1*k11*x1 + (-1)*k1*k13*x3/(k12 + x3))/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file