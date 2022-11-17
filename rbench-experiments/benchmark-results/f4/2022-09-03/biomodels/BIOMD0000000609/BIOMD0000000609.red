load_package f4$
torder({}, revgradlex)$

k1 := 33/250000000000000;
k2 := 2;
k3 := 687/50000000000000000;
k4 := 1600000000000000000;
k5 := 299/100;
k6 := 226000000000000;
k7 := 53/2000000000000000;
k8 := 2;
k9 := 63/200;
k10 := 63/2000;
k11 := 110;
k12 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k12*k6*x4*x1 + (-1)*k12*k8*x1 + 1*k12*k7)/k12,
  diff(x2, t) = ((-1)*k12*k4*x3*x2 + (-1)*k12*k2*x2 + 1*k12*k3)/k12,
  diff(x3, t) = (1*k12*k9*x4 + (-1)*k12*k10*x3 + (-1)*k12*k11*x3 + (-1)*k12*k4*x3*x2)/k12,
  diff(x4, t) = ((-1)*k12*k9*x4 + 1*k12*k10*x3 + (-1)*k12*k5*x4 + (-1)*k12*k6*x4*x1)/k12,
  diff(x5, t) = 1*k12*k11*x3/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file