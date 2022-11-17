load_package f4$
torder({}, revgradlex)$

k1 := 9/100;
k2 := 127/1000;
k3 := 99/1000;
k4 := 9/100;
k5 := 93/1000;
k6 := 49/500;
k7 := 49/500;
k8 := 111/1000;
k9 := 12/125;
k10 := 34/5;
k11 := 0;
k12 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k12*k1*x1 + 1*k12*k10)/k12,
  diff(x2, t) = (1*k12*k1*x1 + (-1)*k12*k2*x2)/k12,
  diff(x3, t) = (1*k12*k2*x2 + (-1)*k12*k3*x3)/k12,
  diff(x4, t) = (1*k12*k3*x3 + (-1)*k12*k4*x4)/k12,
  diff(x5, t) = (1*k12*k4*x4 + (-1)*k12*k5*x5)/k12,
  diff(x6, t) = (1*k12*k5*x5 + (-1)*k12*k6*x6)/k12,
  diff(x7, t) = (1*k12*k6*x6 + (-1)*k12*k7*x7)/k12,
  diff(x8, t) = (1*k12*k7*x7 + (-1)*k12*k8*x8)/k12,
  diff(x9, t) = ((-1)*k12*k9*x9 + 1*k12*k8*x8)/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file