load_package f4$
torder({}, revgradlex)$

k1 := 2;
k2 := 9/10;
k3 := 5/2;
k4 := 3/2;
k5 := 3/5;
k6 := 4/5;
k7 := 1;
k8 := 13/10;
k9 := 3/10;
k10 := 4/5;
k11 := 7/10;
k12 := 49/10;
k13 := 23;
k14 := 9/2;
k15 := 1;
operator diff$
odes := { diff(x1, t) = (1*k1*x2 + (-1)*k2*x1*x3)/k15,
  diff(x2, t) = (1*k13*x7 + (-1)*k14*x2)/k15,
  diff(x3, t) = (1*k3*x4 + (-1)*k4*x3)/k15,
  diff(x4, t) = (1*k9*x1 + (-1)*k10*x6*x4)/k15,
  diff(x5, t) = (1*k5*x2 + (-1)*k6*x3*x5)/k15,
  diff(x6, t) = (1*k7 + (-1)*k8*x5*x6)/k15,
  diff(x7, t) = (1*k11*x1 + (-1)*k12*x7)/k15 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file