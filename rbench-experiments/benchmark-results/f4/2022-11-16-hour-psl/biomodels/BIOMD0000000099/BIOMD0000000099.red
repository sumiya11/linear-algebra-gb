load_package f4$
torder({}, revgradlex)$

k1 := 7/5;
k2 := 9/10;
k3 := 5/2;
k4 := 3/2;
k5 := 3/5;
k6 := 4/5;
k7 := 2;
k8 := 13/10;
k9 := 29/100;
k10 := 1;
k11 := 3/5;
k12 := 31/10;
k13 := 33;
k14 := 9/2;
k15 := 1;
k16 := 1;
operator diff$
odes := { diff(x1, t) = (1*k15*k11*x5 + (-1)*k15*k12*x1)/k15,
  diff(x2, t) = (1*k16*k9*x5 + (-1)*k16*k10*x2*x4)/k16,
  diff(x3, t) = (1*k16*k3*x2 + (-1)*k16*k4*x3)/k16,
  diff(x4, t) = (1*k16*k7 + (-1)*k16*k8*x4*x7)/k16,
  diff(x5, t) = (1*k16*k1*x7 + (-1)*k16*k2*x5)/k16,
  diff(x6, t) = (1*k16*k13*x1 + (-1)*k16*k14*x6*x3)/k16,
  diff(x7, t) = (1*k16*k5*x6 + (-1)*k16*k6*x7*x3)/k16 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file