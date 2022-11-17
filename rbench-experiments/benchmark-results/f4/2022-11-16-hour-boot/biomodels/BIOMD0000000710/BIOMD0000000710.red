load_package f4$
torder({}, revgradlex)$

k1 := 3/1000;
k2 := 7;
k3 := 1/2;
k4 := 3/1000000;
k5 := 1/1000;
k6 := 3/1000000;
k7 := 1/100;
k8 := 1/100;
k9 := 1/100;
k10 := 1/100;
k11 := 4;
k12 := 1/25;
k13 := 26/5;
k14 := 5000000;
k15 := 32000;
k16 := 500000000;
k17 := 800000;
k18 := 1;
operator diff$
odes := { diff(x1, t) = (1*k18*k16*k8 + (-1)*k18*k1*x1*x5 + (-1)*k18*k2*x1*x6 + (-1)*k18*k8*x1)/k18,
  diff(x2, t) = (1*k18*k1*x1*x5 + (-1)*k18*k3*x2 + (-1)*k18*k4*x2*x7)/k18,
  diff(x3, t) = (1*k18*k3*x2 + (-1)*k18*k9*x3 + (-1)*k18*k4*x3*x7)/k18,
  diff(x4, t) = (1*k18*k2*x1*x6 + (-1)*k18*k10*x4)/k18,
  diff(x5, t) = (1*k18*k7*x3 + (-1)*k18*k13*x5)/k18,
  diff(x6, t) = (1*k18*k6*x3 + (-1)*k18*k11*x6)/k18,
  diff(x7, t) = (1*k18*k17*k12 + 1*k18*k5*x3 + (-1)*k18*k12*x7)/k18 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file