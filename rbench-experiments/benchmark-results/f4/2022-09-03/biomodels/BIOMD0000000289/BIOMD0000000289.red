load_package f4$
torder({}, revgradlex)$

k1 := 1/400;
k2 := 1/10000;
k3 := 2000;
k4 := 200;
k5 := 2/125;
k6 := 1000;
k7 := 1/4;
k8 := 1/4;
k9 := 1/4;
k10 := 5;
k11 := 1/4;
k12 := 3/1000000;
k13 := 1/2;
k14 := 2001/400;
k15 := 1600/2001;
k16 := 1;
k17 := 0;
operator diff$
odes := { diff(x1, t) = (1*k2*k1*x4 + (-1)*k7*x1 + (-1)*k11*x1 + (-1)*k12*x1*x2)/k16,
  diff(x2, t) = (1*k4*x1 + 1*k5*x3*x1 + (-1)*k8*x2)/k16,
  diff(x3, t) = (1*k6*x1 + (-1)*k9*x3)/k16,
  diff(x4, t) = ((-1)*k1*x4 + 1*k3*x3 + (-1)*k10*x4)/k16,
  diff(x5, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file