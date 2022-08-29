load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 5;
k3 := 770;
k4 := 420;
k5 := 770;
k6 := 5;
k7 := 27/500;
k8 := 1/50;
k9 := 11/25;
k10 := 0;
k11 := 6;
k12 := 1/50;
k13 := 0;
k14 := 0;
k15 := 20;
k16 := 0;
k17 := 1;
k18 := 170;
k19 := 12/5;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*(k2*x2*x1 - k3*x3) + 1*k1*(k5*x4 - k6*x5*x1) + (-1)*k1*(k9*x1*x7 - k10*x8))/k1,
  diff(x2, t) = ((-1)*k1*(k2*x2*x1 - k3*x3) + 1*k1*(k15*x3*x7 - k16*x2*x9))/k1,
  diff(x3, t) = (1*k1*(k2*x2*x1 - k3*x3) + (-1)*k1*k4*x3 + (-1)*k1*(k15*x3*x7 - k16*x2*x9))/k1,
  diff(x4, t) = (1*k1*k4*x3 + (-1)*k1*(k5*x4 - k6*x5*x1) + (-1)*k1*(k11*x4*x6 - k12*x9))/k1,
  diff(x5, t) = (1*k1*(k5*x4 - k6*x5*x1) + (-1)*k1*(k7*x5*x6 - k8*x7))/k1,
  diff(x6, t) = ((-1)*k1*(k7*x5*x6 - k8*x7) + (-1)*k1*(k11*x4*x6 - k12*x9))/k1,
  diff(x7, t) = (1*k1*(k7*x5*x6 - k8*x7) + (-1)*k1*(k9*x1*x7 - k10*x8) + (-1)*k1*(k15*x3*x7 - k16*x2*x9))/k1,
  diff(x8, t) = (1*k1*(k9*x1*x7 - k10*x8) + 1*k1*(k13*x9 - k14*x8))/k1,
  diff(x9, t) = (1*k1*(k11*x4*x6 - k12*x9) + (-1)*k1*(k13*x9 - k14*x8) + 1*k1*(k15*x3*x7 - k16*x2*x9))/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file