load_package f5$
torder({}, revgradlex)$

k1 := 700;
k2 := 800;
k3 := 850;
k4 := 20;
k5 := 1/5;
k6 := 700;
k7 := 800;
k8 := 850;
k9 := 5;
k10 := 25;
k11 := 8/5;
k12 := 1/50;
k13 := 11/100;
k14 := 1/40;
k15 := 8/125;
k16 := 11/250;
k17 := 1;
k18 := 1/2;
k19 := 14/125;
k20 := 36/125;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = (-1)*k17*(k1*k18*x2 - k2*x3)/k17,
  diff(x3, t) = (1*k17*(k1*k18*x2 - k2*x3) + (-1)*k3*k17*x3)/k17,
  diff(x4, t) = (1*k3*k17*x3 + (-1)*k17*(k4*x4 - k5*x5))/k17,
  diff(x5, t) = 1*k17*(k4*x4 - k5*x5)/k17,
  diff(x6, t) = (-1)*k17*(k6*k18*x6 - k7*x7)/k17,
  diff(x7, t) = (1*k17*(k6*k18*x6 - k7*x7) + (-1)*k8*k17*x7)/k17,
  diff(x8, t) = (1*k8*k17*x7 + (-1)*k17*(k9*k18*x8 - k10*x9))/k17,
  diff(x9, t) = (1*k17*(k9*k18*x8 - k10*x9) + (-1)*k17*(k11*x9 - k12*x10))/k17,
  diff(x10, t) = 1*k17*(k11*x9 - k12*x10)/k17 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file