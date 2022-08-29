load_package groebner$
torder({}, revgradlex)$

k1 := 1/50;
k2 := 1;
k3 := 1/100;
k4 := 4/125;
k5 := 1;
k6 := 15;
k7 := 9/200;
k8 := 1;
k9 := 23/250;
k10 := 1;
k11 := 1/100;
k12 := 1/100;
k13 := 1;
k14 := 1/2;
k15 := 43/500;
k16 := 11/10000;
k17 := 1;
k18 := 500;
k19 := 50;
k20 := 100;
operator diff$
odes := { diff(x1, t) = ((-1)*k17*(k1*x1*x4 - k2*x6) + 1*k17*(k15*x11 - k16*x1*x5))/k17,
  diff(x2, t) = (1*k17*k3*x6 + (-1)*k17*(k4*x2*x4 - k5*x7) + 1*(k10*x9 - k11*x2*x5) + (-1)*k17*(k12*x2*x5 - k13*x10))/k17,
  diff(x3, t) = (1*k17*k6*x7 + (-1)*k17*(k7*x3*x5 - k8*x8))/k17,
  diff(x4, t) = ((-1)*k17*(k1*x1*x4 - k2*x6) + 1*k17*k3*x6 + (-1)*k17*(k4*x2*x4 - k5*x7) + 1*k17*k6*x7)/k17,
  diff(x5, t) = ((-1)*k17*(k7*x3*x5 - k8*x8) + 1*(k10*x9 - k11*x2*x5) + (-1)*k17*(k12*x2*x5 - k13*x10) + 1*k17*(k15*x11 - k16*x1*x5))/k17,
  diff(x6, t) = (1*k17*(k1*x1*x4 - k2*x6) + (-1)*k17*k3*x6)/k17,
  diff(x7, t) = (1*k17*(k4*x2*x4 - k5*x7) + (-1)*k17*k6*x7)/k17,
  diff(x8, t) = (1*k17*(k7*x3*x5 - k8*x8) + (-1)*k17*k9*x8)/k17,
  diff(x9, t) = (1*k17*k9*x8 + (-1)*(k10*x9 - k11*x2*x5))/k17,
  diff(x10, t) = (1*k17*(k12*x2*x5 - k13*x10) + (-1)*k17*k14*x10)/k17,
  diff(x11, t) = (1*k17*k14*x10 + (-1)*k17*(k15*x11 - k16*x1*x5))/k17 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file