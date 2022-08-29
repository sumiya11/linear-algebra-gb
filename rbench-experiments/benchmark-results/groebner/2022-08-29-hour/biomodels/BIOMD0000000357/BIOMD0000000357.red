load_package groebner$
torder({}, revgradlex)$

k1 := 459/5;
k2 := 412/5;
k3 := 303/2;
k4 := 2099/10;
k5 := 129/25;
k6 := 323/10;
k7 := 47/10;
k8 := 213/5;
k9 := 167/5;
k10 := 37/200;
k11 := 109/5;
k12 := 133/5000000;
k13 := 1;
k14 := 3/20000;
k15 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k13*(k1*x1*x3 - k9*x2) + 1*k13*k2*x2 + (-1)*k13*(k3*x1*x4 - k10*x5) + 1*k13*k4*x5 + (-1)*k13*(k5*x1*x3 - k11*x7) + 1*k13*k6*x7 + (-1)*k13*(k7*x1*x8 - k12*x9) + 1*k13*k8*x9)/k13,
  diff(x2, t) = (1*k13*(k1*x1*x3 - k9*x2) + (-1)*k13*k2*x2)/k13,
  diff(x3, t) = ((-1)*k13*(k1*x1*x3 - k9*x2) + (-1)*k13*(k5*x1*x3 - k11*x7))/k13,
  diff(x4, t) = (1*k13*k2*x2 + (-1)*k13*(k3*x1*x4 - k10*x5))/k13,
  diff(x5, t) = (1*k13*(k3*x1*x4 - k10*x5) + (-1)*k13*k4*x5)/k13,
  diff(x6, t) = (1*k13*k4*x5 + 1*k13*k8*x9)/k13,
  diff(x7, t) = (1*k13*(k5*x1*x3 - k11*x7) + (-1)*k13*k6*x7)/k13,
  diff(x8, t) = (1*k13*k6*x7 + (-1)*k13*(k7*x1*x8 - k12*x9))/k13,
  diff(x9, t) = (1*k13*(k7*x1*x8 - k12*x9) + (-1)*k13*k8*x9)/k13 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file