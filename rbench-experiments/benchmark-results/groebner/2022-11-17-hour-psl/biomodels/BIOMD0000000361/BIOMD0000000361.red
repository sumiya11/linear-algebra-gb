load_package groebner$
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
k10 := 33/500;
k11 := 9999997/10000000;
k12 := 1699999/10000;
k13 := 2399999/1000000;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*(k2*x2*x1 - k3*x3) + 1*k1*(k5*x4 - k6*x5*x1) + (-1)*k1*(k9*x1*x7 - k10*x8))/k1,
  diff(x2, t) = (-1)*k1*(k2*x2*x1 - k3*x3)/k1,
  diff(x3, t) = (1*k1*(k2*x2*x1 - k3*x3) + (-1)*k1*k4*x3)/k1,
  diff(x4, t) = (1*k1*k4*x3 + (-1)*k1*(k5*x4 - k6*x5*x1))/k1,
  diff(x5, t) = (1*k1*(k5*x4 - k6*x5*x1) + (-1)*k1*(k7*x5*x6 - k8*x7))/k1,
  diff(x6, t) = (-1)*k1*(k7*x5*x6 - k8*x7)/k1,
  diff(x7, t) = (1*k1*(k7*x5*x6 - k8*x7) + (-1)*k1*(k9*x1*x7 - k10*x8))/k1,
  diff(x8, t) = 1*k1*(k9*x1*x7 - k10*x8)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file