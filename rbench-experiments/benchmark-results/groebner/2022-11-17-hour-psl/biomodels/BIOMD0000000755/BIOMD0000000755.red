load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 121267/1000;
k3 := 94929/200000000000000000000;
k4 := 2569840000000;
k5 := 387701/5;
k6 := 69679400000;
k7 := 472749/100000000;
k8 := 201671/10000000;
k9 := 1/1000000000;
k10 := 1/6250000;
k11 := 7/5000000;
operator diff$
odes := { diff(x1, t) = (-1)*k1*k2*x1*x2/k1,
  diff(x2, t) = (-1)*k1*k2*x1*x2/k1,
  diff(x3, t) = (1*k1*k2*x1*x2 + (-1)*k1*k4*x3*x4 + 1*k1*k5*x6)/k1,
  diff(x4, t) = ((-1)*k1*k3*x1*x2*x4 + (-1)*k1*k4*x3*x4)/k1,
  diff(x5, t) = (1*k1*k3*x1*x2*x4 + 1*k1*k6*x3*x7 + (-1)*k1*k8*x5)/k1,
  diff(x6, t) = (1*k1*k4*x3*x4 + (-1)*k1*k5*x6)/k1,
  diff(x7, t) = (1*k1*k5*x6 + (-1)*k1*k6*x3*x7 + (-1)*k1*k7*x7)/k1,
  diff(x8, t) = 1*k1*k7*x7/k1,
  diff(x9, t) = 1*k1*k8*x5/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file