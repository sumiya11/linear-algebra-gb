load_package f5$
torder({}, revgradlex)$

k1 := 1/10;
k2 := 9/10;
k3 := 1/20;
k4 := 9/1000;
k5 := 27/40;
k6 := 1/100;
k7 := 1/200;
k8 := 1/10;
k9 := 3/10;
k10 := 1/100;
k11 := 1/10;
k12 := 1/2;
k13 := 1/40;
k14 := 1;
operator diff$
odes := { diff(x1, t) = k3/(k11 + x2) - k8*x1,
  diff(x2, t) = k1*((k4 + k5*x1)*x2*x3/(k12 + x1^2) - k9*x2),
  diff(x3, t) = k1*k2*(k6*x1 - (k10*x3 + k7*x1*x3/(k13 + x1))) }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file