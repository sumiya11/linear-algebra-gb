load_package f4$
torder({}, revgradlex)$

k1 := -247/100;
k2 := 1097/50;
k3 := 641/100;
k4 := 7/4;
k5 := 4/5;
k6 := 1/20;
k7 := 693/1000;
k8 := 1/500;
k9 := 7/100;
k10 := 1/5;
k11 := 91/10;
k12 := 1;
operator diff$
odes := { diff(x2, t) = k5*x1 + k6*k7*x2*(1 - k8*(x2 + x3)) - k10*x2,
  diff(x3, t) = (1 - k5)*x1 + k6*k7*x3*(1 - k8*(x2 + x3)) - k10*x3,
  diff(x4, t) = k9*x2 - k11*x4 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file