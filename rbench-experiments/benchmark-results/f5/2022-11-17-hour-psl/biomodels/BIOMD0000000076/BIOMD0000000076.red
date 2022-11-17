load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 47;
k3 := 23/1000;
k4 := 27/50;
k5 := 187/100;
k6 := 29/20;
k7 := 10000;
k8 := 601/100;
k9 := 24/5;
k10 := 237/100;
k11 := 73/100;
k12 := 217/100;
k13 := 2;
k14 := 93/100;
k15 := 6/5;
k16 := 53;
k17 := 7/2;
k18 := 1;
k19 := 1;
k20 := 151/10;
k21 := 59/100;
operator diff$
odes := { diff(x1, t) = (1*k1*k2/(k3*k4)*(k5*k21 - k6*x1/k7)/((1 + k8/k9 + k10/k11 + k12/k13)*(1 + k5/k3 + k6/k14)*(1 + k21/k4 + x1/k15)) + (-1)*k1*k16*x1/k17/((1 + x1/k17)*(1 + k18/k19)))/k1,
  diff(x2, t) = 0,
  diff(x3, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file