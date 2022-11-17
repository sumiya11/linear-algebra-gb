load_package f4$
torder({}, revgradlex)$

k1 := 9/50;
k2 := 23/5;
k3 := 101/1000;
k4 := 1/125;
k5 := 3/2;
k6 := 1/5;
k7 := 3/10;
k8 := 7/5;
k9 := 103/2500;
k10 := 3/10;
k11 := 1/20;
k12 := 9/20;
k13 := 3/100;
k14 := 2/5;
k15 := 3/10;
k16 := 7/20;
k17 := 3/10;
k18 := 1/2;
k19 := 1;
operator diff$
odes := { diff(x1, t) = (1*k19*k1*x1 + (-1)*k19*(k2*x1*x1 + k3*x1*x3 + k4*x1*x5))/k19,
  diff(x2, t) = (1*k19*(k5*x1 + k7*x1*x2) + (-1)*k19*k6*x2)/k19,
  diff(x3, t) = (1*k19*(k8*x1 + k10*x1*x3 + k11*x2*x3) + (-1)*k19*k9*x3)/k19,
  diff(x4, t) = (1*k19*(k12*x1 + k14*x1*x4 + k15*x2*x4) + (-1)*k19*k13*x4)/k19,
  diff(x5, t) = (1*k19*k16*x4 + (-1)*k19*(k17*x5 + k18*x1*x5))/k19 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file