load_package f4$
torder({}, revgradlex)$

k1 := 10;
k2 := 10;
k3 := 10;
k4 := 10;
k5 := 2;
k6 := 1;
k7 := 3/10;
k8 := 1;
k9 := 1;
k10 := 1;
k11 := 1;
k12 := 1/10;
k13 := 1/10;
k14 := 3/10;
k15 := 1/10;
k16 := 1/20;
k17 := 1/20;
k18 := 1/20;
k19 := 3/10;
k20 := 1;
k21 := 4;
k22 := 1;
k23 := 2;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = 0,
  diff(x3, t) = 0,
  diff(x4, t) = ((-1)*k20*(k1*k21*x4 - k2*x5) + 1*k20*(k5*x6 - k6*x9*x4) + (-1)*k20*k7*x4 + 1*k20*(k12*x11 - k13*x4*x8) + 1*k20*(k17*x7 - k18*x4*x8))/k20,
  diff(x5, t) = (1*k20*(k1*k21*x4 - k2*x5) + (-1)*k20*(k3*k23*x5 - k4*x6) + (-1)*k20*(k15*k22*x5 - k16*x7))/k20,
  diff(x6, t) = (1*k20*(k3*k23*x5 - k4*x6) + (-1)*k20*(k5*x6 - k6*x9*x4))/k20,
  diff(x7, t) = (1*k20*(k15*k22*x5 - k16*x7) + (-1)*k20*(k17*x7 - k18*x4*x8))/k20,
  diff(x8, t) = ((-1)*k20*(k8*x9*x8 - k9*x10) + 1*k20*(k12*x11 - k13*x4*x8) + (-1)*k20*k7*x8 + 1*k20*(k17*x7 - k18*x4*x8))/k20,
  diff(x9, t) = (1*k20*(k5*x6 - k6*x9*x4) + (-1)*k20*(k8*x9*x8 - k9*x10) + (-1)*k20*k7*x9)/k20,
  diff(x10, t) = (1*k20*(k8*x9*x8 - k9*x10) + (-1)*k20*(k10*k22*x10 - k11*x11))/k20,
  diff(x11, t) = (1*k20*(k10*k22*x10 - k11*x11) + (-1)*k20*(k12*x11 - k13*x4*x8))/k20 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file