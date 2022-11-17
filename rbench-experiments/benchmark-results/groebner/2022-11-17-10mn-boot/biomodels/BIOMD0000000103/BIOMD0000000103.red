load_package groebner$
torder({}, revgradlex)$

k1 := 1/500;
k2 := 1/10;
k3 := 1/200000;
k4 := 7/20000;
k5 := 1/1000;
k6 := 1/1000;
k7 := 1/1000;
k8 := 1/1000;
k9 := 1/500;
k10 := 1/10;
k11 := 3/1000;
k12 := 1/1000;
k13 := 1/5000;
k14 := 1/5000;
k15 := 1/500;
k16 := 1/10;
k17 := 1/1000;
k18 := 1/1000;
k19 := 1/1000;
k20 := 1/1000;
k21 := 1/500;
k22 := 1/10;
k23 := 1/20000;
k24 := 7/2000;
k25 := 1/1000;
k26 := 1/50;
k27 := 1/1000;
k28 := 1/50;
k29 := 1/1000;
k30 := 1/25;
k31 := 1/1000;
k32 := 1/1000;
k33 := 1/1000;
k34 := 1/1000;
k35 := 1/5;
k36 := 1/1000;
k37 := 1/1000;
k38 := 1/1000;
k39 := 1/1000;
k40 := 1/1000;
k41 := 1/1000;
k42 := 1/1000;
k43 := 1/1000;
k44 := 1/1000;
k45 := 1/1000;
k46 := 1;
k47 := 1;
k48 := 1;
k49 := 1;
k50 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k50*(k1*x1*x2 - k2*x6) + (-1)*k50*(k9*x3*x1 - k10*x5) + (-1)*k50*(k15*x11*x1 - k16*x12) + (-1)*k50*(k21*x10*x1 - k22*x13) + 1*k50*(k26 - k25*x1) + (-1)*k50*(k46*k1*x14*x1 - k47*k2*x15) + (-1)*k50*(k46*k1*x16*x1 - k47*k2*x17))/k50,
  diff(x2, t) = ((-1)*k50*(k1*x1*x2 - k2*x6) + (-1)*k50*(k5*x2*x4 - k6*x3) + (-1)*k50*k13*x2*x8 + 1*k50*(k28 - k27*x2) + (-1)*k50*(k46*k5*x2*x9 - k47*k6*x14))/k50,
  diff(x3, t) = (1*k50*(k5*x2*x4 - k6*x3) + (-1)*k50*(k9*x3*x1 - k10*x5) + (-1)*k50*k31*x3 + (-1)*k50*(k46*k11*x8*x3 - k47*k12*x14))/k50,
  diff(x4, t) = ((-1)*k50*(k5*x2*x4 - k6*x3) + (-1)*k50*(k7*x6*x4 - k8*x5) + (-1)*k50*(k11*x8*x4 - k12*x9) + (-1)*k50*(k17*x11*x4 - k18*x10) + (-1)*k50*(k19*x12*x4 - k20*x13) + 1*k50*(k30 - k29*x4))/k50,
  diff(x5, t) = (1*k50*(k7*x6*x4 - k8*x5) + 1*k50*(k9*x3*x1 - k10*x5) + (-1)*k50*k32*x5 + (-1)*k50*(k46*k11*x8*x5 - k47*k12*x15))/k50,
  diff(x6, t) = (1*k50*(k1*x1*x2 - k2*x6) + (-1)*k50*(k7*x6*x4 - k8*x5) + (-1)*k50*k14*x6*x8 + (-1)*k50*k33*x6 + (-1)*k50*(k46*k5*x6*x9 - k47*k6*x15))/k50,
  diff(x7, t) = ((-1)*k50*k3*x7*x2 + (-1)*k50*k4*x7*x6 + (-1)*k50*k23*x7*x11 + (-1)*k50*k24*x7*x12 + 1*k50*(k35 - k34*x7))/k50,
  diff(x8, t) = (1*k50*k3*x7*x2 + 1*k50*k4*x7*x6 + (-1)*k50*(k11*x8*x4 - k12*x9) + 1*k50*k23*x7*x11 + 1*k50*k24*x7*x12 + (-1)*k50*k36*x8 + (-1)*k50*(k46*k11*x8*x3 - k47*k12*x14) + (-1)*k50*(k46*k11*x8*x5 - k47*k12*x15) + (-1)*k50*(k46*k11*x8*x10 - k47*k12*x16) + (-1)*k50*(k46*k11*x8*x13 - k47*k12*x17))/k50,
  diff(x9, t) = (1*k50*(k11*x8*x4 - k12*x9) + (-1)*k50*k37*x9 + (-1)*k50*(k46*k5*x2*x9 - k47*k6*x14) + (-1)*k50*(k46*k5*x6*x9 - k47*k6*x15) + (-1)*k50*(k46*k5*x11*x9 - k47*k6*x16) + (-1)*k50*(k46*k5*x12*x9 - k47*k6*x17))/k50,
  diff(x10, t) = (1*k50*(k17*x11*x4 - k18*x10) + (-1)*k50*(k21*x10*x1 - k22*x13) + (-1)*k50*k38*x10 + (-1)*k50*(k46*k11*x8*x10 - k47*k12*x16))/k50,
  diff(x11, t) = (1*k50*k13*x2*x8 + (-1)*k50*(k15*x11*x1 - k16*x12) + (-1)*k50*(k17*x11*x4 - k18*x10) + (-1)*k50*k39*x11 + (-1)*k50*(k46*k5*x11*x9 - k47*k6*x16))/k50,
  diff(x12, t) = (1*k50*k14*x6*x8 + 1*k50*(k15*x11*x1 - k16*x12) + (-1)*k50*(k19*x12*x4 - k20*x13) + (-1)*k50*k40*x12 + (-1)*k50*(k46*k5*x12*x9 - k47*k6*x17))/k50,
  diff(x13, t) = (1*k50*(k19*x12*x4 - k20*x13) + 1*k50*(k21*x10*x1 - k22*x13) + (-1)*k50*k41*x13 + (-1)*k50*(k46*k11*x8*x13 - k47*k12*x17))/k50,
  diff(x14, t) = ((-1)*k50*k42*x14 + 1*k50*(k46*k11*x8*x3 - k47*k12*x14) + 1*k50*(k46*k5*x2*x9 - k47*k6*x14) + (-1)*k50*(k46*k1*x14*x1 - k47*k2*x15))/k50,
  diff(x15, t) = (1*k50*(k46*k11*x8*x5 - k47*k12*x15) + 1*k50*(k46*k5*x6*x9 - k47*k6*x15) + 1*k50*(k46*k1*x14*x1 - k47*k2*x15))/k50,
  diff(x16, t) = ((-1)*k50*k44*x16 + 1*k50*(k46*k11*x8*x10 - k47*k12*x16) + 1*k50*(k46*k5*x11*x9 - k47*k6*x16) + (-1)*k50*(k46*k1*x16*x1 - k47*k2*x17))/k50,
  diff(x17, t) = ((-1)*k50*k43*x17 + (-1)*k50*k45*x17 + 1*k50*(k46*k11*x8*x13 - k47*k12*x17) + 1*k50*(k46*k5*x12*x9 - k47*k6*x17) + 1*k50*(k46*k1*x16*x1 - k47*k2*x17))/k50 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file