load_package f4$
torder({}, revgradlex)$

k1 := 60;
k2 := 3333333333/100000000000;
k3 := 14;
k4 := 9/800000000000;
k5 := 9/200000000000;
k6 := 1;
operator diff$
odes := { diff(x1, t) = (1*k6*1/4*k5*(x2 + x3)*x1 + (-1)*k6*1/k3*x1)/k6,
  diff(x2, t) = (1*k6*k2*x2 + (-1)*k6*k5*x2*x1)/k6,
  diff(x3, t) = (-1)*k6*(k5*x3*x1 + 1/k1*x3)/k6 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file