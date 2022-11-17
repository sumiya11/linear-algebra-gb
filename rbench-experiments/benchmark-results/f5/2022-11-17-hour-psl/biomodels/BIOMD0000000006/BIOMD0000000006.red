load_package f5$
torder({}, revgradlex)$

k1 := 3/200;
k2 := 1;
k3 := 180;
k4 := 9/500;
k5 := 1/10000;
k6 := 1;
k7 := 1;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = k3*(x4 - x2)*(k4/k3 + x2^2) - k2*x2,
  diff(x4, t) = k1 - k2*x2 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file