load_package f4$
torder({}, revgradlex)$

k1 := 1;
k2 := 3/10;
k3 := 5;
k4 := 1;
k5 := 1;
k6 := 1;
k7 := 1;
operator diff$
odes := { diff(x1, t) = (1*k4*k1/(k6^k3 + x2^k3) + (-1)*k4*x1/k2/(1 + x1/k2))/k4,
  diff(x2, t) = (1*k4*k1/(k7^k3 + x3^k3) + (-1)*k4*x2/k2/(1 + x2/k2))/k4,
  diff(x3, t) = (1*k4*k1/(k5^k3 + x1^k3) + (-1)*k4*x3/k2/(1 + x3/k2))/k4 }$
odes := for each o in odes collect part(o, 2)$

gb := f4(odes)$

end; % of file