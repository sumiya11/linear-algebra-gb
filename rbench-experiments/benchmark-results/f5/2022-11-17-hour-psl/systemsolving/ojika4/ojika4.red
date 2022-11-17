% ojika-4 system in revgradlex
% characteristic 0

load_package f5;

system := {
        x1+x3*x1^3+x1*x3*x2^2-x1*x3,
        10*x2-2*x2*x3*x1^2-x3*x2^3-x2*x3,
        -6*x3^2*x1^4-3*x1^2*x2^2*x3^2-x3^2*x1^2+28*x3*x1^2 - 3*x3^2*x2^4+2*x3^2*x2^2+7*x3*x2^2+x3^2-11*x3+10
};

vars := {x1, x2, x3};
torder(vars, revgradlex)$

gb := f5(system)$

end;
