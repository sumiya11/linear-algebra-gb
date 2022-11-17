% cyclic-5 system in revgradlex
% characteristic 0
% 0 dim
%

load_package f5;

system := {
        z1 + z2 + z3 + z4 + z5,
        z1*z2 + z1*z5 + z2*z3 + z3*z4 + z4*z5,
        z1*z2*z3 + z1*z2*z5 + z1*z4*z5 + z2*z3*z4 + z3*z4*z5,
        z1*z2*z3*z4 + z1*z2*z3*z5 + z1*z2*z4*z5 + z1*z3*z4*z5 + z2*z3*z4*z5,
        z1*z2*z3*z4*z5 - 1
}$

vars := {z1,z2,z3,z4,z5}$
torder(vars, revgradlex)$

gb := f5(system)$

end;
