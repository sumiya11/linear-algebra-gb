on echo; 

load_package f4;
load_package f5;
load_package groebner;
load_package profile;



system := {
y**2*z+2*x*y*t-2*x-z,
-x**3*z+4*x*y**2*z+4*x**2*y*t+2*y**3*t+4*x**2-10*y**2+4*x*z-10*y*t+2,
2*y*z*t+x*t**2-x-2*z,
-x*z**3+4*y*z**2*t+4*x*z*t**2+2*y*t**3+4*x*z+4*z**2-10*y*t-10*t**2+2
}$

vars := {x, y, z, t}$
torder(vars, revgradlex)$


on backtrace;

profile f4_reduction, f4_select_normal, f4_symbolic_preprocessing, 
            basis_update;

% profile f4_reduction;
% profile matrix_convert_matrix_rows_to_basis_elements, matrix_linear_algebra,
%         sorting_sort_matrix_rows_increasing, sorting_sort_matrix_rows_decreasing,
%         matrix_convert_hashes_to_columns;

on time;
f4(system);

proprint();



end;