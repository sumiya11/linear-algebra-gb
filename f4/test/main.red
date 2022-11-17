
load_package f4;

procedure errorMessage(system);
  {"Wrong Answer", system};

on f4modular$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4({x + y}, {x, y}, revgradlex);

f4({x, x}, {x}, revgradlex);

f4({x + y, x}, {x, y}, revgradlex);

f4({x*y + 1, y*z + 1}, {x, y, z}, revgradlex);

f4({x*y + 1, x + y}, {x, y}, revgradlex);

f4({x^2*y + 1, x*y^2 + 1}, {x, y}, revgradlex);

f4({x^5 - x, x^4 - x}, {x, y}, revgradlex);

f4({x^5 - x, x^4 - x, x^3 - x}, {x}, revgradlex);

f4({x*y - x, x*y^2 + y}, {x, y}, revgradlex);

f4({x*y - x}, {a, b, c, d, e, f, g, x, y}, revgradlex);
f4({x*y - x, x*y^2 + y}, {a, b, c, d, e, f, g, x, y}, revgradlex);

f4({x*y*z - 1, x*y + x*z + y*z}, {x, y, z}, revgradlex);

f4({x*y*z + 1, x*y + x*z + y*z, x + y + z}, {x, y, z}, revgradlex);

f4({10*x1*x2^2 - 11*x1 + 10, 10*x1^2*x2 - 11*x2 + 10}, {x1,x2}, revgradlex);

f4({x0^2 - x0 + 2*x1^2 + 2*x2^2, 2*x0*x1 + 2*x1*x2 - x1, x0 + 2*x1 + 2*x2 - 1}, {x0, x1, x2}, revgradlex);

f4({x*y + 1, x, y}, {x, y, z}, revgradlex);

f4({x*y*z + 1, x, y, z}, {x, y, z}, revgradlex);

f4({x*y*z + 1, x + 1}, {x, y, z}, revgradlex);

f4({x*y*z + 1, x*y + 1}, {x, y, z}, revgradlex);

f4({x*y*z + 1, x*y + x*z}, {x, y, z}, revgradlex);

f4({x*y*z + 1, x*y + x*z + y*z}, {x, y, z}, revgradlex);

f4({x*y*z + 1, x*y + x*z + y*z, x + y + z}, {x, y, z}, revgradlex);

f4({x, x^2, x^3, x*y, x*z, x*y^3, z^4, z^12, x*y*y}, {x,y,z}, revgradlex);

f4({10*x1*x2^2 - 11*x1 + 10, 10*x1^2*x2 - 11*x2 + 10}, {x0,x1,x2}, revgradlex);

f4({x0^2 - x0 + 2*x1^2 + 2*x2^2, 2*x0*x1 + 2*x1*x2 - x1, x0 + 2*x1 + 2*x2 - 1}, {x0,x1,x2}, revgradlex);

f4({x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2,
 2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3,
 2*x0*x2 + x1^2 + 2*x1*x3 - x2,
 x0 + 2*x1 + 2*x2 + 2*x3 - 1}, 
    {x0,x1,x2,x3}, revgradlex);

f4({x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 - x0, 
    2*x0*x1 + 2*x1*x2 + 2*x2*x3 - x1, 
    x1^2 + 2*x0*x2 + 2*x1*x3 - x2, 
    x0 + 2*x1 + 2*x2 + 2*x3 - 1},
    {x0,x1,x2,x3}, revgradlex);

kat := {x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 - x0, 2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 - x1, x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 + 2*x3*x5 - x2, 2*x1*x2 + 2*x0*x3 + 2*x1*x4 + 2*x2*x5 - x3, x2^2 + 2*x1*x3 + 2*x0*x4 + 2*x1*x5 - x4, x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 - 1};
f4(kat, {x0,x1,x2,x3,x4,x5}, revgradlex);

noon := {10*x1*x2^2 + 10*x1*x3^2 + 10*x1*x4^2 + 10*x1*x5^2 - 11*x1 + 10, 10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 + 10*x2*x5^2 - 11*x2 + 10, 10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 + 10*x3*x5^2 - 11*x3 + 10, 10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 + 10*x4*x5^2 - 11*x4 + 10, 10*x1^2*x5 + 10*x2^2*x5 + 10*x3^2*x5 + 10*x4^2*x5 - 11*x5 + 10};
f4(noon, {x1,x2,x3,x4,x5}, revgradlex);

f4({x + 100000000000000000000000000000000000000000000000000000000000000000000000000000000});
f4({x + 1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000});
f4({x + 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000});

% katsura 6
s := {
   x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2,
 2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6,
 2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5 + 2*x4*x6,
 2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 + 2*x3*x6 - x3,
 2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 + 2*x2*x6 - x4,
 2*x0*x5 + 2*x1*x4 + 2*x1*x6 + 2*x2*x3 - x5,
 x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 - 1
};
vars := {x0, x1, x2, x3, x4, x5, x6};
torder(vars, revgradlex)$
gb := f4(s)$

if not (length(gb) = 41) then
  errorMessage("Katsura-6 error");

% katsura 7
s := {
  x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2,
  2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7,
  2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5 + 2*x4*x6 + 2*x5*x7,
  2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 + 2*x3*x6 - x3 + 2*x4*x7,
  2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 + 2*x2*x6 + 2*x3*x7 - x4,
  2*x0*x5 + 2*x1*x4 + 2*x1*x6 + 2*x2*x3 + 2*x2*x7 - x5,
  2*x0*x6 + 2*x1*x5 + 2*x1*x7 + 2*x2*x4 + x3^2 - x6,
  x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 - 1}$
% f4(s, {x0, x1, x2, x3, x4, x5, x6, x7}, revgradlex)$

% noon 6
s := {
  10*x1*x2^2 + 10*x1*x3^2 + 10*x1*x4^2 + 10*x1*x5^2 + 10*x1*x6^2 - 11*x1 + 10,
  10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 + 10*x2*x5^2 + 10*x2*x6^2 - 11*x2 + 10,
  10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 + 10*x3*x5^2 + 10*x3*x6^2 - 11*x3 + 10,
  10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 + 10*x4*x5^2 + 10*x4*x6^2 - 11*x4 + 10,
  10*x1^2*x5 + 10*x2^2*x5 + 10*x3^2*x5 + 10*x4^2*x5 + 10*x5*x6^2 - 11*x5 + 10,
  10*x1^2*x6 + 10*x2^2*x6 + 10*x3^2*x6 + 10*x4^2*x6 + 10*x5^2*x6 - 11*x6 + 10};
f4(s, {x1, x2, x3, x4, x5, x6}, revgradlex);

% eco 7
s := {
   x1*x2*x7 + x1*x7 + x2*x3*x7 + x3*x4*x7 + x4*x5*x7 + x5*x6*x7 - 1,
 x1*x3*x7 + x2*x4*x7 + x2*x7 + x3*x5*x7 + x4*x6*x7 - 2,
 x1*x4*x7 + x2*x5*x7 + x3*x6*x7 + x3*x7 - 3,
 x1*x5*x7 + x2*x6*x7 + x4*x7 - 4,
 x1*x6*x7 + x5*x7 - 5,
 x6*x7 - 6,
 x1 + x2 + x3 + x4 + x5 + x6 + 1
};
vars := {x1, x2, x3, x4, x5, x6, x7};
torder(vars, revgradlex)$
gb := f4(s)$

if not (length(gb) = 32) then
  errorMessage("Eco-7 error");

% eco-10
s := {
  x1*x2*x10 + x1*x10 + x2*x3*x10 + x3*x4*x10 + x4*x5*x10 + x5*x6*x10 + x6*x7*x10 + x7*x8*x10 + x8*x9*x10 - 1,
 x1*x3*x10 + x2*x4*x10 + x2*x10 + x3*x5*x10 + x4*x6*x10 + x5*x7*x10 + x6*x8*x10 + x7*x9*x10 - 2,
 x1*x4*x10 + x2*x5*x10 + x3*x6*x10 + x3*x10 + x4*x7*x10 + x5*x8*x10 + x6*x9*x10 - 3,
 x1*x5*x10 + x2*x6*x10 + x3*x7*x10 + x4*x8*x10 + x4*x10 + x5*x9*x10 - 4,
 x1*x6*x10 + x2*x7*x10 + x3*x8*x10 + x4*x9*x10 + x5*x10 - 5,
 x1*x7*x10 + x2*x8*x10 + x3*x9*x10 + x6*x10 - 6,
 x1*x8*x10 + x2*x9*x10 + x7*x10 - 7,
 x1*x9*x10 + x8*x10 - 8,
 x9*x10 - 9,
 x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 1
};
vars := {x1, x2, x3, x4, x5, x6, x7, x8, x9, x10};
torder(vars, revgradlex)$
% gb := f4(s)$

% if not (length(gb) = 203) then
%   errorMessage("Eco-10 error");

% noon 7
s := {
   10*x1*x2^2 + 10*x1*x3^2 + 10*x1*x4^2 + 10*x1*x5^2 + 10*x1*x6^2 + 10*x1*x7^2 - 11*x1 + 10,
 10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 + 10*x2*x5^2 + 10*x2*x6^2 + 10*x2*x7^2 - 11*x2 + 10,
 10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 + 10*x3*x5^2 + 10*x3*x6^2 + 10*x3*x7^2 - 11*x3 + 10,
 10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 + 10*x4*x5^2 + 10*x4*x6^2 + 10*x4*x7^2 - 11*x4 + 10,
 10*x1^2*x5 + 10*x2^2*x5 + 10*x3^2*x5 + 10*x4^2*x5 + 10*x5*x6^2 + 10*x5*x7^2 - 11*x5 + 10,
 10*x1^2*x6 + 10*x2^2*x6 + 10*x3^2*x6 + 10*x4^2*x6 + 10*x5^2*x6 + 10*x6*x7^2 - 11*x6 + 10,
 10*x1^2*x7 + 10*x2^2*x7 + 10*x3^2*x7 + 10*x4^2*x7 + 10*x5^2*x7 + 10*x6^2*x7 - 11*x7 + 10
};
vars := {x1, x2, x3, x4, x5, x6, x7};
torder(vars, revgradlex)$
% gb := f4(s)$

% if not (length(gb) = 495) then
%   errorMessage("Noon-7 error");

end; % eof