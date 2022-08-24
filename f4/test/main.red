
load_package f4;

procedure errorMessage(system);
  {"Wrong Answer", system};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bugs fixed: 

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

f4({x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 - x0, 
    2*x0*x1 + 2*x1*x2 + 2*x2*x3 - x1, 
    x1^2 + 2*x0*x2 + 2*x1*x3 - x2, 
    x0 + 2*x1 + 2*x2 + 2*x3 - 1},
    {x0,x1,x2,x3}, revgradlex);

kat := {x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 - x0, 2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 - x1, x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 + 2*x3*x5 - x2, 2*x1*x2 + 2*x0*x3 + 2*x1*x4 + 2*x2*x5 - x3, x2^2 + 2*x1*x3 + 2*x0*x4 + 2*x1*x5 - x4, x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 - 1};
f4(kat, {x0,x1,x2,x3,x4,x5}, revgradlex);

noon := {10*x1*x2^2 + 10*x1*x3^2 + 10*x1*x4^2 + 10*x1*x5^2 - 11*x1 + 10, 10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 + 10*x2*x5^2 - 11*x2 + 10, 10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 + 10*x3*x5^2 - 11*x3 + 10, 10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 + 10*x4*x5^2 - 11*x4 + 10, 10*x1^2*x5 + 10*x2^2*x5 + 10*x3^2*x5 + 10*x4^2*x5 - 11*x5 + 10};
f4(noon, {x1,x2,x3,x4,x5}, revgradlex);

end; % eof
