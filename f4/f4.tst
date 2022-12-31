% Simple correctness tests of f4

% Usual computation in Q[x, y]
f4({x^2*y - 10x + 5, 3y^2 + 4x + 7});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different term orders.

% For the detailed description of all supported term orders we refer 
% to the groebner package documentation, subsection "Term Ordering Mode".
% Here, we demonstrate basic use-cases.

% Without torder:

f4({x + z + y});
% Expected {x + y + z}

f4({x + z + y^2}, {z, y, x}, revgradlex);
% Expected {y^2 + z + x}

% With torder:

torder({z, x, y}, lex);
f4({x + y^2 + z});
% Expected {z + x + y^2}

% Weighted order
torder({x,y,z},weighted,{3, 5, 4});
f4({x + y + z});
% Expected {y + z + x}

% Block order
torder({x, y, z}, lexgradlex, 1);
f4({x + y + z^2});
% Expected {x + z^2 + y}

% Graded order
torder({x, y, z}, graded, {1, 1, 2}, lex);
f4({x + y^4 + z^3});
% Expected {z^3 + y^4 + x}

% and Matrix order
torder({x, y, z},matrix, mat(
    (0,3,1),
    (2,0,0),
    (0,0,1)));
f4({x + y + z});
% Expected {y + z + x}
f4({x^2 + y + z^6});
% Expected {z^6 + y + x^2}

% ..or compile the matrix for greater efficiency
torder_compile(m, mat(
  (1,1,1),
  (0,0,1),
  (0,1,0)
));
torder({x,y,z}, m);
f4({x + y^2 + z^2});
% Expected {z^2 + y^2 + x}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test switches.
% Currently, there is only 1 documented switch: f4modular.
% Note that results must be identical for both states of f4modular.

off f4modular;
f4(noon3, {x1, x2, x3}, revgradlex);
f4({5x1 + x2, x1*x2 + 1}, {x1, x2}, lex);
f4({x1 + x2, x1*x2 + 100}, {x1, x2}, lex);
f4({4x1 + x2, 1234x1*x2 + 1e5}, {x1, x2}, lex);
f4({x1 + x2, x1*x2 + 1e10}, {x1, x2}, lex);
f4({x1 + x2, 2347624x1*x2 + 1e100}, {x1, x2}, lex);

on f4modular;
f4(noon3, {x1, x2, x3}, revgradlex);
f4({5x1 + x2, x1*x2 + 1}, {x1, x2}, lex);
f4({x1 + x2, x1*x2 + 100}, {x1, x2}, lex);
f4({4x1 + x2, 1234x1*x2 + 1e5}, {x1, x2}, lex);
% Test that modular works with large numbers
f4({x1 + x2, x1*x2 + 1e10}, {x1, x2}, lex);
f4({x1 + x2, 2347624x1*x2 + 1e100}, {x1, x2}, lex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sanity checks

torder({x1}, lex);

f4({3x1});

f4({4x1}, {x1}, lex);

torder({x1, x2}, revgradlex);

f4({x1, x2^2});

f4({x1 + 1, x1});

f4({x2, x1, x1, x1, x2, x2, x1, x2, x1, x2}, {x1, x2}, lex);

f4({x1 + x2, (x1 + x2)^2, (x1 + x2)^3}, {x1, x2}, lex);

f4({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

f4({x1*x2 + 1, x2*x3 + 1}, {x1, x2, x3}, lex);

f4({x1 + x2 + x3, x1*x2 + x2*x3 + x1*x3, x1*x2*x3 - 1}, {x1, x2, x3}, lex);

f4({10*x1*x2^2 - 11*x1 + 10, 10*x1^2*x2 - 11*x2 + 10}, {x1, x2}, lex);

noon3 := {10*x1*x2^2 + 10*x1*x3^2 - 11*x1 + 10,
          10*x1^2*x2 + 10*x2*x3^2 - 11*x2 + 10,
          10*x1^2*x3 + 10*x2^2*x3 - 11*x3 + 10}$
f4(noon3, {x1, x2, x3}, revgradlex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test from groebner Reduce package

vars := {q1,q2,q3,q4,q5,q6}$
system := {q1,
          q2**2 + q3**2 + q4**2,
          q4*q3*q2,
          q3**2*q2**2 + q4**2*q2**2 + q4**2*q3**2,
          q6**2 + 1/3*q5**2,
          q6**3 - q5**2*q6,
          2*q2**2*q6 - q3**2*q6 - q4**2*q6 + q3**2*q5 - q4**2*q5,
          2*q2**2*q6**2 - q3**2*q6**2 - q4**2*q6**2 - 2*q3**2*q5*q6
          + 2*q4**2*q5*q6 - 2/3*q2**2*q5**2 + 1/3*q3**2*q5**2
          + 1/3*q4**2*q5**2,
          - q3**2*q2**2*q6 - q4**2*q2**2*q6 + 2*q4**2*q3**2*q6 -
          q3**2*q2**2*q5 + q4**2*q2**2*q5,
          - q3**2*q2**2*q6**2 - q4**2*q2**2*q6**2 + 2*q4**2*q3**2*q6**2
          + 2*q3**2*q2**2*q5*q6 - 2*q4**2*q2**2*q5*q6 + 1/3*q3**2*q2**2
          *q5**2 + 1/3*q4**2*q2**2*q5**2 - 2/3*q4**2*q3**2*q5**2,
          - 3*q3**2*q2**4*q5*q6**2 + 3*q4**2*q2**4*q5*q6**2
          + 3*q3**4*q2**2*q5*q6**2 - 3*q4**4*q2**2*q5*q6**2
          - 3*q4**2*q3**4*q5*q6**2 + 3*q4**4*q3**2*q5*q6**2
          + 1/3*q3**2*q2**4*q5**3 - 1/3*q4**2*q2**4*q5**3
          - 1/3*q3**4*q2**2*q5**3 + 1/3*q4**4*q2**2*q5**3 + 1/3*q4**2
            *q3**4*q5**3 - 1/3*q4**4*q3**2*q5**3}$

f4(system, vars, lex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests from Groebner.jl Julia package
%   https://github.com/sumiya11/Groebner.jl

% Root-5
system := {
  x1 + x2 + x3 + x4 + x5,
  x1*x2 + x1*x3 + x1*x4 + x1*x5 + x2*x3 + x2*x4 + x2*x5 + x3*x4 + x3*x5 + x4*x5,
  x1*x2*x3 + x1*x2*x4 + x1*x2*x5 + x1*x3*x4 + x1*x3*x5 + x1*x4*x5 + x2*x3*x4 + x2*x3*x5 + x2*x4*x5 + x3*x4*x5,
  x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5,
  x1*x2*x3*x4*x5 - 1
}$
vars := {x1, x2, x3, x4, x5}$
f4(system, vars, lex);

% Sparse-5
system := {
  x1^2*x2^2*x3^2*x4^2*x5^2 + 3*x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + x4^2 + x5^2 + 5,
  x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + 3*x2^2 + x3^2 + x4^2 + x5^2 + 5,
  x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + 3*x3^2 + x4^2 + x5^2 + 5,
  x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + 3*x4^2 + x5^2 + 5,
  x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + x4^2 + 3*x5^2 + 5
}$
vars := {x1, x2, x3, x4, x5}$
f4(system, vars, revgradlex);

% Katsura-5
system := {
  x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2,
  2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5,
  2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5,
  2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 - x3,
  2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 - x4,
  x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 - 1
}$
vars := {x0,x1,x2,x3,x4,x5}$
f4(system, vars, revgradlex);

% Katsura-6
system := {
 x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2,
 2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6,
 2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5 + 2*x4*x6,
 2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 + 2*x3*x6 - x3,
 2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 + 2*x2*x6 - x4,
 2*x0*x5 + 2*x1*x4 + 2*x1*x6 + 2*x2*x3 - x5,
 x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 - 1
}$
vars := {x0,x1,x2,x3,x4,x5,x6}$
f4(system, vars, revgradlex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test symbolic mode entry point

expr1 := (x + y^2);
expr2 := 2*x - 1;
share expr1, expr2;

lisp;

torder({{'list, 'y, 'z, 'x}, 'lex});

f4_groebnerq({simp 'x, simp 'y});
f4_groebnerq({simp expr1, simp expr2});

f4_groebnerf({numr simp 'x, numr simp 'z});
f4_groebnerf({numr simp expr1, numr simp expr2});
f4_groebnerf({numr simp expr1, denr simp expr2}); % this GB is just (1)

algebraic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ad-hoc, test that consecutive torder's change the order correctly 

torder({x, y, z},matrix, mat(
    (0,3,1),
    (2,0,0),
    (0,0,1)));
f4({x + y + z});
% Expected y + z + x

torder({x,y,z},weighted,{3, 5, 4});
f4({x + y + z});
% Expected {y + z + x}

torder_compile(m, mat(
  (1,1,1),
  (0,0,1),
  (0,1,0)
));
torder({x,y,z}, m);
f4({x + y^2 + z^2});
% Expected {z^2 + y^2 + x}

% Graded order
torder({x, y, z}, graded, {1, 1, 2}, lex);
f4({x + y^4 + z^3});
% Expected {z^3 + y^4 + x}

% Block order
torder({x, y, z}, lexgradlex, 1);
f4({x + y + z^2});
% Expected {x + z^2 + y}

f4({x + z + y^2}, {z, y, x}, revgradlex);
% Expected {y^2 + z + x}

end;  % of file
