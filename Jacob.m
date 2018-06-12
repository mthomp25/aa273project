%compute jacobian
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 p
xdot = [x4;
    x5;
    x6;
    x1*x10^2 *(1 + 2*x7/p) + 2*x10 *(x5-x2*x8/x7);
    -2*x10 *(x4-x1*x8/x7) + x2*x10^2 *(1-x7/p); 
    -x7*x10^2 *x3/p;
    x8;
    x7*x10^2 *(1-x7/p);
    x10;
    -2*x8*x10/x7];

A = simplify(jacobian(xdot, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10]));