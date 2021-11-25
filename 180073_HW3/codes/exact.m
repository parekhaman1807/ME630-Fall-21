% AMAN PAREKH - 180073 - ME630 - Monsoon 2021
% Disclaimer: This code can take some time to run due to symbolic integration

clear;
clc;
syms eta expr1 expr2 expr3 expr4
x = -1:0.01:1;
re = 50;
t = 0.36;
for i = 1:201
    expr1 = exp(-0.5*re*(eta + ((0.5*(x(i)-eta)^2)/t)));
    expr2 = exp(-0.5*re*((0.5*(x(i)-eta)^2)/t));
    expr3 = ((x(i)-eta)/t)*exp(-0.5*re*(eta + ((0.5*(x(i)-eta)^2)/t)));
    expr4 = ((x(i)-eta)/t)*exp(-0.5*re*((0.5*(x(i)-eta)^2)/t));

    n = int(expr3, eta, -inf, 0) + int(expr4, eta, 0, inf);
    d = int(expr1, eta, -inf, 0) + int(expr2, eta, 0, inf);

    val(i) = n/d;
end
val = double(val);
