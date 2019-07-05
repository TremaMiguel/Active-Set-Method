function [x, lambda] = pcnulo(Q,A,c,b)
% Método del espacio nulo para el problema
% Min: (1/2)*x'*Q*x + c'*x
% s.a A*x=b
%
% A es mxn con rango igual a m
% Q es simétrica positiva definida en el espacio nulo de A

xp = A\b;
Z = null(A);

% Problema cuadrático
% Min (1/2)*xz'*G*xz + op'*xz
G = Z'*Q*Z;
op = Z'*(Q*xp + c);
% Sistema lineal a resolver
% G*xz = -op
xz = -G\op;
x = xp + Z*xz;
%Se calcula el multiplicador de lagrange a partir de la formula
%Q*x+c+A'*lambda=0 (primer condición de KKT)
lambda = -A'\(Q*x+c);
end