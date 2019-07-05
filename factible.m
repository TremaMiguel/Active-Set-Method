function [x,reac] = factible(F,d)
% Problema de programación lineal para encontrar un punto factible de
%  F*x>=d
%
% Se resuleve el problema
% Min 0'*x + e'*z
% s.a F*x + z = d
%           z >= 0
%In
% F.- matriz de orden rxn.
% d.- vector columna de orden r.
%
%Out
% x.- vector columna de tamaño n, punto factible del conjunto.
% reac.- vector con las restricciones activas de F*x - d >= 0
%       reac(i) ==   0    indica que la restricción #i no es activa en x
%       reac(i) ==  i   inidca que la restricción # i es activa en x.
%-----------------------------------------------------------------------

%Definición de función obejtivo y de las restricciones
[r,n] = size(F);
f = [zeros(n,1); ones(r,1)];
Ain = [-F -eye(r)]; bin = -d;
A = zeros(1,n+r); b = 0;
lb = [-inf*ones(n,1);zeros(r,1)];
ub = inf*ones(n+r,1);

%Se resuelve el problema de P.L.
[w,fw] = linprog(f,Ain,bin,A,b,lb,ub);

%Inicialización de los vectores respuesta
x = [];
reac = zeros(r,1);

%Resultados de la solución del problema de P.L.
if(abs(fw) < 1.e-06) %Conjunto no vacío     %MODIFICARLO
    x = w(1:n);
    for i=1:r
        if (abs(F(i,:)*x - d(i)) < 1.e-12)
            reac(i) = i;
        end
    end
end
reac(reac == 0) = [];
end