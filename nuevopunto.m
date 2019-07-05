function [xnuevo, reac] = nuevopunto(x,p,F,d,ind)
%Calcula el nuevo punto en conjunto activo dependiendo de si el
%conjunto de restricciones activas son todas las restricciones de
%desigualdad o no. Además calcula el nuevo conjunto de restricciones
%activas.
%
% In
% x.- punto actual en R^n
% p.- paso en R^n
% F.- Matriz de restricciones de desigualdad de rxn
% d.- Vector de cotas en R^r
% ind.- Vector de índices.
%
%Out
% alpha.- valor considerado.
%--------------------------------------------------------------------------

%Escalares y vectores constantes 
[r,~] = size(F);
reac = zeros(r,1);
tol = 1.e-12;

if(isempty(ind))
    alpha = 1;
else
    dd = d(ind);
    FF = F(ind,:);
    FFNeg = FF*p;
    vecInd = find(FFNeg < 0);
    if(isempty(vecInd))
       alpha = 1; 
    else
       alphaAux = min((dd(vecInd)-FF(vecInd,:)*x)./(FFNeg(vecInd)));
       alpha = min(alphaAux,1);
    end 
end
xnuevo = x + alpha*p;
%Cálculo del nuevo conjunto de restricciones activas
for i=1:r
    if (abs(F(i,:)*x - d(i)) < tol)
        reac(i) = i;
    end
end
reac(reac == 0) = [];
end