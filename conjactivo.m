function [x, mu, iter] = conjactivo(Q,F,c,d)
% Metodo de conjunto activo para el problema cuadrático
% Min   (0.5)* x' * Q * x + c´* x
% s.a.  F*x >= d
%
%  Llamado: [x, mu] = conjactivo(Q,F,c,d);
% In
% Q.- matriz nxn simétrica y positiva definida
% F.- matriz rxn.
% d.- vector columna en R^r .
% c.- vector columna en R^n .
%
%Out
% x.- vector en R^n con la aproximación del mínimo local.
% mu.- vector en R^r con la aproximación al multiplicador
%          de Lagrange asociado a la restricción de desigualdad.
% iter.- número de iteraciones.
%--------------------------------------------------------------------------

%Definción de constantes y de variables
r = length(d);
tol = 1.e-06; %Tolerancia para la norma del vector de paso
flag = 1; %con esta variable se controla el while
maxiter = 100;
iter = 0;
I = [];
mu = [];
eta = [];

%Se calcula si el conjunto factible es vacío o no
[x,reac] = factible(F,d);

if (isempty(x))%Es vacío
    fprintf('El conjunto factible es vacío')
else %No es vacío
    mu = zeros(r,1);
    while(flag == 1 && iter < maxiter)
        if (isempty(I)) %I es vacío
            %El subproblema es
            % min (1/2)p'*Q*p + (Q*x+c)'*p
            %Se resuelve derivando la función objetivo e igualando a cero
            aux = -Q\c;
            p = aux - x;
        else
            %Se resuelve el subproblema
            % min (1/2)p'*Q*+gradq'*p
            % s.a FF*x=0   esta FF es la F con los índices en I
            gradq = Q*x + c; %Gradiente de la función objetivo
            FF = F(I,:);
            [p, eta] = pcnulo(Q,FF,gradq,0); %Método del espacio nulo
        end
        if(norm(p) > tol)
            if(~isempty(reac))
                % Se calcula si existe j en reac-I tal que f_j'*p < 0
                VecDiff = setdiff(reac,I); %reac-I
                noreac = setdiff(1:r,reac); %[1,...,r]-reac
                if(~isempty(VecDiff))
                    %Se hace la búsqueda de algún producto F_j'*p <0
                    naux = 0;
                    contaux = 0;
                    while(contaux < length(VecDiff) && naux == 0)
                        contaux = contaux + 1;
                        if(F(VecDiff(contaux),:)*p < 0)
                            naux = VecDiff(contaux);
                        end
                    end
                    if(naux > 0) %Si existió j en reac-I tal que el producto punto fuera <0
                        I = [I;naux]; 
                    else %Para toda j en reac-I el prodcuto punto fue >=0
                        [x, reac] = nuevopunto(x,p,F,d,noreac);
                    end
                else
                    [x, reac] = nuevopunto(x,p,F,d,noreac);
                end
            else %No hay restricciones activas
                [x, reac] = nuevopunto(x,p,F,d,1:r);
            end
        else
            if(isempty(eta)) %El paso p fue cero con I siendo vacío
                flag = 0; %Se sale del while. Punto óptimo en el interior
            else %I no fue vacío
                if(all(eta<=0)) % Se tiene un mínimo local
                    mu(I) = -eta;
                    flag = 0; %Se sale del while
                else %Se busca el min índice para el cual eta tiene una entrada positiva
                    indPos = find(eta > 0,1);
                    I(I==indPos) = [];% I <- I - {indPos}
                end
            end
        end
        iter = iter + 1;
    end
end
end



