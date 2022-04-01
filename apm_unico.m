function [fit, idp] = apm_unico(f, cons)

% apm - Adaptive Penalty Function
% 
% fit - funcao penalizada
% idp - identificador de panalizacao (identifica se o individuo foi
%       penalizado ou nao
% 
% f - valor da funcao objetivo (a ser normalizado)
% cons - valor das restricoes (individuo = coluna)

pops = length(f);
fit = zeros(size(f));
idp = fit; multi = 1;

fmed = mean(f);
viola = cons;
viola(cons<0) = 0;
vmed = mean(viola,2);     % media de todas as restricoes
sqvmed = sum(vmed.^2);    % denominador do parametro de penalizacao
if sqvmed == 0
    sqvmed = 1;
    multi = 0;
end

penal = multi*abs(fmed)*vmed/sqvmed;  % parametro de penalizacao

for jj = 1:pops
    if find(cons(:,jj)>0)
        if f(jj) > fmed
            fcor = f(jj);
        else
            fcor = fmed;
        end
        fit(jj) = fcor + 1000000;
        % fit(jj) = fcor + viola(:,jj)'*penal;
        idp(jj) = 1;
    else
        fit(jj) = f(jj);
    end
end