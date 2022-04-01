% restricao de tempo nova 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTE 1 - NOVAS RESTRICOES DE TEMPO.
% DESCRICAO: verificar se os tempos de troca de cada poco estao
%            sequenciados,por exemplo, verificar se o valor inicial de um 
%            poco nao vai
%            para o valor numero 5, ao inves de ir para o numero 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ntcp ~= 1;
   t = zeros(ntcp,(npp+npi));
   for i = 1:linhas_x1;
    for k = 1:ntcp;
      t(k,:) = x(i,((k-1)*(npp+npi))+((npp+npi)*(ntcp+1)+1):...
          ((k-1)*(npp+npi))+((npp+npi)*(ntcp+1)+1)+(npp+npi-1));
    end  
      T = sort(t);
    if t ~= T;
       x1(i,:) = 0;  % zerar o individuo
    end
   end
end