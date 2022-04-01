function [A_ind, b_ind] = Abcc2Abind(A,b,vdg_rpl,npp,npi,ntcp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRICAO:
% Transforma a matriz A e o vetor b, provenientes do 'lincons.m' e baseado
% em ciclo de controle, em nova matriz A e vetor b baseado em pocos 
% independentes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DADOS:
%   A_ind: matriz de restricoes dependente das variáveis no caso de pocos 
%          independentes                                              (out)
%   b_ind: vetor de restricoes independente das variaveis no caso de pocos 
%          independentes                                              (out)
%       A: matriz de restricoes dependente das variáveis no caso de pocos 
%          com ciclos de controle                                      (in)                              
%       b: VETOR de restricoes independente das variaveis no caso de pocos 
%          com ciclos de controle                                      (in) 
% vdg_rpl: flag com valor de 0 ou 1 para consideracao do voidage 
%          replacement, 0 = NAO cansidera , 1 = considera              (in)
%     npp: numero de pocos produtores                                  (in)
%     npi: numero de pocos injetores                                   (in)
%    ntcp: numero de trocas de cada poco                               (in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBSERVACOES:
% - Esse programa tem o intuito de crias a matriz de restricoes de um dado
%   problema, para a otimizacao baseada em pocos independentes seguindo a
%   formulacao abaixo:
%       A_ind*x_ind = b_ind
% onde:
%       A_ind e b_ind = sao calculados aqui.
%               x_ind = deve ser calculado em 'xcc2xind.m'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTOR:
% Danilo Rudrigues de Almeida Lira
% VERSAO: 2017.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quantidade de trocas totais:
qdtt = (npp+npi)*ntcp;

if vdg_rpl == 0;
   A_ind = zeros(((qdtt+1)*2)+(npp+npi),((npp+npi)*(qdtt+1))+qdtt);
   b_ind = zeros(((qdtt+1)*2)+(npp+npi),1); 
   Acc(1:2,1:(npp+npi)) = A(1:2,1:(npp+npi)); 
   bcc(1:2,1) = b(1:2,1);
   for i = 1:(qdtt+1);
       A_ind((2*i)-1:2*i,(i*(npp+npi))-(npp+npi-1):i*(npp+npi)) = Acc;
       b_ind((2*i)-1:2*i,1) = bcc;
   end
else
   A_ind = zeros(((qdtt+1)*3)+(npp+npi),((npp+npi)*(qdtt+1))+qdtt);
   b_ind = zeros(((qdtt+1)*3)+(npp+npi),1);  
   Acc(1:3,1:(npp+npi)) = A(1:3,1:(npp+npi)); 
   bcc(1:3,1) = b(1:3,1);
   for i = 1:(qdtt+1);
       A_ind((3*i)-2:3*i,(i*(npp+npi))-(npp+npi-1):i*(npp+npi)) = Acc;
       b_ind((3*i)-2:3*i,1) = bcc;
   end
end
 v_ones = ones(1,(npp+npi));
rest_tempo = diag(v_ones);
for ii = 1:ntcp
    A_ind((end-(npp+npi)+1):end,(end-qdtt+1+((ii-1)*(npp+npi))):(end-qdtt+((ii-1)*(npp+npi))+(npp+npi))) = rest_tempo;
end
b_ind((end-(npp+npi)+1):end,1) = 1;
A_ind = sparse(A_ind);
end