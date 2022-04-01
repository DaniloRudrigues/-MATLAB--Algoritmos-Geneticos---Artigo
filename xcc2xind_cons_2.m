function x_ind = xcc2xind_cons_2(x,npp,npi,ntcp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRICAO:
% Transforma populção de individuos baseado em ciclo de controle, em
% populacao de individuos baseado em pocos independentes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DADOS:
% x_ind: matriz de individuos para otimizacao com pocos independentes (out)
%     x: matriz de individuos para otimizacao com ciclos de controle   (in)
%   npp: numero de pocos produtores                                    (in)
%   npi: numero de pocos injetores                                     (in)
%  ntcp: numero de trocas de cada poco                                 (in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBSERVACOES:
% - O numero de ciclos de controle e equivalente ao numero de vazoes
%   normalizadas que cada poco vai trocar durante a otimizacao, ou seja, o
%   numero de troca de cada poco (ntcp) e igual a (ncc-1).
%   ntcp = ncc-1
% - e importante alterar o IMEX para se adaptar ao conceito de ntcp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTOR:
% Danilo Rudrigues de Almeida Lira
% VERSAO: 2017.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = x;
linhas_x1 = size(x1,1);
qdtt = (npp+npi)*ntcp; 
x_ind = zeros(linhas_x1,(npp+npi)+((npp+npi)*qdtt)+qdtt);

for i = 1:linhas_x1;
    m = zeros((npp+npi)*ntcp,3);
    % tempos:
    m(:,1) = x1(i,((npp+npi)*(ntcp+1)+1):end);   
    
    % indice do poco que vai ser trocado:
    for o = 1:(qdtt/(npp+npi));
        m(((o-1)*(npp+npi))+1:(o-1)*(npp+npi)+(npp+npi),2) = 1:(npp+npi);
    end
    
    % novo valor do poco
    m(:,3) = x1(i,(npp+npi)+1:((npp+npi)+1)+((npp+npi)*ntcp-1)); 
    
    % Ordenacao crescente:
    M = sortrows(m);
    
    % criacao da matriz de trocas:
    M_troc = zeros(qdtt,(npp+npi));
    for o = 1:qdtt;
        M_troc(o,:) = x1(i,1:(npp+npi));
    end
    for o = 1:qdtt;
        M_troc(o,M(o,2)) = M(o,3);
        M_troc(o+1,:) = M_troc(o,:);
    end
    M_troc(end,:) = [];
    d = M(:,1)';
    sd = size(d,2);
    M_troc2 = M_troc;
    for ik = 1:sd
        id = find(d==d(sd-ik+1));
        sid = size(id,2);
        for i2 = 1:sid
             M_troc2(id(sid-i2+1),:) = M_troc(id(sid),:);
        end
    end
    % criacao do x_ind:
    
        % primeira parte do x_ind:
    x_ind(i,1:(npp+npi)) = x1(i,1:(npp+npi));
        % segunda parte do x_ind:
    for o = 2:qdtt+1;
        x_ind(i,((o-1)*(npp+npi))+1:(o-1)*(npp+npi)+(npp+npi)) = M_troc2(o-1,:);
    end
        % ultima parte do x_ind:
    %x_ind(i,(npp+npi)+((npp+npi)*qdtt)+1:end) = M(:,1);
    x_ind(i,(npp+npi)+((npp+npi)*qdtt)+1:end) = x(i,((npp+npi)*(ntcp+1)+1):end);
end



end