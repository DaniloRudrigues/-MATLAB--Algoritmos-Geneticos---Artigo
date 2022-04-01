function [tags] = ManipParamNp_indep(x_manip,ac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRICAO:
%  Cria os cards de insercao para o IMEX. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DADOS:
%    tags: estrutura de valores para ser inserido no IMEX             (out)
%      ac: estrutura de dados coletados no arquivo '.in'               (in)
% x_manip: vetor linha igual ao xcc do 'xcc2xind.m'             (in-global)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBSERVACOES:
% - O numero de ciclos de controle e equivalente ao numero de vazoes
%   normalizadas que cada poco vai trocar durante a otimizacao, ou seja, o
%   numero de troca de cada poco (ntcp) e igual a (ncc-1).
%   ntcp = ncc-1
% - Para utilizacao do programa, no arquivo '.in', no lugar de inserir o
%   'ncc' deve-se inserir o 'ntcp' que deve ter o valor de (ntcp = ncc-1),
%   exemplo: se o 'ncc' planejado para o problema for 5, deve-se inserir o
%   valor 4, pois ntcp = ncc-1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTOR:
% Danilo Rudrigues de Almeida Lira
% VERSAO: 2017.2 - TEMPO NORMAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npp = ac.npp;	
npi = ac.npi;
ntcp = ac.ncc;   % veja as OBSERVACOES
Tc = ac.T;
Mqp = ac.Qp;
Mqi = ac.Qi;
tags = struct([]);
x_tempo = zeros(1,((npp+npi)*(ntcp+1)));
x_tempo(1,((npp+npi)+1):end) = x_manip(1,((npp+npi)*(ntcp+1)+1):end);

% ordenação do x_manip
x_ord = zeros(size(x_manip));
x_ord_tempo = zeros(size(x_tempo));
x_ord(1,1:(npp+npi)) = x_manip(1,1:(npp+npi));
pc = zeros(ntcp,ntcp);
for ip = 1:(npp+npi);
    for it = 2:ntcp+1;
        pc((it-1),1) = x_tempo(((npp+npi)*(it-1))+ip);
        pc((it-1),2) = x_manip(((npp+npi)*(it-1))+ip);
    end
    pc_ord = sortrows(pc);
    for it = 2:ntcp+1;
              x_ord(((npp+npi)*(it-1))+ip) = pc_ord((it-1),2);
        x_ord_tempo(((npp+npi)*(it-1))+ip) = pc_ord((it-1),1);
    end
end
x_ord(((npp+npi)*(ntcp+1))+1:end) =  x_ord_tempo(npp+npi+1:end);
x_tempo_acumulado = zeros(size(x_ord_tempo));
for ita = (npp+npi+1):size(x_ord_tempo,2);
    x_tempo_acumulado(ita) = x_ord_tempo(ita) + x_tempo_acumulado(ita-(npp+npi));
end


% criacao dos tags para o IMEX
for it = 1:ntcp+1;
    for ip = 1:npp;
        ivarp = ip + (npp + npi)*(it-1);
        vazao_absoluta = x_manip(ivarp)*Mqp;
        tags(ivarp).name = sprintf('$GP%i_%i',it,ip);
        tags(ivarp).val = vazao_absoluta;
        tags(ivarp).type = 1;      % ---> t_var
        tags(ivarp).number = ip;   % (p_prod)
        tags(ivarp).time = Tc*365*x_tempo_acumulado(ivarp);
    end        

    for ip = 1:npi;
        ivari = ip + (npp + npi)*(it-1) + (npp);
        vazao_absoluta = x_manip(ivari)*Mqi;
        tags(ivari).name = sprintf('$GI%i_%i',it,ip);
        tags(ivari).val = vazao_absoluta;
        tags(ivari).type = 2;      % ---> t_var
        tags(ivari).number = ip;   % (p_inj)
        tags(ivari).time = Tc*365*x_tempo_acumulado(ivari);
    end
end
end