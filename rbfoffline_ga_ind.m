function [Stf,smod,lambda_r,gamma_r] = rbfoffline_ga_ind(vlb_sr, vub_sr, smod, funcobj)

global ac ppc
%
% Geracao do modelo substituto utilizando krigagem
%
% Chamada
%     [dmodel,perf] = rbfoffline(vlb_sr, vub_sr, itersao)
%
% Entrada
%     vlb_sr, vub_sr : limites inferior e superior da subregiao de confianca
%     itersao        : iteracao atual do SAO
%
% Saida
%     dmodel  : modelo de krigagem do DACE, estrutura com os elementos
%       regr   : function handle para o modelo de regressao
%       corr   : function handle para a funcao de correlacao
%       theta  : parametros da funcao de correlacao
%       beta   : estimativa dos minimos quadrados generalisados
%       gamma  : fatores de correlacao
%       sigma2 : estimativa de máxima verossimilhança do processo de variancia
%       S      : amostra em escala
%       Ssc    : fatores de escala para os argumentos da amostra
%       Ysc    : fatores de escala para as imagens da amostra
%       C      : fator deCholesky da matriz de correlacao
%       Ft     : matriz de regressao descorrelacionada
%       G      : Fatorizacao QR: Ft = Q*G' .
%     perf    : estrutura com informacoes da performance do modelo
%       nv     : Numero de avaliacoes da funcao objetivo
%       perf   : vetor (q+2)*nv, onde q eh o numero de elementos em teta
%                e as colunas dao os valores atuaos de
%                [theta;  psi(theta);  type]
%                |type| = 1, 2 or 3, indicando 'start', 'explore' ou 'move'
%                Um valor negativo indica um passo para cima

%Modificado 29/10/2011 -
%Reseta o seed

global pfeas fcount

ndvab = length(vlb_sr);
%%%%%%%%%%%%%%  criando amostra  %%%%%%%%%%%%%%%%%%%%%%%%%
npt =  1 %10*ndvab; % 10n pontos

Sl  = lhsdesign(npt,smod.dim,'criterion','maximin','iteration',10);
%Ajustando amostra para os limites sup e inf
St = ones(npt,1)*vlb_sr + Sl.*(ones(npt,1)*(vub_sr - vlb_sr));
[smod.npt n] = size(St);
% Introducao de pontos viaveis
St_feas = feaspop_lump_solo(vlb_sr,vub_sr,npt);
if pfeas > 0.5
    idSt = randi(1);
    while length(idSt) < round((1-pfeas)*npt)
        compl = round((1-pfeas)*npt) - length(idSt);
        idSt = unique([idSt randi(npt,1,compl)]);
    end
    idSt_feas = setdiff(1:npt,idSt);
else
    idSt_feas = randi(1);
    while length(idSt_feas) < round(pfeas*npt)
        compl = round(pfeas*npt) - length(idSt_feas);
        idSt_feas = unique([idSt_feas randi(npt,1,compl)]);
    end
    idSt = setdiff(1:npt,idSt_feas);
end

Stf = [St_feas(idSt_feas,:); St(idSt,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%  Avaliando a funcao nos pontos dados %%%%%

if numlabs == 1 || npt < numlabs
    for ii=1:npt    
        smod.Y(ii,:) = feval(funcobj,Stf(ii,:));
    end
else
    for nn = 1:(numlabs - 1)                %
        id = numlabs - (nn - 1);            % indica aos outros labs que
        labSend(1,id);                      % havera funcao a ser executada
        labSend(1,id);                      %
    end                                     %
    clear pops ppl xs ppl_id init fin
    pops = size(Stf,1);
    ppl = fix(pops/numlabs);
    xs = rem(pops,numlabs);
    ppl_id = ppl*ones(1,numlabs);
    ppl_id(1:xs) = ppl + 1;
    init = ones(1,numlabs);
    fin = init*(ppl+1) - 1;
    for nn = 1:(numlabs - 1)
        id = numlabs - (nn - 1);
        init(id) = 1 + sum(ppl_id(1:id-1));
        fin(id) = init(id) + ppl_id(id) - 1;
        labSend(Stf(init(id):fin(id),:),id);
        labSend(func2str(funcobj),id);
    end
    for ii = init(1):fin(1)
        smod.Y(ii,:) = feval(funcobj,Stf(ii,:));
    end
    for nn = 1:(numlabs - 1)
        id = numlabs - (nn - 1);
        data_sqp = labReceive(id);
        init(id);
        fin(id);
        smod.Y(init(id):fin(id),:) = data_sqp;
    end
end

fcount = fcount + size(smod.Y,1);
Y = smod.Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%  Criando Modelo Aproximado %%%%%%%%%%%%

rbftype = smod.reg;

for ii=1:size(Y,2)
    [lambda1, gamma1]=RBF(Stf,Y,rbftype);
    lambda_r(:,ii)=lambda1;
    gamma_r(:,ii)=gamma1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Salvando o modelo %%%%%%%%%%%%%%%%%%%
save tape1g Stf Y lambda_r gamma_r smod  % Usa em myfun_krig