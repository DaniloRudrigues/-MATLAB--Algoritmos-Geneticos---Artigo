%###############################################################%
%#                                                             #%
%#                    Gglobal search - GA                      #%
%#                                                             #%
%###############################################################%

% Algoritmo genetico para problemas de gerenciamento da producao
% de reservatorios de petroleo
% versao - 2014.1

clear variables global
clc

global ac infile fcount pfeas id_apm cons_bnd A b std_fit fobj fcons lb ub ntcp npp npi vdg_rpl id_hyb scale shrink ppc
 
myname = labindex;

addpath('dace','RBF');

% Dados de entrada do problema:

% 1 - cota das variaveis

lb = 0.001*ones(1,12);
ub = [0.18 0.18 0.18 0.18 0.18 0.18 0.18 1500/5750 1500/5750 1500/5750 1500/5750 1500/5750];     % inicial
     % 0.75 0.75 0.75 0.75 0.75 0.75 0.75 1 1 1 1 1 ...            trocas 1
     % 0.75 0.75 0.75 0.75 0.75 0.75 0.75 1 1 1 1 1 ...            trocas 2
     % 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6...          tempos 1
      % 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6];         % tempos 2
dim = length(lb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2 - funcao objetivo e restricoes

% 2.1* - caso particular de otimizacao de reservatorio

infile = 'BCO_nt_12q0t_ind.in';     % arquivo com os dados da otimizacao
ac = ac_data(infile);
npp = ac.npp;
npi = ac.npi;
ntcp = ac.ncc;
if ac.opera == 1          % definicao do parametro de "overinjection" para
    ac.delta = 1.15;      % manutencao de pressao no reservatorio
end                       %

obj = @fvpl_ind;
fcount = 0;

ac.vdg_rpl = 0;           % identifica a consideracao do voidage replacement
                          % como restricao na oimizacao
                          % 0 - nao considera
                          % 1 - considera
vdg_rpl = ac.vdg_rpl;
[A, b] = lincons_mod(ub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3 - modelos substitutos

% 3.1 - tipo de modelo considerado

sur_model = 4;   % 0 - high fidelity
                 % 1 - kriging como regpoly 0
                 % 2 - kriging como regpoly 1
                 % 3 - RBF linear
                 % 4 - RBF cubico
                 % 5 - RBF tps

switch sur_model
    case 0
        fobj.ga = obj;
    case {1,2}
        fobj.ga = @myfun_krig_ga;
    case {3,4,5}
        fobj.ga = @myfun_rbf_ga;
end

if sur_model ~= 0
    smod = data_smod(length(lb),sur_model);
end

if myname == 1

% 3.2 - atualizacao do modelo substituto para busca global

    npts_up = 5;       % numero de pontos acrescentados a amostra inicial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4 - configuracao dos otimizadores

% 4.1 - busca global - (GA)

    
% %  GA settings
% %
% %     'InitialPopulation': matrix
% %     'PopulationType': 'bitstring'
% %     'FitnessScalingFcn': {@fitscalingshiftlinear,n}
% %                          @fitscalingprop
% %                          {@fitscalingtop,n} <--- no funcionou
% %     'SelectionFcn': @selectionremainder
% %                     @selectionuniform
% %                     @selectionroulette
% %                     {@selectiontournament,size}
% %     'EliteCount',n
% %     'CrossoverFraction',n
% %     'CrossoverFcn': @crossoversinglepoint
% %                     @crossovertwopoint
% %                     {@crossoverintermediate,n}
% %                     {@crossoverheuristic,n}
% %                     @crossoverarithmetic
% %     'MutationFcn': {@mutationuniform,n}
% %                    @mutationadaptfeasible

% as opcoes estao definidas posteriormente devido a criacao da populacao
% inicial diferenciada!

    id_apm = 1;     % indica se o metodo apm e utilizado
                    % 0 - nao
                    % 1 - sim

    cons_bnd = 1;   % indica se as restricoes de limites das varizaveis sao 
                    % consideradas como restricao do problema
                    % 0 - nao
                    % 1 - sim

    pfeas = 0.6;    % indica o percentual de cromossomos viaveis na criacao do
                    % substituto (abs: 1 = 100%)
    ps = 100;
    std_fit = 1;
    x_med = (ub + lb)/2;
    std_vpl = feval(obj,x_med);
    std_fit = std_vpl;
    for nn = 1:(numlabs - 1)
        id = numlabs - (nn - 1);
        labSend(std_fit,id);
    end
    fcount = fcount + 1;

    if sur_model ~= 0
        switch sur_model
            case {1, 2}
                [dmodel, perf] = krigoffline_ga(lb,ub,smod,obj);
            case {3, 4, 5}
                [Stf,smod,lambda_r,gamma_r] = rbfoffline_ga_ind(lb,ub,smod,obj);
        end
        load tape1g
        switch sur_model
            case {1, 2}
                save tape2g Stf Y dmodel perf smod
            case {3, 4, 5}
                save tape2g Stf Y lambda_r gamma_r smod
        end
        doe0 = [Stf Y];
        poolfeas = feaspop(lb,ub,round(0.3*ps));
        maxpt = 2*length(Stf);
        solfound = 0;
        positsf = [];
        fstr_old = [];
    end
    % ps = 10;               % tamanho da populacao
    generation = 50;
    PopInitRange = [0.001;0.200];
    CrossoverFraction = 0.2;
    scale = 1;
    shrink = 0.5;
    ppc = 1;
    conv = 0;
    %tent = [0.486100000000000 0.493200000000000 0.891400000000000 0.434300000000000 0.447300000000000 0.812500000000000 0.445600000000000 0.423881250000000 0.819175000000000 0.567400000000000 0.241300000000000 0.365100000000000 0.214100000000000 0.486200000000000 0.455326953125000];
    %tent2 = [0.5309 0.4072 0.9365 0.5468 0.3339 0.8532 0.5468 0.3339 0.8532 0.3040 0.1701 0.4092 0.3940 0.4920 0.1295];
    popi = feaspop(lb,ub,ps);
    optionsH = gaoptimset('CrossoverFraction', CrossoverFraction,...
        'PopInitRange', PopInitRange,...
        'Generation', generation,...
        'InitialPopulation', popi,...
        'PopulationSize', ps,...
        'MutationFcn', @mut_gaussiana_ind_2,...{@mutationgaussian,scale,shrink},...
        'Vectorized', 'on');
    
    if sur_model == 0
        [xstr_ga,fit,exitflag,output] = ga(@f_ga_SP_ind_PAR_2,length(lb),[],[],[],[],[],[],[],optionsH);
    else
        while conv == 0 && isempty(poolfeas) == 0
            if exist('xstr_ga')
                popi = [feaspop(lb,ub, ps-1); xstr_ga];
            end
            
            [xstr_ga,fit,exitflag,output] = ga(@f_ga_SP_ind_PAR_2,length(lb),[],[],[],[],[],[],[],optionsH);
            
            fsbo = fit;
            if isempty(xstr_ga)
                xsurup = zeros(npts_up,dim);
                fsurup = zeros(npts_up,size(Y,2));
                for pp = 1:npts_up
                    if isempty(poolfeas)
                        fprintf('\n The pool of feasible samples is empty!\n\n');
                        xsurup(pp:end,:) = [];
                        fsurup(pp:end,:) = [];
                        break
                    end
                    iisel = randi(size(poolfeas,1));
                    xsurup(pp,:) = poolfeas(iisel,:);
                    poolfeas(iisel,:) = [];
                    fsurup(pp,:) = feval(obj, xsurup(pp,:));
                    fcount = fcount + 1;
                end
                xnew = xsurup(end,:);
                xsurup(end,:) = [];
            else
                solfound = solfound + 1;
                positsf = [positsf (size(Stf,1)+1)];
                xnew = xstr_ga;
                xsurup = zeros(npts_up-1,dim);
                fsurup = zeros(npts_up-1,size(Y,2));
                for pp = 1:npts_up-1
                    if isempty(poolfeas)
                        fprintf('\n The pool of feasible samples is empty!\n\n');
                        xsurup(pp:end,:) = [];
                        fsurup(pp:end,:) = [];
                        break
                    end
                    iisel = randi(size(poolfeas,1));
                    xsurup(pp,:) = poolfeas(iisel,:);
                    poolfeas(iisel,:) = [];
                    fsurup(pp,:) = feval(obj, xsurup(pp,:));
                    fcount = fcount + 1;
                end
                fstr = feval(obj, xnew);
                fcount = fcount + 1;
                if isempty(fstr_old)
                    fstr_old = fstr;
                else
                    if fstr == fstr_old
                        xnew = [];
                        fstr = [];
                        solfound = solfound - 1;
                        positsf(end) = [];
                    else
                        fstr_new = min([fstr_old fstr]);
                        fstr_old = fstr;
                        fstr = fstr_new;
                    end
                end
            end
            if abs((fstr - fsbo)/fstr) < 1e-6 % se n convergir troca o modelo
                conv = 1;
            else
                Stf = [Stf; xnew; xsurup];
                Y = [Y; fstr; fsurup];
            end
            
            switch sur_model
                case {1, 2}
                    Stnew = Stf;
                    Ynew = Y;
                    if smod.lb == 0 & smod.ub == 0
                        [dmodel, perf] = dacefit(Stnew,Ynew,smod.reg,smod.cor,smod.theta);
                    else
                        [dmodel, perf] = dacefit(Stnew,Ynew,smod.reg,smod.cor,...
                            smod.theta,smod.lb,smod.ub);
                    end
                    save tape2g Stf Y dmodel perf smod
                case {3, 4, 5}
                    Yrbf = Y;
                    rbftype = smod.reg;
                    [lambda_r, gamma_r] = RBF(Stf,Yrbf,rbftype);
                    save tape2g Stf Y lambda_r gamma_r smod
            end
            save backup_temp.mat
        end
    end
    
    fstr_ga = std_fit*feval(obj,xstr_ga)
    fcount = fcount + 1;
    
    for nn = 1:(numlabs - 1)
        id = numlabs - (nn - 1);
        labSend(0,id);
    end
    
    if ac.opera == 0
        opera = 'top';
    else
        opera = 'ntop';
    end
    
    if ac.time == 0
        tempo = 'tcte';
    else
        tempo = 'tvar';
    end
    d_rod = clock;
    comp_name = ['d' num2str(d_rod(3)) 'm' num2str(d_rod(2)) 'a' num2str(d_rod(1)) 'h' num2str(d_rod(4)) '_' num2str(d_rod(5))];
    name = ['res_hyb_' num2str(ac.npp) 'p' num2str(ac.npi) 'i_' num2str(ac.ncc) 'cc_' opera '_' tempo '_' comp_name '.mat'];
    save(name)
    fprintf('Optimization completed! See file: ')
    disp(name)
else
    std_fit = labReceive(1);
    id_run = labReceive(1);
    while id_run == 1
        clear id_fun x_test f_par fobj_par f_test
        id_fun = labReceive(1);
        x_test = labReceive(1);
        f_par = labReceive(1);
        fobj_par = str2func(f_par);
        for ii = 1:size(x_test,1)
            f_test(ii,:) = feval(fobj_par,x_test(ii,:));
        end
        labSend(f_test,1);
        id_run = labReceive(1);
    end
    fprintf('Optimization completed! See Lab 1\n ')
end