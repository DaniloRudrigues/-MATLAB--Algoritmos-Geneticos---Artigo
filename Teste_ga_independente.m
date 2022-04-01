%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%% 
%%%%                         GA_Independente                           %%%%
%%%%                                                                   %%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRICAO: 
%           Projeto de otimizacao de VPL simplificado para reservatorios
%           de petroleo, considerando:
%           - Algoritmo genetico (GA).
%           - Restricoes LINERES.
%           - Pocos independentes (trocas de vazoes independentes por poco)
%           - Metodo da penalizacao adaptativa (metodo da barreira)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTOR:
%       Danilo Rudrigues de Almeida Lira
% VERSAO:                                                     DATA:
%       2017.3 - SOMATORIA DOS TEMPOS DE CADA POCO <= 1             07/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBSERVACOES:
%       - No arquivo '.in' o NCC vira NTCP, ou seja, ntcp = ncc-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear,clc
global ac infile std_fit id_apm A b fobj cons_bnd lb ub vdg_rpl npp npi ntcp
%--------------------------------------------------------------------------
% 1 -  PRE-PROCESSAMENTO
%--------------------------------------------------------------------------
% 1.1 - DADOS IMPORTADOS/EXPORTADOS:
infile = 'D1_9q6t_ind.in';     
    ac = ac_data(infile);
    ac.delta = 1.1;
    
% 1.2 - DADOS DE ENTRADA:
    % 1.2.1 - FLAGS [zero(0) ou um(1)]:
  id_apm = 1;
cons_bnd = 0;
 vdg_rpl = 1;
    % 1.2.2 - FORNECIDOS:
    fobj = @fvpl_ind; 
      lb = [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001];
      %ub = [0.75 0.75 1 0.75 0.75 1 0.75 0.75 1 0.6 0.6 0.6 0.6 0.6 0.6];
      ub = [0.5 0.5 1 0.5 0.5 1 0.5 0.5 1 0.5 0.5 0.5 0.5 0.5 0.5];
      ps = 10; 
    [A, b] = lincons_mod(ub);
    fat_viavel = 1;
    popi = feaspop_ind(lb,ub,ps,ac,fat_viavel);     
    opti = gaoptimset('CrossoverFraction', 0.2,...
        'Generation', 5,...
        'InitialPopulation', popi,... 
        'PopulationSize', ps,...
        'MutationFcn', {@mutationadaptfeasible},... @mutationuniform,... %{@mutationgaussian,0.05,0.5},...
        'Vectorized', 'on');
    

    % 1.3 - DADOS CALCULADOS:
nv = length(lb);
npp = ac.npp;
npi = ac.npi;
ntcp = ac.ncc;
ac.vdg_rpl = vdg_rpl;
x_med = (ub + lb)/2;
std_fit = 1;
std_fit = feval(fobj,x_med);
%__________________________________________________________________________


%--------------------------------------------------------------------------
% PROCESSAMENTO
%--------------------------------------------------------------------------
[x_ot,fit,exitflag,output] = ga(@f_ga_ind,nv,[],[],[],[],lb,ub,[],opti);
%__________________________________________________________________________


%--------------------------------------------------------------------------
% POS - PROCESSAMENTO
%--------------------------------------------------------------------------
VPL = std_fit*feval(fobj,x_ot);
x_ind_ot = xcc2xind(x_ot,npp,npi,ntcp);
%__________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('__________________Fim do Processo__________________')
toc