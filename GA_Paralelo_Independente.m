%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%% 
%%%%                      GA_Paralelo_Independente                     %%%%
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
%           - otimização em paralelo no Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTOR:
%       Danilo Rudrigues de Almeida Lira
% VERSAO:                                                     DATA:
%       2017.1                                                      09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBSERVACOES:
%       - No arquivo '.in' o NCC vira NTCP, ou seja, ntcp = ncc-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc
global id_apm A b fobj cons_bnd lb ub vdg_rpl npp npi ntcp ac infile std_fit scale shrink
myname = labindex;
tempo_inicial = clock;
%%%%%%%%%%%%%%%%%%%%%%%%%% Dados de entrada: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ativar o método de penalizacao adaptativa: (1-ativa ; 0-nao ativa)
id_apm = 1;
% Ativar restrições dos limites: (1-ativa ; 0-nao ativa)
cons_bnd = 1;
% Ativar restricao de conservacao de energia: (1-ativa ; 0-nao ativa)
vdg_rpl = 1;
% limites:
   % inferior
lb = [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001];
   % superior
ub = [0.75 0.75 1 0.75 0.75 1 0.75 0.75 1 0.6 0.6 0.6 0.6 0.6 0.6];
% Dados do reservatorio:
infile = 'D1_9q6t_ind.in';     
    ac = ac_data(infile);
    ac.delta = 1.1;
    npp = ac.npp;
    npi = ac.npi;
    ntcp = ac.ncc;
    ac.vdg_rpl = vdg_rpl;
% matriz de restricoes baseado em ciclo de controle:
[A, b] = lincons_mod(ub);
if myname == 1
% funcao objetivo
fobj = @fvpl_ind;
x_med = (ub + lb)/2;
std_fit = 1;
std_fit = feval(fobj,x_med);
for nn = 1:(numlabs - 1)
        id = numlabs - (nn - 1);
        labSend(std_fit,id);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%           GA           %%%%%%%%%%%%%%%%%%%%%%%%%
nv = length(lb);
ps = 100;
generation = 100;
parcela_dominio = 1;
CrossoverFraction = 0.2;
PopInitRange = [0.001;0.200];
pt_ex = [0.54 0.46 0.926 0.092 0.076 0.152 0.082 0.616 0.665 0.432 0.432 0.432 0.546 0.546 0.546;...
         0.54 0.46 0.926 0.092 0.076 0.152 0.082 0.616 0.665 0.132 0.232 0.432 0.746 0.646 0.546];
popi = [pt_ex;feaspop_ind(lb,ub,ps-2,ac,parcela_dominio)];
scale = 0.05;
shrink = 0.5;
   opti = gaoptimset('CrossoverFraction', CrossoverFraction,...
        'PopInitRange', PopInitRange,...
        'Generation', generation,...
        'InitialPopulation', popi,... 
        'PopulationSize', ps,...
        'MutationFcn',  @mut_gaussiana_ind,...{@mutationgaussian,0.05,0.5},...,...@mutationuniform,...  
        'Vectorized', 'on');
[x_ot,fit,exitflag,output] = ga(@f_ga_SP_ind_PAR_3,nv,[],[],[],[],[],[],[],opti);
%%%%%%%%%%%%%%%%%%%%%%%%%%%          VPL          %%%%%%%%%%%%%%%%%%%%%%%%%
VPL = std_fit*feval(fobj,x_ot)
x_ind_ot = xcc2xind(x_ot,npp,npi,ntcp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%   pos processamento   %%%%%%%%%%%%%%%%%%%%%%%%% 
for nn = 1:(numlabs - 1)
        id = numlabs - (nn - 1);
        labSend(0,id);
end
tempo_final = clock;
tempo_total = etime(tempo_final,tempo_inicial);
name = ['VPL_' num2str(VPL) '_tempo_' num2str(tempo_total) '.mat'];
save(name);
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
end    
disp('__________________Fim do Processo__________________')