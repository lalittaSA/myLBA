% adapted from Danial's LBA
% created 20170106


% TODO organize recorded data & run through
clear; clc; close all;

tic
%% set paths
script_path = '/modelling/myLBA/';
data_path = '/dataAnalysis/data/';

filename = '';

%% parameters
RT_bins = 0:30:1000;
threshold = 1;
max_it = 1000;
noise_sd = 0; % mean noise = 0
leakageCoef = 0;

Sim_nTrials = 504; % cond*n
freeParams= {'V1', 'V2', 'Sp1', 'Sp2','B12', 'B21'};
fixedParams= {'Dn', 'Sv1', 'Sv2'};
Prim_Sim_data= {NaN, NaN, NaN};
Opt_init_val= nan(1, size(PUBPAR.FreeParams, 2));

%% preallocating a big matrix (list) of all parameter-sets

range_of_parameters= [
    1.50e-3, 0.35e-3 ,4.50e-3;
    1.50e-3, 0.35e-3 ,4.50e-3;
    100e-3, 100e-3 , 500e-3;
    100e-3, 100e-3 , 500e-3;
    0.0e-3, 0.3e-3 ,4.00e-3;
    0.0e-3, 0.3e-3 ,4.00e-3];

A= range_of_parameters(1,1):range_of_parameters(1,2): range_of_parameters(1,3);
B= range_of_parameters(2,1):range_of_parameters(2,2): range_of_parameters(2,3);
C= range_of_parameters(3,1):range_of_parameters(3,2): range_of_parameters(3,3);
D= range_of_parameters(4,1):range_of_parameters(4,2): range_of_parameters(4,3);
E= range_of_parameters(5,1):range_of_parameters(5,2): range_of_parameters(5,3);
F= range_of_parameters(6,1):range_of_parameters(6,2): range_of_parameters(6,3);

Numel_paramSet= numel(A)*numel(B)*numel(C)*numel(D)*numel(E)*numel(F);
repetition_inx= 15;
objfun_param_List= nan(Numel_paramSet*repetition_inx, 6);

ctr=0;
for AA=A
    
    for BB=B
        
        for CC=C
            
            for DD=D
                
                for EE=E
                    
                    for FF=F
                        
                        for G=1:repetition_inx
                            ctr=ctr+1;
                            objfun_param_List(ctr,:)= [AA BB CC DD EE FF ];
                        end
                    end
                end
            end
        end
    end
end
%% simulation variables
t=1;
dt=1;
Tau=1;
coef= dt/Tau;
%% load data

load([data_path filename])

nSubjects = length(unique(subjects));
for ss = 1:nSubjects % loop over subjects
    
    % get data
    data.RT
    data.CP
    
    % preparing for simulations - initializing variables
    [noise_Go1, noise_Go2, noise_Stop] = deal(0);
    nonDecisionDelay = round(0.050*1000);
	nonDecisionDelay(2) = nonDecisionDelay(1);
    
    Sim_rStop = [];
    
    Sim_List_rGo1 = [zeros(Sim_nTrials, nonDecisionDelay(1)) ones(Sim_nTrials,  max_it-nonDecisionDelay(1))];
    Sim_List_rGo2 = [zeros(Sim_nTrials, nonDecisionDelay(2)) ones(Sim_nTrials,  max_it-nonDecisionDelay(2))];

    Sim_List_rStop = [];
    
    activity_Go1 = repmat(rand(Sim_nTrials,1), 1, max_it);
    activity_Go2 = repmat(rand(Sim_nTrials,1), 1, max_it);

    activity_Stop = zeros(Sim_nTrials, max_it);
    
    % list of simulated RTs
    Reaction_time = NaN(Sim_nTrials, 1);
    idx_overThreshold1 = NaN(Sim_nTrials, max_it);
    idx_overThreshold2 = NaN(Sim_nTrials, max_it);
    % list of Winners of the race of each trial
    Winner = NaN(Sim_nTrials,1);
    
    %% simulation
    objfun_param_results= nan(Numel_paramSet*repetition_inx, 14);
    
    for DA = 1: size(objfun_param_List,1)
        objfun_param = objfun_param_List(DA,:);
        
        
        Sim_rGo= [normrnd(objfun_param(find(strcmp('V1', freeParams))), 0.7e-3, Sim_nTrials, 1), ...
            normrnd(objfun_param(find(strcmp('V2', freeParams))), 0.7e-3, Sim_nTrials, 1)];
        
        Sim_List_rGo1= Sim_List_rGo1 .* repmat(Sim_rGo(:,1), 1, max_it);
        Sim_List_rGo2= Sim_List_rGo2 .* repmat(Sim_rGo(:,2), 1, max_it);
        
        % setting starting point
        activity_Go1= 0.1 * activity_Go1 + (objfun_param(find(strcmp('Sp1', freeParams)))-0.1);
        activity_Go1(activity_Go1<0)= -activity_Go1(activity_Go1<0);
        activity_Go2= 0.1 * activity_Go2 + (objfun_param(find(strcmp('Sp2', freeParams)))-0.1);
        activity_Go2(activity_Go2<0)= -activity_Go2(activity_Go2<0);
        
        %% calculate activities per step for trials (discretized stochastic differential equation from paper Ramakrishnan 2012)
        for t = min(nonDecisionDelay):max_it

            activity_Go1(:,t) = activity_Go1(:,t-1) + coef*(Sim_List_rGo1(:,t) ...
                - leakageCoef * activity_Go1(:,t-1) ...
                - (objfun_param(find(strcmp('B21', freeParams)))) * activity_Go2(:,t-1) ); %+ noise_Go1(:,t);
            
            activity_Go2(:,t)= activity_Go2(:,t-1)+ coef*(Sim_List_rGo2(:,t) ...
                - leakageCoef * activity_Go2(:,t-1) ...
                - (objfun_param(find(strcmp('B12', freeParams)))) * activity_Go1(:,t-1) ); %+ noise_Go2(:,t);
            
            activity_Go1(activity_Go1(:,t)<0,t) = 0; % reset the activity to zero if it's negative value!!!
            activity_Go2(activity_Go2(:,t)<0,t) = 0; % reset the activity to zero if it's negative value!!!
            idx_overThreshold1(activity_Go1(:,t)>=1,t) = t; % mark when the activity goes beyond the threshold
            idx_overThreshold2(activity_Go2(:,t)>=1,t) = t; % mark when the activity goes beyond the threshold
            
            %  end if any accumulator exceeds the threshold
            if sum(all(isnan([min(idx_overThreshold1,[],2) min(idx_overThreshold2,[],2)])==1,2))==0
                break
            end
        end
        
        %% kstest
        [Reaction_time,Winner] = min([min(idx_overThreshold1,[],2)  min(idx_overThreshold2,[],2)], [],2);
        if ~isempty(Obs_Data(1,1).RT) && ~isempty(Reaction_time(Winner==1,1)) && ~all(isnan(Reaction_time(Winner==1,1)))
            [~, Mdl_KS2_P(1), Mdl_KS2_stat(1)] = kstest2(Obs_Data(1,1).RT, Reaction_time(Winner==1,1));
            
        else
            [~, Mdl_KS2_P(1), Mdl_KS2_stat(1)] = deal(nan);
        end
        
        if ~isempty(Obs_Data(1,2).RT) && ~isempty(Reaction_time(Winner==2,1)) && ~all(isnan(Reaction_time(Winner==2,1)))
            [~, Mdl_KS2_P(2), Mdl_KS2_stat(2)] = kstest2(Obs_Data(1,2).RT, Reaction_time(Winner==2,1));
        else
            [~, Mdl_KS2_P(2), Mdl_KS2_stat(2)] = deal(nan);
        end
        
        Mdl_KS2_H = nan(1,2);
        
        sim_choicepref = [sum(Winner==1) sum(Winner==2)]/(sum(Winner==1)+sum(Winner==2));
        sim_choicepref_op = abs(sim_choicepref(1)-Obs_Data(1).SR_ovrl)*100;
        
        %% G2 test... QMLE # CHECK with DATA
        % prepare quantile probabilities based on the proportion of correct/incorrect or left/right choice
        obs_ovrl_SR_diff = abs(Obs_Data(1).SR_ovrl-Obs_Data(2).SR_ovrl)*100;
        switch true
            case obs_ovrl_SR_diff<30;
                qntl_probs = {[0.1 0.3 0.5 0.7 0.9 ]; [0.1 0.3 0.5 0.7 0.9  ]}; % for left; for right
                
            case obs_ovrl_SR_diff>=30 & obs_ovrl_SR_diff<=70;
                qntl_probs{[Obs_Data(1).SR_ovrl Obs_Data(2).SR_ovrl]<0.5, 1} = [0.3 0.7 ];
                qntl_probs{[Obs_Data(1).SR_ovrl Obs_Data(2).SR_ovrl]>0.5, 1} = [0.1 0.3 0.5 0.7 0.9 ];
                
            case obs_ovrl_SR_diff>70;
                qntl_probs{[Obs_Data(1).SR_ovrl Obs_Data(2).SR_ovrl]<0.5, 1} = [0.5 ];
                qntl_probs{[Obs_Data(1).SR_ovrl Obs_Data(2).SR_ovrl]>0.5, 1} = [0.1 0.3 0.5 0.7 0.9 ];
        end
        
        % quantile calculations
        [Sim_qntl_output] = DA_count_quantiles({Reaction_time(Winner==1,1); Reaction_time(Winner==2,1)}, qntl_probs); % simulated RT left right
        Nt_sim = sum([Sim_qntl_output.Nt]);
        Sim_qntl_prob_L = 100*Sim_qntl_output(1,1).QuantileBinCounts ./Nt_sim;
        Sim_qntl_prob_R = 100*Sim_qntl_output(1,2).QuantileBinCounts ./Nt_sim;
        
        
        for jl = 1:2 % 1 left   2 right # CHECK with DATA
            tmpqntl = [0 Sim_qntl_output(1,jl).dataQuantiles];
            Obs_qntl_output(1,jl).dataQuantiles= Sim_qntl_output(1,jl).dataQuantiles;
            for jk = 1:numel(tmpqntl)
                if jk==numel(tmpqntl)
                    Obs_qntl_output(1,jl).quantiled_data{jk} = Obs_Data(1,jl).RT(Obs_Data(1,jl).RT>=tmpqntl(jk));
                else
                    Obs_qntl_output(1,jl).quantiled_data{jk} = Obs_Data(1,jl).RT(Obs_Data(1,jl).RT>=tmpqntl(jk) & Obs_Data(1,jl).RT<tmpqntl(jk+1));
                end
                Obs_qntl_output(1,jl).QuantileBinCounts(jk) = numel(Obs_qntl_output(1,jl).quantiled_data{jk});
            end
            Obs_qntl_output(1,jl).Nt = sum(Obs_qntl_output(1,jl).QuantileBinCounts);
        end
        
        Nt_obs = numel(Obs_Data(1,1).RT)+ numel(Obs_Data(1,2).RT);
        Obs_qntl_prob_L = 100*Obs_qntl_output(1,1).QuantileBinCounts ./Nt_obs;
        Obs_qntl_prob_R = 100*Obs_qntl_output(1,2).QuantileBinCounts ./Nt_obs;
           
        
        % G2 test
        G_left  = 2 * nansum(Obs_qntl_prob_L .* abs(log(Obs_qntl_prob_L ./ Sim_qntl_prob_L)));
        G_right = 2 * nansum(Obs_qntl_prob_R .* abs(log(Obs_qntl_prob_R ./ Sim_qntl_prob_R)));
        
        % conditions # CHECK with DATA
        if sum(Sim_qntl_prob_L)==0; G_left  = Inf; end
        if sum(Sim_qntl_prob_R)==0; G_right = Inf; end
        
        if sum(Sim_qntl_output(1,1).QuantileBinCounts)<=5; G_left  = Inf; end
        if sum(Sim_qntl_output(1,2).QuantileBinCounts)<=5; G_right = Inf; end

        % CHECK with DATA
        switch  Sim_Task_cond
            case 'LL'
                G_2_test = [G_left G_right];
                ObjFun_op = G_left ; % objective function output is the sum of G_2_test for each side
                
            case 'RR'
                G_2_test = [ G_left G_right];
                ObjFun_op =  G_right; % objective function output is the sum of G_2_test for each side
                
            otherwise
                G_2_test = [G_left G_right];
                ObjFun_op = G_left + G_right; % objective function output is the sum of G_2_test for each side
        end
        objfun_param_results(DA,:) = [objfun_param, nan, nan, ObjFun_op,  G_2_test, Mdl_KS2_stat, sim_choicepref_op];
    end
    
    disp('10000 times repetition')
    
    G2fit_rep_param = reshape(ObjFun_op, repetition_inx, numel(ObjFun_op)/repetition_inx);
    median_G2fit_replaced = reshape(repmat(nanmedian(G2fit_rep_param, 1), repetition_inx,1), numel(ObjFun_op),1);
    [sort_ObjFun_op, sort_ObjFun_opInx] = sort(median_G2fit_replaced);
    
    tmp_objfun_param_sorted = objfun_param_results(sort_ObjFun_opInx(1:repetition_inx*20), 1:6);
    
    tmp_objfun_param = tmp_objfun_param_sorted(1:repetition_inx:repetition_inx*20, 1:6);
    tmp_repeated_objfun_param = [];
    med_repetition_10000_G2sum = [];
    N_itr_Loop = 10000;
    
    for jj = 1:20
        objfun_param = tmp_objfun_param(jj,:);
        for DA = 1: N_itr_Loop 
        
        
        Sim_rGo = [normrnd(objfun_param(find(strcmp('V1', freeParams))), 0.7e-3, Sim_nTrials, 1), ...
            normrnd(objfun_param(find(strcmp('V2', freeParams))), 0.7e-3, Sim_nTrials, 1)];
        
        Sim_List_rGo1 = Sim_List_rGo1 .* repmat(Sim_rGo(:,1), 1, max_it);
        Sim_List_rGo2 = Sim_List_rGo2 .* repmat(Sim_rGo(:,2), 1, max_it);
        
        % setting starting point
        activity_Go1 = 0.1 * activity_Go1 + (objfun_param(find(strcmp('Sp1', freeParams)))-0.1);
        activity_Go1(activity_Go1<0) = -activity_Go1(activity_Go1<0);
        activity_Go2 = 0.1 * activity_Go2 + (objfun_param(find(strcmp('Sp2', freeParams)))-0.1);
        activity_Go2(activity_Go2<0) = -activity_Go2(activity_Go2<0);
        
        %% calculate activities per step for trials (discretized stochastic differential equation from paper Ramakrishnan 2012)
        for t = min(nonDecisionDelay):max_it

            activity_Go1(:,t) = activity_Go1(:,t-1) + coef*(Sim_List_rGo1(:,t) ...
                - leakageCoef * activity_Go1(:,t-1) ...
                - (objfun_param(find(strcmp('B21', freeParams)))) * activity_Go2(:,t-1) ); %+ noise_Go1(:,t);
            
            activity_Go2(:,t)= activity_Go2(:,t-1)+ coef*(Sim_List_rGo2(:,t) ...
                - leakageCoef * activity_Go2(:,t-1) ...
                - (objfun_param(find(strcmp('B12', freeParams)))) * activity_Go1(:,t-1) ); %+ noise_Go2(:,t);
            
            activity_Go1(activity_Go1(:,t)<0,t) = 0; % reset the activity to zero if it's negative value!!!
            activity_Go2(activity_Go2(:,t)<0,t) = 0; % reset the activity to zero if it's negative value!!!
            idx_overThreshold1(activity_Go1(:,t)>=1,t) = t; % mark when the activity goes beyond the threshold
            idx_overThreshold2(activity_Go2(:,t)>=1,t) = t; % mark when the activity goes beyond the threshold
            
            %  end if any accumulator exceeds the threshold
            if sum(all(isnan([min(idx_overThreshold1,[],2) min(idx_overThreshold2,[],2)])==1,2))==0
                break
            end
        end
        
        %% kstest
        [Reaction_time,Winner] = min([min(idx_overThreshold1,[],2)  min(idx_overThreshold2,[],2)], [],2);
        if ~isempty(Obs_Data(1,1).RT) && ~isempty(Reaction_time(Winner==1,1)) && ~all(isnan(Reaction_time(Winner==1,1)))
            [~, Mdl_KS2_P(1), Mdl_KS2_stat(1)] = kstest2(Obs_Data(1,1).RT, Reaction_time(Winner==1,1));
            
        else
            [~, Mdl_KS2_P(1), Mdl_KS2_stat(1)] = deal(nan);
        end
        
        if ~isempty(Obs_Data(1,2).RT) && ~isempty(Reaction_time(Winner==2,1)) && ~all(isnan(Reaction_time(Winner==2,1)))
            [~, Mdl_KS2_P(2), Mdl_KS2_stat(2)] = kstest2(Obs_Data(1,2).RT, Reaction_time(Winner==2,1));
        else
            [~, Mdl_KS2_P(2), Mdl_KS2_stat(2)] = deal(nan);
        end
        
        Mdl_KS2_H = nan(1,2);
        
        sim_choicepref = [sum(Winner==1) sum(Winner==2)]/(sum(Winner==1)+sum(Winner==2));
        sim_choicepref_op = abs(sim_choicepref(1)-Obs_Data(1).SR_ovrl)*100;
        
        %% G2 test... QMLE # CHECK with DATA
        % prepare quantile probabilities based on the proportion of correct/incorrect or left/right choice
        obs_ovrl_SR_diff = abs(Obs_Data(1).SR_ovrl-Obs_Data(2).SR_ovrl)*100;
        switch true
            case obs_ovrl_SR_diff<30;
                qntl_probs = {[0.1 0.3 0.5 0.7 0.9 ]; [0.1 0.3 0.5 0.7 0.9  ]}; % for left; for right
                
            case obs_ovrl_SR_diff>=30 & obs_ovrl_SR_diff<=70;
                qntl_probs{[Obs_Data(1).SR_ovrl Obs_Data(2).SR_ovrl]<0.5, 1} = [0.3 0.7 ];
                qntl_probs{[Obs_Data(1).SR_ovrl Obs_Data(2).SR_ovrl]>0.5, 1} = [0.1 0.3 0.5 0.7 0.9 ];
                
            case obs_ovrl_SR_diff>70;
                qntl_probs{[Obs_Data(1).SR_ovrl Obs_Data(2).SR_ovrl]<0.5, 1} = [0.5 ];
                qntl_probs{[Obs_Data(1).SR_ovrl Obs_Data(2).SR_ovrl]>0.5, 1} = [0.1 0.3 0.5 0.7 0.9 ];
        end
        
        % quantile calculations
        [Sim_qntl_output] = DA_count_quantiles({Reaction_time(Winner==1,1); Reaction_time(Winner==2,1)}, qntl_probs); % simulated RT left right
        Nt_sim = sum([Sim_qntl_output.Nt]);
        Sim_qntl_prob_L = 100*Sim_qntl_output(1,1).QuantileBinCounts ./Nt_sim;
        Sim_qntl_prob_R = 100*Sim_qntl_output(1,2).QuantileBinCounts ./Nt_sim;
        
        
        for jl = 1:2 % 1 left   2 right # CHECK with DATA
            tmpqntl = [0 Sim_qntl_output(1,jl).dataQuantiles];
            Obs_qntl_output(1,jl).dataQuantiles= Sim_qntl_output(1,jl).dataQuantiles;
            for jk = 1:numel(tmpqntl)
                if jk==numel(tmpqntl)
                    Obs_qntl_output(1,jl).quantiled_data{jk} = Obs_Data(1,jl).RT(Obs_Data(1,jl).RT>=tmpqntl(jk));
                else
                    Obs_qntl_output(1,jl).quantiled_data{jk} = Obs_Data(1,jl).RT(Obs_Data(1,jl).RT>=tmpqntl(jk) & Obs_Data(1,jl).RT<tmpqntl(jk+1));
                end
                Obs_qntl_output(1,jl).QuantileBinCounts(jk) = numel(Obs_qntl_output(1,jl).quantiled_data{jk});
            end
            Obs_qntl_output(1,jl).Nt = sum(Obs_qntl_output(1,jl).QuantileBinCounts);
        end
        
        Nt_obs = numel(Obs_Data(1,1).RT)+ numel(Obs_Data(1,2).RT);
        Obs_qntl_prob_L = 100*Obs_qntl_output(1,1).QuantileBinCounts ./Nt_obs;
        Obs_qntl_prob_R = 100*Obs_qntl_output(1,2).QuantileBinCounts ./Nt_obs;
           
        
        % G2 test
        G_left  = 2 * nansum(Obs_qntl_prob_L .* abs(log(Obs_qntl_prob_L ./ Sim_qntl_prob_L)));
        G_right = 2 * nansum(Obs_qntl_prob_R .* abs(log(Obs_qntl_prob_R ./ Sim_qntl_prob_R)));
        
        % conditions # CHECK with DATA
        if sum(Sim_qntl_prob_L)==0; G_left  = Inf; end
        if sum(Sim_qntl_prob_R)==0; G_right = Inf; end
        
        if sum(Sim_qntl_output(1,1).QuantileBinCounts)<=5; G_left  = Inf; end
        if sum(Sim_qntl_output(1,2).QuantileBinCounts)<=5; G_right = Inf; end

        % CHECK with DATA
        switch  Sim_Task_cond
            case 'LL'
                G_2_test = [G_left G_right];
                ObjFun_op = G_left ; % objective function output is the sum of G_2_test for each side
                
            case 'RR'
                G_2_test = [ G_left G_right];
                ObjFun_op =  G_right; % objective function output is the sum of G_2_test for each side
                
            otherwise
                G_2_test = [G_left G_right];
                ObjFun_op = G_left + G_right; % objective function output is the sum of G_2_test for each side
        end
        tmp_repeated_objfun_param(DA,:)= [objfun_param, nan, nan, ObjFun_op,  G_2_test, Mdl_KS2_stat, sim_choicepref_op];
        end
        repeated_objfun_param{1,jj} = tmp_repeated_objfun_param;
        med_repetition_10000_G2sum(jj,1) = median(tmp_repeated_objfun_param(:,9));
    end
    [min_ObjFun, min_ObjFunInx]= min(med_repetition_10000_G2sum);
    min_objfun_param_final= repeated_objfun_param{min_ObjFunInx}(1:6);
    min_ObjFun_final= min_ObjFun(1);
    
    SIM_OUTPUT.Obs_Data= Obs_Data;
    SIM_OUTPUT.range_of_parameters=range_of_parameters;
    SIM_OUTPUT.param_results_fulllist= objfun_param_results(sort_ObjFun_opInx(1:1500), :);
    SIM_OUTPUT.repeated_objfun_param= repeated_objfun_param;
    SIM_OUTPUT.parameters= min_objfun_param_final;
    SIM_OUTPUT.Opt_init_val= Opt_init_val;
    SIM_OUTPUT.FMINOUTPUT.fmin_x=min_objfun_param_final;
    SIM_OUTPUT.FMINOUTPUT.fmin_feval=min_ObjFun_final;
    SIM_OUTPUT.FMINOUTPUT.fminEXITFLAG=NaN;
    SIM_OUTPUT.FMINOUTPUT.fminoutput=NaN;
    SIM_OUTPUT.FMINOUTPUT.maxfuneval= NaN;
    SIM_OUTPUT.FMINOUTPUT.maxit= NaN;
    SIM_OUTPUT.FMINOUTPUT.tolfun= NaN;
    SIM_OUTPUT.FMINOUTPUT.tolx= NaN;
    SIM_OUTPUT.FMINOUTPUT.disp= NaN;
    SIM_OUTPUT.FMINOUTPUT.N_itr_Loop= [repetition_inx N_itr_Loop];
    
    SIM_OUTPUT.Reaction_time= [];
    SIM_OUTPUT.Winner= [];
    SIM_OUTPUT.sim_choicepref= [];
    SIM_OUTPUT.ObjFun_op= [];
    SIM_OUTPUT.KStest= [];
    SIM_OUTPUT.G_2_test= [];
    
    SIM_OUTPUT= MDa_addfieldsstruct2struct(SIM_OUTPUT, tmpStruct);
    Simulation_output(Subject,1)= SIM_OUTPUT;
    
    clc

end
eval(['Simulation_output_' Sim_Task_cond '_' Sim_model_type '=' 'Simulation_output']);
eval([' save Simulation_output_' Sim_Task_cond '_' Sim_model_type ' Simulation_output_' Sim_Task_cond '_' Sim_model_type]);

toc/3600;

