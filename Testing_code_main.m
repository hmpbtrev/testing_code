%% Modeling the Effectiveness of RT-PCR, RT-LAMP, and Antigen Testing Strategies for COVID-19 Control
%  
% Description:  
% This MATLAB script implements a mathematical model to evaluate the  
% effectiveness of RT-PCR, RT-LAMP, and antigen testing strategies for  
% controlling COVID-19 transmission.  
%  
% Usage:  
% - Run the main script to perform simulations and generate results.  
%  
% *Corresponding Author: Charin Modchang (charin.mod@mahidol.edu)  

%% Clear work space
clear
clc
close all

%% Batch
Batch = 100; % number of batches

for iBatch=1:Batch
    
    %% Clear work space
    clearvars -except iBatch
    close all;

    %% Random seed
    rng('shuffle')

    %% Interation
    Nrepeat = 500; % # of simulation in a batch
    %% CC on/off
    cc_par = 1; % 1 = on / 0 = off
    %% Disease parameter
    %  Incubation period parameters
    TG_a = 3.6428;     
    TG_b = 1.5907;

    % Proportion of cases that are eventually symptomatic
    k_sym_S = 0.573;  
    k_sym_V = 0.431;
    % k_sym_S = 0.0;  
    % k_sym_V = 0.0;

    % Reduction of infectiousness of asymptomatic infection
    k_reduce = 0.58;

    % IFPChoice = [IFPS IFPQ];
    IChoice = [5 6];        % I => R or I => X

    %% Distribution of " the number of secondary infectious caused by each case (Z)"
    Dist = 3;                  	%  1) Poisson
                                %  2) Geometric                           
                                %  3) Negative binomial
    R00 = 5.08; k=0.08; 
%     R00 = 2.2; k = 0.1;  
%     ZSSE = 9; % * k Correlate to R0. If R0 change, k change >> CHECK 
    % R00 = 10;

    %% Vaccine-induced protection
    %  Vaccine efficiency
    eS = 0.79;              % efficinecy against infection
    eI = 1-0.7545;          % efficiency against transmission
    %  Vaccine coverage
    vc = 0/100;

    %% Infectiousness profile parameters
    load('VL_profile.mat')
    load('VL_lin_profile.mat')

    Inf_model = 3;      % 1) Linear 2)Power law 3) Saturation
    pow_law_para1 = 2.4*10^-5;
    pow_law_para2 = 0.53;
    sat_mod_para1 = 8.9*10^6;
    sat_mod_para2 = 0.51;
    sat_mod_para3 = 0.2;
    if Inf_model == 2
        inf_profile_UV  = [UV_all_lin(:,1) pow_law_para1.*UV_all_lin(:,2).^pow_law_para2];
        inf_profile_V   = [V_all_lin(:,1) pow_law_para1.*V_all_lin(:,2).^pow_law_para2];
    elseif Inf_model == 3
        inf_profile_UV = [UV_all_lin(:,1) (sat_mod_para3*UV_all_lin(:,2).^sat_mod_para2)./((UV_all_lin(:,2).^sat_mod_para2)+(sat_mod_para1^sat_mod_para2))];
        inf_profile_V = [V_all_lin(:,1) (sat_mod_para3*V_all_lin(:,2).^sat_mod_para2)./((V_all_lin(:,2).^sat_mod_para2)+(sat_mod_para1^sat_mod_para2))];
    else
    end
    inf_profile_UV  = inf_profile_UV(1:3310,:);
    inf_profile_V   = inf_profile_V(1:3310,:);

    %% percent of infection from close contact
    infToCC = 0.15;

    %% Isolation
    % Time delay to isolation since detection
    tdelay = 1;

    %% Time window (days)
    t0 = 0;
    tf = 301; % simulation time
    dt = 1e-2;
    dtsave = 1;    
    Nstep = round((tf-t0)/dt);
    Nsave = (round((tf-t0)/dtsave));

    %% Testing parameters 
    % symptom testing
    symptom_test_pop            = 1; % [0-1], proportion of symptom-based testing 
    % routine testing
    routine_test_pc             = 0; % [0-1] to be varied, proportion of daily screened
    routine_test_freq           = 1; % in days
    % random testing
    random_test_capa            = 0; % x test per day
    % close contact testing
    close_contact_test_pc       = 0; % percentage of close contacts to be tested                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ; % [0-1]
    close_contact_test_freq     = 1; % in days
    close_contact_test_period   = 7; % period of cc testing
    % isolation
    isolation_test_freq         = 4; % frequency of testing in isolation
    %
    LOD_pcr             = log10(10^3);     % Limit of detection (log10(Vmin)), PCR
    LOD_lamp            = log10(3*10^3);   % LAMP
    LOD_atk             = log10(10^5);     % ATK
    %--- LOD ---
    LOD_symp            = LOD_lamp;     
    LOD_rou             = LOD_lamp;   
    LOD_ran             = LOD_atk;    
    LOD_cc              = LOD_lamp;   
    FN              = 0;     % False negative
    FP              = 0;     % False positive
    % Positive result since infection
    %PCR
    UV_all_pos_pcr       = UV_all((UV_all(:,2) > LOD_pcr),1);
    V_all_pos_pcr        = V_all((V_all(:,2) > LOD_pcr),1);
    %LAMP
    UV_all_pos_lamp      = UV_all((UV_all(:,2) > LOD_lamp),1);
    V_all_pos_lamp       = V_all((V_all(:,2) > LOD_lamp),1);
    %ATK
    UV_all_pos_atk       = UV_all((UV_all(:,2) > LOD_atk),1);
    V_all_pos_atk        = V_all((V_all(:,2) > LOD_atk),1);
    %--
    UV_all_pos_symp     = UV_all_pos_pcr;
    V_all_pos_symp      = V_all_pos_pcr;
    UV_all_pos_rou      = UV_all_pos_pcr;
    V_all_pos_rou       = V_all_pos_pcr;
    UV_all_pos_ran      = UV_all_pos_atk;
    V_all_pos_ran       = V_all_pos_atk;
    UV_all_pos_cc       = UV_all_pos_pcr;
    V_all_pos_cc        = V_all_pos_pcr;
    % turnaround time
    turnAround_time_pcr  = 1;
    turnAround_time_lamp = 0;
    turnAround_time_atk  = 0;
    % turnaround time
    turnAround_time_symp    = turnAround_time_pcr;
    turnAround_time_rou     = turnAround_time_pcr;
    turnAround_time_ran     = turnAround_time_atk;
    turnAround_time_cc      = turnAround_time_lamp;
    % import case (per week)
    num_import_case = 0;

    %% Storage Matrix
    tsave = NaN(Nsave,Nrepeat);
    Ssave = NaN(Nsave,Nrepeat);
    LSsave = NaN(Nsave,Nrepeat);
    ISSsave = NaN(Nsave,Nrepeat);
    IASsave = NaN(Nsave,Nrepeat);

    Vsave = NaN(Nsave,Nrepeat);
    LVsave = NaN(Nsave,Nrepeat);
    ISVsave = NaN(Nsave,Nrepeat);
    IAVsave = NaN(Nsave,Nrepeat);

    Iallsave = NaN(Nsave,Nrepeat);

    RSsave = NaN(Nsave,Nrepeat);
    RVsave = NaN(Nsave,Nrepeat);
    Rsave = NaN(Nsave,Nrepeat);
    Rallsave = NaN(Nsave,Nrepeat);
    QSSsave = NaN(Nsave,Nrepeat);
    QTSsave = NaN(Nsave,Nrepeat);
    QSVsave = NaN(Nsave,Nrepeat);
    QTVsave = NaN(Nsave,Nrepeat);
    Qallsave = NaN(Nsave,Nrepeat);
    % Qsave = NaN(Nsave,Nrepeat);
    PLSsave = NaN(Nsave,Nrepeat);
    PLVsave = NaN(Nsave,Nrepeat);

    CumISSsave = NaN(Nsave,Nrepeat);
    CumIASsave = NaN(Nsave,Nrepeat);
    CumISVsave = NaN(Nsave,Nrepeat);
    CumIAVsave = NaN(Nsave,Nrepeat);
    CumISsave = NaN(Nsave,Nrepeat);
    CumIAsave = NaN(Nsave,Nrepeat);
    CumIallsave = NaN(Nsave,Nrepeat);

    CumQSSsave = NaN(Nsave,Nrepeat);
    CumQTSsave = NaN(Nsave,Nrepeat);
    CumQSVsave = NaN(Nsave,Nrepeat);
    CumQTVsave = NaN(Nsave,Nrepeat);
    CumXsave = NaN(Nsave,Nrepeat);
    CumRSsave = NaN(Nsave,Nrepeat);
    CumRVsave = NaN(Nsave,Nrepeat);
    CumRsave = NaN(Nsave,Nrepeat);

    LS_count_save   = NaN(Nstep,Nrepeat);
    LV_count_save   = NaN(Nstep,Nrepeat);
    Lall_count_save = NaN(Nstep,Nrepeat);

    Rt_save = zeros(Nsave,Nrepeat);
    Rt2_save = zeros(Nstep,Nrepeat);
    meanDetectionTime = zeros(Nrepeat,1);
    
    CCtest_matrix_save = nan(Nrepeat,3);
    num_outbreak = 0;
    
    Test_count = nan(Nrepeat,1);
    
    for irepeat = 1:Nrepeat

        %  Initial Population 
        N0 = 1e4; 

        NLS0 = 1;         
        NISS0 = 0;       
        NIAS0 = 0; 
        NRS0 = 0.0*N0;
        NQSS0 = 0;
        NQTS0 = 0;
        NPLS0 = 0; 

        NLV0 = 0;        
        NISV0 = 0;      
        NIAV0 = 0; 
        NRV0 = 0;
        NQSV0 = 0;
        NQTV0 = 0;
        NPLV0 = 0; 

        % N = #Vaccinated population + Unvaccinated population
        % V 
        NV0 = (round(vc*N0)) - NLV0 - NISV0 - NIAV0 - NRV0 - NQSV0 - NQTV0 - NPLS0; 
        % UV
        NS0 = N0 - NV0 - NLS0 - NISS0 - NIAS0 - NRS0 - NQSS0 - NQTS0 - NPLS0;

        % Count number of cases
        CumISS = NISS0;      
        CumIAS = NIAS0;
        CumISV = NISV0;      
        CumIAV = NIAV0;
        CumIS = CumISS+CumISV;      
        CumIA = CumIAS+CumIAV;

        CumQTS = NQTS0; 
        CumQSS = NQSS0; 
        CumQSV = NQSV0;
        CumQTV = NQTV0;
        CumQ = CumQTS + CumQSS + CumQTV + CumQSV; 

        CumRV = NRV0;
        CumRS = NRS0;
        CumR = CumRV + CumRS; 

        %% State matrix
        % State =       1)ID            2)CurrentStatus     3)NextStatus   
        %               4)IS/IA         5)Q/R               6)Speader'sID  
        %               7)Z             8)Secondary infection 

        % Infection state: 1)S   2)LS   3)ISS   4)IAS   5)RS   6)QSS   7)QTS 
        %                  10)PLS (Pre-infected state)
        %                  11)V  12)LV  13)ISV  14)IAV  15)RV  16)QSV  17)QTV    
        %                  20)PLV 

        % Speader's ID:     NaN for not yet infected individual (S) 
        %                   0 for gettiing infection fron index case

        State = NaN(N0,8); 
        State(:,1) = 1:N0;
        State(:,2) = 1;                 % All initially assiged as S

        %% Alarm matrix
        %  Alarm = [ 1)ID       2)t in current state   3)t in S         
        %            4)t in E   5)t in I]

        Alarm = NaN(N0,5);
        Alarm(:,1) = 1:N0;
        Alarm(:,2) = 0;

        % initial time
        t = t0;
        day = 0;

        %% Initial population
        int = randperm(N0,N0-NS0);
        idxIniLS = int(1:NLS0);
        idxIniLV = int(NLS0+1:NLS0+NLV0);

        idxIniISS = int(NLS0+NLV0+1:NLS0+NLV0+NISS0);
        idxIniISV = int(NLS0+NLV0+NISS0+1:NLS0+NLV0+NISS0+NISV0);

        idxIniIAS = int(NLS0+NLV0+NISS0+NISV0+1:NLS0+NLV0+NISS0+NISV0+NIAS0);
        idxIniIAV = int(NLS0+NLV0+NISS0+NISV0+NIAS0+1:NLS0+NLV0+NISS0+NISV0+...
            NIAS0+NIAV0);

        idxIniRS = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+1:NLS0+NLV0+NISS0+...
            NISV0+NIAS0+NIAV0+NRS0);
        idxIniRV = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+1:NLS0+NLV0+...
            NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0);

        idxIniQSS = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+1:NLS0+...
            NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0);
        idxIniQSV = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+...
            1:NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+NQSV0);

        idxIniQTS = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+...
            NQSV0+1:NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+NQSV0+...
            NQTS0);
        idxIniQTV = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+...
            NQSV0+NQTS0+1:NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+...
            NQSV0+NQTS0+NQTV0);

        idxIniPLS = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+...
            NQSV0+NQTS0+NQTV0+1:NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+...
            NQSS0+NQSV0+NQTS0+NQTV0+NPLS0);
        idxIniPLV = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+...
            NQSV0+NQTS0+NQTV0+NPLS0+1:NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+...
            NRS0+NRV0+NQSS0+NQSV0+NQTS0+NQTV0+NPLS0+NPLV0);

        idxIniV = int(NLS0+NLV0+NISS0+NISV0+NIAS0+NIAV0+NRS0+NRV0+NQSS0+...
            NQSV0+NQTS0+NQTV0+NPLS0+NPLV0+1:NLS0+NLV0+NISS0+NISV0+NIAS0+...
            NIAV0+NRS0+NRV0+NQSS0+NQSV0+NQTS0+NQTV0+NPLS0+NPLV0+NPLV0);

        % routine tester
        idxRouT = randperm(N0,round(routine_test_pc*N0))'; % id of routine tester

        % close contact testing
        TT_CCT                 = [];
        CloseContact_data      = [];
        Iso_data               = [];
        Testing_data           = [];
        idxCCT                 = [nan nan];
        DetectionTime_all      = inf(N0,1);
        S_Q_all                = cell(N0,2);
        Testing_time           = nan(N0,1);

        %-----close contact matrix-----
        if cc_par==1
        load('sfnetwork1e4link17_R0508.mat')

            for cc_ind=1:N0        
                cc_temp2 = randperm(length(CC_matrix{cc_ind}), round(close_contact_test_pc*length(CC_matrix{cc_ind}))); 
                CC_matrix_toTest{cc_ind} = CC_matrix{cc_ind}(cc_temp2);
                cc_temp2 = [];
            end
            CC_matrix_toTest=CC_matrix_toTest';
        end

        State(idxIniV,2) = 11;
        State(idxIniLS,2) = 2;      State(idxIniLV,2) = 12;
        State(idxIniISS,2) = 3;     State(idxIniISV,2) = 13;
        State(idxIniIAS,2) = 4;     State(idxIniIAV,2) = 14;
        State(idxIniRS,2) = 5;      State(idxIniRV,2) = 15;
        State(idxIniQSS,2) = 6;     State(idxIniQSV,2) = 16;
        State(idxIniQTS,2) = 7;     State(idxIniQTV,2) = 17;
        State(idxIniPLS,2) = 10;    State(idxIniPLV,2) = 20;

        % Spreader ID 
        % O = iniPL -1 = iniI -2 = iniQ -3 = iniR -4 = imported
        State(idxIniLS,6)  = 0;     State(idxIniLV,6) = 0;
        State(idxIniISS,6) = -1;    State(idxIniISV,6) = -1;
        State(idxIniIAS,6) = -1;    State(idxIniIAV,6) = -1;
        State(idxIniQTS,6) = -2;    State(idxIniQTV,6) = -2;
        State(idxIniQSS,6) = -2;    State(idxIniQSV,6) = -2;
        State(idxIniRS,6 ) = -3;    State(idxIniRV,6)  = -3;

        Rt = zeros(Nsave,2);
        Rt(:,1) = 1:tf;
        Rt2 = zeros(Nstep,2);
        Rt2(:,1) = [dt:dt:tf]';

        %% Next state
        % LS
        for i = 1:NLS0
            SorA = rand();
            RorQ = rand();

            if SorA <= k_sym_S     % ISS
                State(idxIniLS(i),4) = 3;
                R0 = R00;
            else                % IA
                State(idxIniLS(i),4) = 4;
                R0 = R00*k_reduce;
            end

            % Record next state 
            State(idxIniLS(i),3) = State(idxIniLS(i),4);  

            % Incubation period
            tc = 5.1; % 5.1 for both unvaccinated and vaccinated

            p   = (1+(R0)/k)^-1;   
            Z  = nbinrnd(k,p);
            if cc_par==1
                idxCC = CC_matrix{idxIniLS(i)};
                idxCC_toTest = CC_matrix_toTest{idxIniLS(i)};
                if State(idxIniLS(i),4) == 3 % IS
                    Z = Z_all(idxIniLS(i),1);             
                elseif State(idxIniLS(i),4) == 4 % IA
                    Z = Z_all(idxIniLS(i),2);   
                end
            end

            % detection time (time=0 since infection)
            PosRouTestTime  = inf;
            RanTestTime     = inf;  
            PosCCTestTime   = inf;
    %         IsoP            = 0;

            DetectionTime = inf;
            
            VL_since_iso            = [];
            num_test_till_first_neg = [];

            if State(idxIniLS(i),4) == 3 % IS
                if rand() < symptom_test_pop
                    DetectionTime = tc + turnAround_time_symp;  % latest detection time for IS is when symptom onset 
                
                    VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                    num_test_till_first_neg = min(find(VL_since_iso < LOD_symp))-1;  % number of test until first negative result
             
                end
            end 

            %---------routine testing--------
            if ismember(idxIniLS(i),idxRouT) == 1  
                PosRouTestTime = round(routine_test_freq*rand(1) + min(UV_all_pos_rou),2) + turnAround_time_rou;    
                DetectionTime = min([DetectionTime PosRouTestTime]);   
            
                VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_rou))-1;  % number of test until first negative result
            
            end
            %---------------------------------

            %----------random testing---------
            RanTT = rand(size(UV_all,1),1) < random_test_capa/N0*dt; 
            RanTT_dt = find(RanTT==1);
            for RanTT_ii=1:length(RanTT_dt)
                if UV_all(RanTT_dt(RanTT_ii),2) > LOD_ran
                    RanTestTime = UV_all(RanTT_dt(RanTT_ii),1) + turnAround_time_ran;
                    DetectionTime = min([DetectionTime RanTestTime]);
                    
                    VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                    num_test_till_first_neg = min(find(VL_since_iso < LOD_ran))-1;  % number of test untill first negative result
            
                end
                
            end
            %----------------------------------

            %----------CC testing---------
            if cc_par==1
                if (Testing_time(idxIniLS(i)) > 0) && (Testing_time(idxIniLS(i))~=inf)
                    if UV_all(UV_all(:,1)==Testing_time(idxIniLS(i)),2) > LOD_cc
                        LastCCTestTime = Testing_time(idxIniLS(i));
                        while UV_all(UV_all(:,1) == round(LastCCTestTime-close_contact_test_freq,2),2) > LOD_cc
                            LastCCTestTime = LastCCTestTime - close_contact_test_freq;
                        end
                        PosCCTestTime = LastCCTestTime + turnAround_time_cc;
                        DetectionTime = min([DetectionTime PosCCTestTime]);
                
                        VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                        num_test_till_first_neg = min(find(VL_since_iso < LOD_cc))-1;  % number of test untill first negative result
            
                    end
                end
            end
            %----------------------------------

            %-------------isolation------------
            IsoP = num_test_till_first_neg*isolation_test_freq; % Isolation period = number of test till first negative * test freq.
            date_out_Q = DetectionTime + IsoP;  % reference time since infection (=0 since infection) 
            %----------------------------------

            Iso_data = [Iso_data;idxIniLS(i) DetectionTime~=inf DetectionTime {date_out_Q} {IsoP} {num_test_till_first_neg}];

            Testing_data = [Testing_data;idxIniLS(i) ismember(idxIniLS(i),idxRouT) PosRouTestTime...
               ~isempty(find(RanTT==1)) RanTestTime ...
                PosCCTestTime DetectionTime];

            State(idxIniLS(i),7) = Z; 

            % Making infection

            if cc_par==1
                idxInfected = idxCC(randperm(length(idxCC),Z));       % random infectees    
            else
                idxInfected = randperm(N0,Z);
            end

            TS_ZQ = [];         TS_ZR = [];         % for ts of each infectee
            iNZQ_S = 0;         iNZR_S = 0;         % # successful infection
            idexS_Q = [];       idexS_R = [];       % index of infectee (successful infection)

            %--------------------------
            % If I => RS , No detection 
            for iZR = 1:Z
                TSrand = randpdf(inf_profile_UV(:,2),inf_profile_UV(:,1),[1,1]);
                if State(idxInfected(iZR),2) == 1   % S     % (infectee is Susceptible)
                    iNZR_S = iNZR_S+1;
                    idexS_R(iNZR_S) = idxInfected(iZR);
                    TS_ZR(iNZR_S) = TSrand;
                elseif State(idxInfected(iZR),2) == 11  % V
                    pv = rand();
                    if pv > eS
                        iNZR_S = iNZR_S+1; 
                        idexS_R(iNZR_S) = idxInfected(iZR);        
                        TS_ZR(iNZR_S) = TSrand;
                    end
                elseif State(idxInfected(iZR),2) == 10  % PLS
                    iNZR_S = iNZR_S+1; 
                    idexS_R(iNZR_S) = idxInfected(iZR);
                    TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]); 
                elseif State(idxInfected(iZR),2) == 20  % PLV
                    pv = rand();
                    if pv > eS
                        iNZR_S = iNZR_S+1; 
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]);
                    end
                end
            end

    %         %------------------------%
    %         % If I => Q (QSS/QTS) 
            for iZQ = 1:Z
                TSrand = randpdf(inf_profile_UV(:,2),inf_profile_UV(:,1),[1,1]);
                if State(idxInfected(iZQ),2) == 1   % S 
                    if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                        iNZQ_S = iNZQ_S+1;
                        idexS_Q(iNZQ_S) = idxInfected(iZQ);
                        TS_ZQ(iNZQ_S) = TSrand;
                    end

                elseif State(idxInfected(iZQ),2) == 11   % V 
                    pv = rand();
                    if pv > eS
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S+1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = TSrand;
                        end
                    end
                elseif State(idxInfected(iZQ),2) == 10   % PLS
                    if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                        iNZQ_S = iNZQ_S+1;
                        idexS_Q(iNZQ_S) = idxInfected(iZQ);
                        TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                    end
                elseif State(idxInfected(iZQ),2) == 20  % PLV
                    pv = rand();
                    if pv > eS
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S+1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                        end
                    end   
                end
            end
            %----------------------------------%

            DetectionTime_all(idxIniLS(i))  = DetectionTime;

            %---------------------------------%
            if cc_par==1
                idxCC_TestedandInfected = idxCC_toTest(ismember(idxCC_toTest,idexS_Q));
                for ii1=1:length(idxCC_TestedandInfected) 
                    Testing_time_temp = inf;
                    if DetectionTime ~= inf
                        Testing_time_temp = round(DetectionTime-TS_ZQ(ii1)+close_contact_test_period-mod(close_contact_test_period,close_contact_test_freq),2);
                    end
                    Testing_time(idxCC_TestedandInfected(ii1)) = min([Testing_time(idxCC_TestedandInfected(ii1)) Testing_time_temp]);
                end
                CloseContact_data = [CloseContact_data;idxIniLS(i) {idxCC} {idexS_Q} {TS_ZQ} {idxCC_toTest} {idxCC_TestedandInfected}];
            end
            %---------------------------------%

            % Secondary infection    

            if cc_par==1
                if DetectionTime == inf 
                    State(idxIniLS(i),5) = 5;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniLS(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) ||...
                (isnan(Testing_time(idxIniLS(i)))==0 && DetectionTime==PosCCTestTime+turnAround_time_cc)  % QT
                    State(idxIniLS(i),5) = 7;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniLS(i),5) = 6;
                    Z2 = length(TS_ZQ); 
                end
            else
                if DetectionTime == inf        % not quarantined
                State(idxIniLS(i),5) = 5;
                Z2 = length(TS_ZR); 
                elseif (ismember(idxIniLS(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) % QT
                    State(idxIniLS(i),5) = 7;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniLS(i),5) = 6;
                    Z2 = length(TS_ZQ); 
                end
            end

            State(idxIniLS(i),8) = Z2;

            % Update next state

            if (State(idxIniLS(i),5) == 5) % (UV) L >> R
                idv = find(State(idexS_R,2)==11);   
                State(idexS_R,2) = 10;          % State(:,2) = current state
                State(idexS_R,3) = 2;           % State(:,3) = next state
                State(idexS_R(idv),3) = 12;
                State(idexS_R,6) = idxIniLS(i); % State(:,6) = Spreader
                Alarm(idexS_R,2) = 0;           % Alarm(:,2) = t in current state
                Alarm(idexS_R,3) = TS_ZR;       % Alarm(:,3) = t until infection    % time being in S, before infected 
                tL = min(TS_ZR);                    % latent period of infector
                tI = max(TS_ZR)-min(TS_ZR);         % infection period of infector
                tRt = ceil(t+TS_ZR);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZR+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            else    % Q
                idv = find(State(idexS_Q,2)==11);
                State(idexS_Q,2) = 10;
                State(idexS_Q,3) = 2;
                State(idexS_Q(idv),3) = 12;
                State(idexS_Q,6) = idxIniLS(i);
                Alarm(idexS_Q,2) = 0; 
                Alarm(idexS_Q,3) = TS_ZQ;
                tL = min(TS_ZQ);
                tI = max(TS_ZQ) - min(TS_ZQ);
                tRt = ceil(t+TS_ZQ);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZQ+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            end

            % Time to progress
            if length(tL) == 0
                Alarm(idxIniLS(i),4) = 0;
            else
                Alarm(idxIniLS(i),4) = tL; 
            end 

            if length(tI) == 0
                Alarm(idxIniLS(i),5) = 0;
            else
                Alarm(idxIniLS(i),5) = tI;
            end  

        end  

        % LV
        for i = 1:NLV0
            SorA = rand();

            if SorA <= k_sym_V      % ISV
                State(idxIniLV(i),4) = 13;
                R0 = R00*(1-eI);
            else                    % IAV
                State(idxIniLV(i),4) = 14;
                R0 = R00*k_reduce*(1-eI);
            end

            % Record next state 
            State(idxIniLV(i),3) = State(idxIniLV(i),4);  

            tc = 5.1;

            p   = (1+(R0)/k)^-1; 
            Z  = nbinrnd(k,p);
            if cc_par==1
                idxCC = CC_matrix{idxIniLV(i)};
                idxCC_toTest = CC_matrix_toTest{idxIniLV(i)};
                if State(idxIniLV(i),4) == 3 % IS
                    Z = Z_all(idxIniLV(i),3);             
                elseif State(idxIniLV(i),4) == 4 % IA
                    Z = Z_all(idxIniLV(i),4);   
                end
            end

            PosRouTestTime  = inf;
            RanTestTime     = inf;  
            PosCCTestTime   = inf;
    %         IsoP            = 0;
            
            VL_since_iso            = [];
            num_test_till_first_neg = [];

            if State(idxIniLV(i),4) == 13 % ISV
                if rand() < symptom_test_pop
                    DetectionTime = tc + turnAround_time_symp;  
            
                    VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                    num_test_till_first_neg = min(find(VL_since_iso < LOD_symp))-1;  
             
                end
            end % 

            %---------routine testing--------
            if ismember(idxIniLV(i),idxRouT) == 1  
                PosRouTestTime = round(routine_test_freq*rand(1) + min(V_all_pos_rou),2) + turnAround_time_rou; 
                DetectionTime = min([DetectionTime PosRouTestTime]); 
            
                VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_rou))-1;  % number of test untill first negative result
             
            end
            %---------------------------------

            %----------random testing---------
            RanTT = rand(size(V_all,1),1) < random_test_capa/N0*dt; 
            RanTT_dt = find(RanTT==1);
            for RanTT_ii=1:length(RanTT_dt)
                if V_all(RanTT_dt(RanTT_ii),2) > LOD_ran
                    RanTestTime = V_all(RanTT_dt(RanTT_ii),1) + turnAround_time_ran;
                    DetectionTime = min([DetectionTime RanTestTime]);
            
                    VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                    num_test_till_first_neg = min(find(VL_since_iso < LOD_ran))-1;  % number of test untill first negative result
            
                end
            end
            %----------------------------------

            %-----------CC testing-------------%
            if cc_par==1
                if (Testing_time(idxIniLV(i)) > 0) && (Testing_time(idxIniLV(i))~=inf)
                    if V_all(V_all(:,1)==Testing_time(idxIniLV(i)),2) > LOD_cc
                        LastCCTestTime = Testing_time(idxIniLV(i));
                        while V_all(V_all(:,1) == round(LastCCTestTime-close_contact_test_freq,2),2) > LOD_cc
                            LastCCTestTime = LastCCTestTime - close_contact_test_freq;
                        end
                        PosCCTestTime = LastCCTestTime + turnAround_time_cc;
                        DetectionTime = min([DetectionTime PosCCTestTime]);
            
                        VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                        num_test_till_first_neg = min(find(VL_since_iso < LOD_cc))-1;  % number of test untill first negative result
            
                    end
                end
            end
            %----------------------------------

            %-------------isolation------------
            IsoP = num_test_till_first_neg*isolation_test_freq;
            date_out_Q = DetectionTime + IsoP;   
            %----------------------------------

            Iso_data = [Iso_data;idxIniLV(i) DetectionTime~=inf DetectionTime {date_out_Q} {IsoP} {num_test_till_first_neg}];

            Testing_data = [Testing_data;idxIniLV(i) ismember(idxIniLV(i),idxRouT) PosRouTestTime ...
               ~isempty(find(RanTT==1)) RanTestTime PosCCTestTime DetectionTime];

            State(idxIniLV(i),7) = Z; 

            % Making infection
            if cc_par==1
                idxInfected = idxCC(randperm(length(idxCC),Z)); 
            else
                idxInfected = randperm(N0,Z);  
            end

            TS_ZQ = [];         TS_ZR = [];      
            iNZQ_S = 0;         iNZR_S = 0;      
            idexS_Q = [];       idexS_R = [];     

            %--------------------------
            % If I => RS , No detection 
            for iZR = 1:Z
                TSrand = randpdf(inf_profile_V(:,2),inf_profile_V(:,1),[1,1]);
                if State(idxInfected(iZR),2) == 1   % S  
                    iNZR_S = iNZR_S+1;
                    idexS_R(iNZR_S) = idxInfected(iZR);
                    TS_ZR(iNZR_S) = TSrand;
                elseif State(idxInfected(iZR),2) == 11  % V
                    pv = rand();
                    if pv > eS
                        iNZR_S = iNZR_S+1; 
                        idexS_R(iNZR_S) = idxInfected(iZR);        
                        TS_ZR(iNZR_S) = TSrand;
                    end
                elseif State(idxInfected(iZR),2) == 10  % PLS
                    iNZR_S = iNZR_S+1; 
                    idexS_R(iNZR_S) = idxInfected(iZR);
                    TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]); 
                elseif State(idxInfected(iZR),2) == 20  % PLV
                    pv = rand();
                    if pv > eS
                        iNZR_S = iNZR_S+1; 
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]);
                    end
                end
            end

    %         %------------------------%
    %         % If I => Q (QSS/QTS)  
            for iZQ = 1:Z
                TSrand = randpdf(inf_profile_V(:,2),inf_profile_V(:,1),[1,1]);
                if State(idxInfected(iZQ),2) == 1   % S  % (infectee is Susceptible)
                    if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                        iNZQ_S = iNZQ_S +1;
                        idexS_Q(iNZQ_S) = idxInfected(iZQ);
                        TS_ZQ(iNZQ_S) = TSrand;
                    end

                elseif State(idxInfected(iZQ),2) == 11   % V an on
                    pv = rand();
                    if pv > eS
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S+1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = TSrand;
                        end
                    end

                elseif State(idxInfected(iZQ),2) == 10   % PLS
                    if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                        iNZQ_S = iNZQ_S +1;
                        idexS_Q(iNZQ_S) = idxInfected(iZQ);
                        TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                    end

                elseif State(idxInfected(iZQ),2) == 20  % PLV
                    pv = rand();
                    if pv > eS
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                        end
                    end   
                end
            end
            %----------------------------------%
            DetectionTime_all(idxIniLV(i))  = DetectionTime;

            if cc_par==1
                idxCC_TestedandInfected = idxCC_toTest(ismember(idxCC_toTest,idexS_Q));
                for ii1=1:length(idxCC_TestedandInfected)
                    Testing_time_temp = inf;
                    if DetectionTime ~= inf
                        Testing_time_temp = round(DetectionTime-TS_ZQ(ii1)+close_contact_test_period-mod(close_contact_test_period,close_contact_test_freq),2);
                    end
                    Testing_time(idxCC_TestedandInfected(ii1)) = min([Testing_time(idxCC_TestedandInfected(ii1)) Testing_time_temp]);
                end     
                CloseContact_data = [CloseContact_data;idxIniLV(i) {idxCC} {idexS_Q} {TS_ZQ} {idxCC_toTest} {idxCC_TestedandInfected}];
            end

            % Secondary infection   
            if cc_par==1
                if DetectionTime == inf 
                    State(idxIniLV(i),5) = 15;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniLV(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                        (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) ||...
                        (isnan(Testing_time(idxIniLV(i)))==0 && DetectionTime==PosCCTestTime+turnAround_time_cc)  % QT
                    State(idxIniLV(i),5) = 17;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniLV(i),5) = 16;
                    Z2 = length(TS_ZQ); 
                end
            else
                if DetectionTime == inf        % not quarantined
                    State(idxIniLV(i),5) = 15;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniLV(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran)
                    State(idxIniLV(i),5) = 17;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniLV(i),5) = 16;
                    Z2 = length(TS_ZQ); 
                end
            end

            State(idxIniLV(i),8) = Z2;

            % Update next state

            if (State(idxIniLV(i),5) == 5) % (UV) L >> R
                idv = find(State(idexS_R,2)==11);   
                State(idexS_R,2) = 10;          % State(:,2) = current state
                State(idexS_R,3) = 2;           % State(:,3) = next state
                State(idexS_R(idv),3) = 12;
                State(idexS_R,6) = idxIniLV(i); % State(:,6) = Spreader
                Alarm(idexS_R,2) = 0;           % Alarm(:,2) = t in current state % still confuse here
                Alarm(idexS_R,3) = TS_ZR;       % Alarm(:,3) = t until infection    % time being in S, before infected 
                tL = min(TS_ZR);                    % latent period of infector
                tI = max(TS_ZR)-min(TS_ZR);         % infection period of infector
                tRt = ceil(t+TS_ZR);
                Rt(tRt,2) = Rt(tRt,2)+1;
                t                = ceil((TS_ZR+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            else    % Q
                idv = find(State(idexS_Q,2)==11);
                State(idexS_Q,2) = 10;
                State(idexS_Q,3) = 2;
                State(idexS_Q(idv),3) = 12;
                State(idexS_Q,6) = idxIniLV(i);
                Alarm(idexS_Q,2) = 0; 
                Alarm(idexS_Q,3) = TS_ZQ;
                tL = min(TS_ZQ);
                tI = max(TS_ZQ) - min(TS_ZQ);
                tRt = ceil(t+TS_ZQ);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZQ+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            end

            % Time to progress
            if length(tL) == 0
                Alarm(idxIniLV(i),4) = 0;
            else
                Alarm(idxIniLV(i),4) = tL; 
            end 

            if length(tI) == 0
                Alarm(idxIniLV(i),5) = 0;
            else
                Alarm(idxIniLV(i),5) = tI;
            end  

        end

        %ISS
        for i = 1:NISS0

            State(idxIniISS(i),4) = 3;
            R0 = R00;

            % Incubation period
            tc = 5.1;

    %         time_since_infect = round(tc*rand(1),2);
            time_since_infect = round((UV_all(end,1)-UV_all(1,1))*rand(1),2);

            p  = (1+(R0)/k)^-1;   
            Z  = nbinrnd(k,p);
            if cc_par==1        
                idxCC = CC_matrix{idxIniISS(i)};
                idxCC_toTest = CC_matrix_toTest{idxIniISS(i)};
                Z = Z_all(idxIniISS(i),1);  
            end

            PosRouTestTime  = inf;
            RanTestTime     = inf;  
            PosCCTestTime   = inf;
    %         IsoP            = 0;
            
            VL_since_iso            = [];
            num_test_till_first_neg = [];

            if rand() < symptom_test_pop
                DetectionTime = tc + turnAround_time_symp; 
            
                VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_symp))-1;  % number of test untill first negative result
              
            end     

            %---------routine testing--------
            if ismember(idxIniISS(i),idxRouT) == 1  
                PosRouTestTime = round(routine_test_freq*rand(1) + ...
                    min(UV_all_pos_rou(UV_all_pos_rou > time_since_infect)),2) + turnAround_time_rou;
                DetectionTime  = min([DetectionTime PosRouTestTime]);
            
                VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_rou))-1;  % number of test untill first negative result
            
            end
            %---------------------------------

            %----------random testing---------
            RanTT = rand(size(UV_all,1),1) < random_test_capa/N0*dt; 
            RanTT_dt = find(RanTT==1 & UV_all(:,1)>time_since_infect);
            for RanTT_ii=1:length(RanTT_dt)
                if UV_all(RanTT_dt(RanTT_ii),2) > LOD_ran
                    RanTestTime = UV_all(RanTT_dt(RanTT_ii),1) + turnAround_time_ran;
                    DetectionTime = min([DetectionTime RanTestTime]);
                end
            
                VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_ran))-1;  % number of test untill first negative result
            
            end
            %----------------------------------

            %----------cc testing---------
            if cc_par==1
                if (Testing_time(idxIniISS(i)) > 0) && (Testing_time(idxIniISS(i))~=inf)
                    if UV_all(UV_all(:,1)==Testing_time(idxIniISS(i)),2) > LOD_cc
                        LastCCTestTime = Testing_time(idxIniISS(i));
                        while UV_all(UV_all(:,1) == round(LastCCTestTime-close_contact_test_freq,2),2) > LOD_cc
                            LastCCTestTime = LastCCTestTime - close_contact_test_freq;
                        end
                        PosCCTestTime = LastCCTestTime + turnAround_time_cc;
                        DetectionTime = min([DetectionTime PosCCTestTime]);
            
                        VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                        num_test_till_first_neg = min(find(VL_since_iso < LOD_cc))-1;  % number of test untill first negative result
            
                    end
                end
            end
            %----------------------------------

            %-------------isolation------------
            IsoP = num_test_till_first_neg*isolation_test_freq; % Isolation period = number of test till first negative * test freq.
            date_out_Q = DetectionTime + IsoP;   
            %----------------------------------

            Iso_data = [Iso_data;idxIniISS(i) DetectionTime~=inf DetectionTime {date_out_Q} {IsoP} {num_test_till_first_neg}];

            Testing_data = [Testing_data;idxIniISS(i) ismember(idxIniISS(i),idxRouT) PosRouTestTime...
               ~isempty(find(RanTT==1)) RanTestTime ...
               PosCCTestTime DetectionTime];

            State(idxIniISS(i),7) = Z; 

            % Making infection
            if cc_par==1
                idxInfected = idxCC(randperm(length(idxCC),Z)); 
            else
                idxInfected = randperm(N0,Z);
            end

            TS_ZQ = [];         TS_ZR = [];      
            iNZQ_S = 0;         iNZR_S = 0;      
            idexS_Q = [];       idexS_R = [];     

            %--------------------------
            % If I => RS 
            for iZR = 1:Z
                TSrand = randpdf(inf_profile_UV(:,2),inf_profile_UV(:,1),[1,1]);
                if TSrand > time_since_infect
                    if State(idxInfected(iZR),2) == 1   % S  
                        iNZR_S = iNZR_S+1;
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = TSrand;
                    elseif State(idxInfected(iZR),2) == 11  % V
                        pv = rand();
                        if pv > eS
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);        
                            TS_ZR(iNZR_S) = TSrand;
                        end
                    elseif State(idxInfected(iZR),2) == 10  % PLS
                        iNZR_S = iNZR_S+1; 
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]); 
                    elseif State(idxInfected(iZR),2) == 20  % PLV
                        pv = rand();
                        if pv > eS
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);
                            TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]);
                        end
                    end
                end
            end

    %         %------------------------%
    %         % If I => Q (QSS/QTS)  
            for iZQ = 1:Z
                TSrand = randpdf(inf_profile_UV(:,2),inf_profile_UV(:,1),[1,1]); 
                if TSrand > time_since_infect
                    if State(idxInfected(iZQ),2) == 1   % S 
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = TSrand;
                        end
                    elseif State(idxInfected(iZQ),2) == 11   % V an on
                        pv = rand();
                        if pv > eS
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S+1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = TSrand;
                            end
                        end
                    elseif State(idxInfected(iZQ),2) == 10   % PLS
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                        end
                    elseif State(idxInfected(iZQ),2) == 20  % PLV
                        pv = rand();
                        if pv > eS
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S +1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                            end
                        end   
                    end
                end
            end
            %----------------------------------%
            DetectionTime_all(idxIniISS(i))  = DetectionTime;

            if cc_par==1
                idxCC_TestedandInfected = idxCC_toTest(ismember(idxCC_toTest,idexS_Q));
                for ii1=1:length(idxCC_TestedandInfected) 
                    Testing_time_temp = inf;
                    if DetectionTime ~= inf
                        Testing_time_temp = round(DetectionTime-TS_ZQ(ii1)+close_contact_test_period-mod(close_contact_test_period,close_contact_test_freq),2);
                    end
                    Testing_time(idxCC_TestedandInfected(ii1)) = min([Testing_time(idxCC_TestedandInfected(ii1)) Testing_time_temp]);
                end     
                CloseContact_data = [CloseContact_data;idxIniISS(i) {idxCC} {idexS_Q} {TS_ZQ} {idxCC_toTest} {idxCC_TestedandInfected}];
            end

            % Secondary infection   
            if cc_par==1
                if DetectionTime == inf 
                    State(idxIniISS(i),5) = 5;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniISS(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) ||...
                (isnan(Testing_time(idxIniISS(i)))==0 && DetectionTime==PosCCTestTime+turnAround_time_cc)  % QT
                    State(idxIniISS(i),5) = 7;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniISS(i),5) = 6;
                    Z2 = length(TS_ZQ); 
                end
            else
                if DetectionTime == inf        % not quarantined
                    State(idxIniISS(i),5) = 5;
                    Z2 = length(TS_ZR);
                elseif (ismember(idxIniISS(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) % QT
                    State(idxIniISS(i),5) = 7;
                    Z2 = length(TS_ZQ);
                else % QS
                    State(idxIniISS(i),5) = 6;
                    Z2 = length(TS_ZQ); 
                end
            end

            % Update next state
            State(idxIniISS(i),3) = State(idxIniISS(i),5);

            % Z2
            State(idxIniISS(i),8) = Z2;

            % Update next state

            if (State(idxIniISS(i),5) == 5) % (UV) L >> R
                idv = find(State(idexS_R,2)==11);   
                State(idexS_R,2) = 10;          % State(:,2) = current state
                State(idexS_R,3) = 2;           % State(:,3) = next state
                State(idexS_R(idv),3) = 12;
                State(idexS_R,6) = idxIniISS(i); % State(:,6) = Spreader
                Alarm(idexS_R,2) = 0;           % Alarm(:,2) = t in current state 
                Alarm(idexS_R,3) = TS_ZR;       % Alarm(:,3) = t until infection    % time being in S, before infected 
                tL = min(TS_ZR);                    % latent period of infector
                tI = max(TS_ZR)-min(TS_ZR);         % infection period of infector
                tRt = ceil(t+TS_ZR);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZR+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            else    % Q
                idv = find(State(idexS_Q,2)==11);
                State(idexS_Q,2) = 10;
                State(idexS_Q,3) = 2;
                State(idexS_Q(idv),3) = 12;
                State(idexS_Q,6) = idxIniISS(i);
                Alarm(idexS_Q,2) = 0; 
                Alarm(idexS_Q,3) = TS_ZQ;
                tL = min(TS_ZQ);
                tI = max(TS_ZQ) - min(TS_ZQ);
                tRt = ceil(t+TS_ZQ);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZQ+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            end

            % Time to progress
            if length(tL) == 0
                Alarm(idxIniISS(i),4) = 0;
            else
                Alarm(idxIniISS(i),4) = tL; 
            end 

            if length(tI) == 0
                Alarm(idxIniISS(i),5) = 0;
            else
                Alarm(idxIniISS(i),5) = tI;
            end  

        end

        %IAS
        for i = 1:NIAS0

            State(idxIniIAS(i),4) = 3;
            R0 = R00*k_reduce;

            % Incubation period
            tc = 5.1;

            time_since_infect = round((UV_all(end,1)-UV_all(1,1))*rand(1),2);

            p  = (1+(R0)/k)^-1;   
            Z  = nbinrnd(k,p);
            if cc_par==1
                idxCC = CC_matrix{idxIniIAS(i)};
                idxCC_toTest = CC_matrix_toTest{idxIniIAS(i)};
                Z = Z_all(idxIniIAS(i),2);  
            end

            PosRouTestTime  = inf;
            RanTestTime     = inf;  
            PosCCTestTime   = inf;
    %         IsoP            = 0;

            DetectionTime = inf;      
            
            VL_since_iso            = [];
            num_test_till_first_neg = []; 

            %---------routine testing--------
            if ismember(idxIniIAS(i),idxRouT) == 1  
                PosRouTestTime = round(routine_test_freq*rand(1) + ...
                    min(UV_all_pos_rou(UV_all_pos_rou > time_since_infect)),2) + turnAround_time_rou;
                if isempty(PosRouTestTime) == 1
                    PosRouTestTime = inf;
                end
                DetectionTime  = min([DetectionTime PosRouTestTime]);
            
                VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_rou))-1;  % number of test untill first negative result
            
            end
            %---------------------------------

            %----------random testing---------
            RanTT = rand(size(UV_all,1),1) < random_test_capa/N0*dt; 
            RanTT_dt = find(RanTT==1 & UV_all(:,1)>time_since_infect);
            for RanTT_ii=1:length(RanTT_dt)
                if UV_all(RanTT_dt(RanTT_ii),2) > LOD_ran
                    RanTestTime = UV_all(RanTT_dt(RanTT_ii),1) + turnAround_time_ran;
                    DetectionTime = min([DetectionTime RanTestTime]);
            
                VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_ran))-1;  % number of test untill first negative result
            
                end
            end
            %----------------------------------

            %----------cc testing---------
            if cc_par==1
                if (Testing_time(idxIniIAS(i)) > 0) && (Testing_time(idxIniIAS(i))~=inf)
                    if UV_all(UV_all(:,1)==Testing_time(idxIniIAS(i)),2) > LOD_cc
                        LastCCTestTime = Testing_time(idxIniIAS(i));
                        while UV_all(UV_all(:,1) == round(LastCCTestTime-close_contact_test_freq,2),2) > LOD_cc
                            LastCCTestTime = LastCCTestTime - close_contact_test_freq;
                        end
                        PosCCTestTime = LastCCTestTime + turnAround_time_cc;
                        DetectionTime = min([DetectionTime PosCCTestTime]);
            
                        VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                        num_test_till_first_neg = min(find(VL_since_iso < LOD_cc))-1;  % number of test untill first negative result
            
                    end
                end
            end
            %----------------------------------

            %-------------isolation------------
            IsoP = num_test_till_first_neg*isolation_test_freq; % Isolation period = number of test till first negative * test freq.
            date_out_Q = DetectionTime + IsoP;   
            %----------------------------------

            Iso_data = [Iso_data;idxIniIAS(i) DetectionTime~=inf DetectionTime {date_out_Q} {IsoP} {num_test_till_first_neg}];

            Testing_data = [Testing_data;idxIniIAS(i) ismember(idxIniIAS(i),idxRouT) PosRouTestTime...
               ~isempty(find(RanTT==1)) RanTestTime ...
               PosCCTestTime DetectionTime];

            State(idxIniIAS(i),7) = Z; 

            % Making infection
            if cc_par==1
                idxInfected = idxCC(randperm(length(idxCC),Z)); 
            else
                idxInfected = randperm(N0,Z);
            end

            TS_ZQ = [];         TS_ZR = [];      
            iNZQ_S = 0;         iNZR_S = 0;      
            idexS_Q = [];       idexS_R = [];     

            %--------------------------
            % If I => RS 
            for iZR = 1:Z
                TSrand = randpdf(inf_profile_UV(:,2),inf_profile_UV(:,1),[1,1]);
                if TSrand > time_since_infect
                    if State(idxInfected(iZR),2) == 1   % S  
                        iNZR_S = iNZR_S+1;
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = TSrand;
                    elseif State(idxInfected(iZR),2) == 11  % V
                        pv = rand();
                        if pv > eS
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);        
                            TS_ZR(iNZR_S) = TSrand;
                        end
                    elseif State(idxInfected(iZR),2) == 10  % PLS
                        iNZR_S = iNZR_S+1; 
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]); 
                    elseif State(idxInfected(iZR),2) == 20  % PLV
                        pv = rand();
                        if pv > eS
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);
                            TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]);
                        end
                    end
                end
            end

    %         %------------------------%
    %         % If I => Q (QSS/QTS)  
            for iZQ = 1:Z
                TSrand = randpdf(inf_profile_UV(:,2),inf_profile_UV(:,1),[1,1]); 
                if TSrand > time_since_infect
                    if State(idxInfected(iZQ),2) == 1   % S 
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = TSrand;
                        end
                    elseif State(idxInfected(iZQ),2) == 11   % V an on
                        pv = rand();
                        if pv > eS
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S+1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = TSrand;
                            end
                        end
                    elseif State(idxInfected(iZQ),2) == 10   % PLS
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                        end
                    elseif State(idxInfected(iZQ),2) == 20  % PLV
                        pv = rand();
                        if pv > eS
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S +1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                            end
                        end   
                    end
                end
            end
            %----------------------------------%
            DetectionTime_all(idxIniIAS(i))  = DetectionTime;

            if cc_par==1
                idxCC_TestedandInfected = idxCC_toTest(ismember(idxCC_toTest,idexS_Q));
                for ii1=1:length(idxCC_TestedandInfected) 
                    Testing_time_temp = inf;
                    if DetectionTime ~= inf
                        Testing_time_temp = round(DetectionTime-TS_ZQ(ii1)+close_contact_test_period-mod(close_contact_test_period,close_contact_test_freq),2);
                    end
                    Testing_time(idxCC_TestedandInfected(ii1)) = min([Testing_time(idxCC_TestedandInfected(ii1)) Testing_time_temp]);
                end     
                CloseContact_data = [CloseContact_data;idxIniIAS(i) {idxCC} {idexS_Q} {TS_ZQ} {idxCC_toTest} {idxCC_TestedandInfected}];
            end

            % Secondary infection  
            if cc_par==1
                if DetectionTime == inf 
                    State(idxIniIAS(i),5) = 5;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniIAS(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) ||...
                (isnan(Testing_time(idxIniIAS(i)))==0 && DetectionTime==PosCCTestTime+turnAround_time_cc)  % QT
                    State(idxIniIAS(i),5) = 7;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniIAS(i),5) = 6;
                    Z2 = length(TS_ZQ); 
                end
            else
                if DetectionTime == inf        % not quarantined
                    State(idxIniIAS(i),5) = 5;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniIAS(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) % QT
                    State(idxIniIAS(i),5) = 7;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniIAS(i),5) = 6;
                    Z2 = length(TS_ZQ); 
                end   
            end
            % Update next state
            State(idxIniIAS(i),3) = State(idxIniIAS(i),5);

            % Z2
            State(idxIniIAS(i),8) = Z2;

            % Update next state

            if (State(idxIniIAS(i),5) == 5) % (UV) L >> R
                idv = find(State(idexS_R,2)==11);   
                State(idexS_R,2) = 10;          % State(:,2) = current state
                State(idexS_R,3) = 2;           % State(:,3) = next state
                State(idexS_R(idv),3) = 12;
                State(idexS_R,6) = idxIniIAS(i); % State(:,6) = Spreader
                Alarm(idexS_R,2) = 0;           % Alarm(:,2) = t in current state % still confuse here
                Alarm(idexS_R,3) = TS_ZR;       % Alarm(:,3) = t until infection    % time being in S, before infected 
                tL = min(TS_ZR);                    % latent period of infector
                tI = max(TS_ZR)-min(TS_ZR);         % infection period of infector
                tRt = ceil(t+TS_ZR);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZR+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            else    % Q
                idv = find(State(idexS_Q,2)==11);
                State(idexS_Q,2) = 10;
                State(idexS_Q,3) = 2;
                State(idexS_Q(idv),3) = 12;
                State(idexS_Q,6) = idxIniIAS(i);
                Alarm(idexS_Q,2) = 0; 
                Alarm(idexS_Q,3) = TS_ZQ;
                tL = min(TS_ZQ);
                tI = max(TS_ZQ) - min(TS_ZQ);
                tRt = ceil(t+TS_ZQ);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZQ+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            end

            % Time to progress
            if length(tL) == 0
                Alarm(idxIniIAS(i),4) = 0;
            else
                Alarm(idxIniIAS(i),4) = tL; 
            end 

            if length(tI) == 0
                Alarm(idxIniIAS(i),5) = 0;
            else
                Alarm(idxIniIAS(i),5) = tI;
            end  

        end

        %ISV
        for i = 1:NISV0

            State(idxIniISV(i),4) = 13;
            R0 = R00*(1-eI);

            % Record next state 
    %         State(idxIniISS(i),3) = State(idxIniISS(i),4);  

            % Incubation period
            tc = 5.1;

    %         time_since_infect = round(tc*rand(1),2);
            time_since_infect = round((V_all(end,1)-V_all(1,1))*rand(1),2);

            p  = (1+(R0)/k)^-1;   
            Z  = nbinrnd(k,p);
            if cc_par==1
                idxCC = CC_matrix{idxIniISV(i)};
                idxCC_toTest = CC_matrix_toTest{idxIniISV(i)};
                Z = Z_all(idxIniISV(i),3);     
            end

            PosRouTestTime  = inf;
            RanTestTime     = inf;  
            PosCCTestTime   = inf;
    %         IsoP            = 0;
            
            VL_since_iso            = [];
            num_test_till_first_neg = [];

    %         DetectionTime = tc;   
            if rand() < symptom_test_pop
                DetectionTime = tc + turnAround_time_symp; 
            
                VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_symp))-1;  % number of test untill first negative result
                        
            end        

            %---------routine testing--------
            if ismember(idxIniISV(i),idxRouT) == 1  
                PosRouTestTime = round(routine_test_freq*rand(1) + ...
                    min(V_all_pos_rou(V_all_pos_rou > time_since_infect)),2) + turnAround_time_rou;
                DetectionTime  = min([DetectionTime PosRouTestTime]);
            
                VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_rou))-1;  % number of test untill first negative result
            
            end
            %---------------------------------

            %----------random testing---------
            RanTT = rand(size(V_all,1),1) < random_test_capa/N0*dt; 
            RanTT_dt = find(RanTT==1 & V_all(:,1)>time_since_infect);
            for RanTT_ii=1:length(RanTT_dt)
                if V_all(RanTT_dt(RanTT_ii),2) > LOD_ran
                    RanTestTime = V_all(RanTT_dt(RanTT_ii),1) + turnAround_time_ran;
                    DetectionTime = min([DetectionTime RanTestTime]);
            
                    VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                    num_test_till_first_neg = min(find(VL_since_iso < LOD_ran))-1;  % number of test untill first negative result
            
                end
            end
            %----------------------------------

            %----------cc testing---------
            if cc_par==1
                if (Testing_time(idxIniISV(i)) > 0) && (Testing_time(idxIniISV(i))~=inf)
                    if V_all(V_all(:,1)==Testing_time(idxIniISV(i)),2) > LOD_cc
                        LastCCTestTime = Testing_time(idxIniISV(i));
                        while V_all(V_all(:,1) == round(LastCCTestTime-close_contact_test_freq,2),2) > LOD_cc
                            LastCCTestTime = LastCCTestTime - close_contact_test_freq;
                        end
                        PosCCTestTime = LastCCTestTime + turnAround_time_cc;
                        DetectionTime = min([DetectionTime PosCCTestTime]);
            
                        VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                        num_test_till_first_neg = min(find(VL_since_iso < LOD_cc))-1;  % number of test untill first negative result
            
                    end
                end
            end
            %----------------------------------

            %-------------isolation------------
            IsoP = num_test_till_first_neg*isolation_test_freq; 
            date_out_Q = DetectionTime + IsoP;   
            %----------------------------------

            Iso_data = [Iso_data;idxIniISV(i) DetectionTime~=inf DetectionTime {date_out_Q} {IsoP} {num_test_till_first_neg}];

            Testing_data = [Testing_data;idxIniISV(i) ismember(idxIniISV(i),idxRouT) PosRouTestTime...
               ~isempty(find(RanTT==1)) RanTestTime ...
               PosCCTestTime DetectionTime];

            State(idxIniISV(i),7) = Z; 

            % Making infection
            if cc_par==1
                idxInfected = idxCC(randperm(length(idxCC),Z)); 
            else    
                idxInfected = randperm(N0,Z); 
            end

            TS_ZQ = [];         TS_ZR = [];      
            iNZQ_S = 0;         iNZR_S = 0;      
            idexS_Q = [];       idexS_R = [];     

            %--------------------------
            % If I => RS 
            for iZR = 1:Z
                TSrand = randpdf(inf_profile_V(:,2),inf_profile_V(:,1),[1,1]);
                if TSrand > time_since_infect
                    if State(idxInfected(iZR),2) == 1   % S  
                        iNZR_S = iNZR_S+1;
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = TSrand;
                    elseif State(idxInfected(iZR),2) == 11  % V
                        pv = rand();
                        if pv > eS
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);        
                            TS_ZR(iNZR_S) = TSrand;
                        end
                    elseif State(idxInfected(iZR),2) == 10  % PLS
                        iNZR_S = iNZR_S+1; 
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]); 
                    elseif State(idxInfected(iZR),2) == 20  % PLV
                        pv = rand();
                        if pv > eS
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);
                            TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]);
                        end
                    end
                end
            end

    %         %------------------------%
    %         % If I => Q (QSS/QTS)  
            for iZQ = 1:Z
                TSrand = randpdf(inf_profile_V(:,2),inf_profile_V(:,1),[1,1]); 
                if TSrand > time_since_infect
                    if State(idxInfected(iZQ),2) == 1   % S 
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = TSrand;
                        end
                    elseif State(idxInfected(iZQ),2) == 11   % V an on
                        pv = rand();
                        if pv > eS
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S+1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = TSrand;
                            end
                        end
                    elseif State(idxInfected(iZQ),2) == 10   % PLS
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                        end
                    elseif State(idxInfected(iZQ),2) == 20  % PLV
                        pv = rand();
                        if pv > eS
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S +1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                            end
                        end   
                    end
                end
            end
            %----------------------------------%
            DetectionTime_all(idxIniISV(i))  = DetectionTime;

            if cc_par==1
                idxCC_TestedandInfected = idxCC_toTest(ismember(idxCC_toTest,idexS_Q));
                for ii1=1:length(idxCC_TestedandInfected) 
                    Testing_time_temp = inf;
                    if DetectionTime ~= inf
                        Testing_time_temp = round(DetectionTime-TS_ZQ(ii1)+close_contact_test_period-mod(close_contact_test_period,close_contact_test_freq),2);
                    end
                    Testing_time(idxCC_TestedandInfected(ii1)) = min([Testing_time(idxCC_TestedandInfected(ii1)) Testing_time_temp]);
                end 
                CloseContact_data = [CloseContact_data;idxIniISV(i) {idxCC} {idexS_Q} {TS_ZQ} {idxCC_toTest} {idxCC_TestedandInfected}];
            end

            % Secondary infection    
            if cc_par==1
                if DetectionTime == inf 
                    State(idxIniISV(i),5) = 15;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniISV(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) ||...
                (isnan(Testing_time(idxIniISV(i)))==0 && DetectionTime==PosCCTestTime+turnAround_time_cc)  % QT
                    State(idxIniISV(i),5) = 17;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniISV(i),5) = 16;
                    Z2 = length(TS_ZQ); 
                end
            else
                if DetectionTime == inf        % not quarantined
                    State(idxIniISV(i),5) = 15;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniISV(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) % QT
                    State(idxIniISV(i),5) = 17;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniISV(i),5) = 16;
                    Z2 = length(TS_ZQ); 
                end
            end

            % Update next state
            State(idxIniISV(i),3) = State(idxIniISV(i),5);

            % Z2
            State(idxIniISV(i),8) = Z2;

            % Update next state

            if (State(idxIniISV(i),5) == 15) % (UV) L >> R
                idv = find(State(idexS_R,2)==11);   
                State(idexS_R,2) = 10;          % State(:,2) = current state
                State(idexS_R,3) = 2;           % State(:,3) = next state
                State(idexS_R(idv),3) = 12;
                State(idexS_R,6) = idxIniISV(i); % State(:,6) = Spreader
                Alarm(idexS_R,2) = 0;           % Alarm(:,2) = t in current state % still confuse here
                Alarm(idexS_R,3) = TS_ZR;       % Alarm(:,3) = t until infection    % time being in S, before infected 
                tL = min(TS_ZR);                    % latent period of infector
                tI = max(TS_ZR)-min(TS_ZR);         % infection period of infector
                tRt = ceil(t+TS_ZR);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZR+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            else    % Q
                idv = find(State(idexS_Q,2)==11);
                State(idexS_Q,2) = 10;
                State(idexS_Q,3) = 2;
                State(idexS_Q(idv),3) = 12;
                State(idexS_Q,6) = idxIniISV(i);
                Alarm(idexS_Q,2) = 0; 
                Alarm(idexS_Q,3) = TS_ZQ;
                tL = min(TS_ZQ);
                tI = max(TS_ZQ) - min(TS_ZQ);
                tRt = ceil(t+TS_ZQ);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZQ+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            end

            % Time to progress
            if length(tL) == 0
                Alarm(idxIniISV(i),4) = 0;
            else
                Alarm(idxIniISV(i),4) = tL; 
            end 

            if length(tI) == 0
                Alarm(idxIniISV(i),5) = 0;
            else
                Alarm(idxIniISV(i),5) = tI;
            end  
        end

        % IAV
        for i = 1:NIAV0

            State(idxIniIAV(i),4) = 13;
            R0 = R00*k_reduce*(1-eI);

            % Incubation period
            tc = 5.1;

            time_since_infect = round((V_all(end,1)-V_all(1,1))*rand(1),2);

            p  = (1+(R0)/k)^-1;   
            Z  = nbinrnd(k,p);
            if cc_par==1
                idxCC = CC_matrix{idxIniIAV(i)};
                idxCC_toTest = CC_matrix_toTest{idxIniIAV(i)};
                Z = Z_all(idxIniIAV(i),4);             
            end

            PosRouTestTime  = inf;
            RanTestTime     = inf;  
            PosCCTestTime   = inf;
    %         IsoP            = 0;

            DetectionTime = inf;       
            
            VL_since_iso            = [];
            num_test_till_first_neg = [];

            %---------routine testing--------
            if ismember(idxIniIAV(i),idxRouT) == 1  
                PosRouTestTime = round(routine_test_freq*rand(1) + ...
                    min(V_all_pos_rou(V_all_pos_rou > time_since_infect)),2) + turnAround_time_rou;
                if isempty(PosRouTestTime) == 1
                    PosRouTestTime = inf;
                end
                DetectionTime  = min([DetectionTime PosRouTestTime]);
            
                VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                num_test_till_first_neg = min(find(VL_since_iso < LOD_rou))-1;  % number of test untill first negative result
            
            end
            %---------------------------------

            %----------random testing---------
            RanTT = rand(size(V_all,1),1) < random_test_capa/N0*dt; 
            RanTT_dt = find(RanTT==1 & V_all(:,1)>time_since_infect);
            for RanTT_ii=1:length(RanTT_dt)
                if V_all(RanTT_dt(RanTT_ii),2) > LOD_ran
                    RanTestTime = V_all(RanTT_dt(RanTT_ii),1) + turnAround_time_ran;
                    DetectionTime = min([DetectionTime RanTestTime]);
            
                    VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                    num_test_till_first_neg = min(find(VL_since_iso < LOD_ran))-1;  % number of test untill first negative result
            
                end
            end
            %----------------------------------

            %----------cc testing---------
            if cc_par==1
                if (Testing_time(idxIniIAV(i)) > 0) && (Testing_time(idxIniIAV(i))~=inf)
                    if V_all(V_all(:,1)==Testing_time(idxIniIAV(i)),2) > LOD_cc
                        LastCCTestTime = Testing_time(idxIniIAV(i));
                        while V_all(V_all(:,1) == round(LastCCTestTime-close_contact_test_freq,2),2) > LOD_cc
                            LastCCTestTime = LastCCTestTime - close_contact_test_freq;
                        end
                        PosCCTestTime = LastCCTestTime + turnAround_time_cc;
                        DetectionTime = min([DetectionTime PosCCTestTime]);
            
                        VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                        num_test_till_first_neg = min(find(VL_since_iso < LOD_cc))-1;  % number of test untill first negative result
            
                    end
                end
            end
            %----------------------------------

            %-------------isolation------------
            IsoP = num_test_till_first_neg*isolation_test_freq; % Isolation period = number of test till first negative * test freq.
            date_out_Q = DetectionTime + IsoP;   
            %----------------------------------

            Iso_data = [Iso_data;idxIniIAV(i) DetectionTime~=inf DetectionTime {date_out_Q} {IsoP} {num_test_till_first_neg}];

            Testing_data = [Testing_data;idxIniIAV(i) ismember(idxIniIAV(i),idxRouT) PosRouTestTime...
               ~isempty(find(RanTT==1)) RanTestTime ...
               PosCCTestTime DetectionTime];

            State(idxIniIAV(i),7) = Z; 

            % Making infection
            if cc_par==1
                idxInfected = idxCC(randperm(length(idxCC),Z)); 
            else
                idxInfected = randperm(N0,Z); 
            end

            TS_ZQ = [];         TS_ZR = [];      
            iNZQ_S = 0;         iNZR_S = 0;      
            idexS_Q = [];       idexS_R = [];     

            %--------------------------
            % If I => RS 
            for iZR = 1:Z
                TSrand = randpdf(inf_profile_V(:,2),inf_profile_V(:,1),[1,1]);
                if TSrand > time_since_infect
                    if State(idxInfected(iZR),2) == 1   % S  
                        iNZR_S = iNZR_S+1;
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = TSrand;
                    elseif State(idxInfected(iZR),2) == 11  % V
                        pv = rand();
                        if pv > eS
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);        
                            TS_ZR(iNZR_S) = TSrand;
                        end
                    elseif State(idxInfected(iZR),2) == 10  % PLS
                        iNZR_S = iNZR_S+1; 
                        idexS_R(iNZR_S) = idxInfected(iZR);
                        TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]); 
                    elseif State(idxInfected(iZR),2) == 20  % PLV
                        pv = rand();
                        if pv > eS
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);
                            TS_ZR(iNZR_S) = min([TSrand  Alarm(idxInfected(iZR),3)]);
                        end
                    end
                end
            end

    %         %------------------------%
    %         % If I => Q (QSS/QTS)  
            for iZQ = 1:Z
                TSrand = randpdf(inf_profile_V(:,2),inf_profile_V(:,1),[1,1]); 
                if TSrand > time_since_infect
                    if State(idxInfected(iZQ),2) == 1   % S 
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = TSrand;
                        end
                    elseif State(idxInfected(iZQ),2) == 11   % V an on
                        pv = rand();
                        if pv > eS
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S+1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = TSrand;
                            end
                        end
                    elseif State(idxInfected(iZQ),2) == 10   % PLS
                        if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                            iNZQ_S = iNZQ_S +1;
                            idexS_Q(iNZQ_S) = idxInfected(iZQ);
                            TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                        end
                    elseif State(idxInfected(iZQ),2) == 20  % PLV
                        pv = rand();
                        if pv > eS
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S +1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = min([TSrand Alarm(idxInfected(iZQ),3)]);
                            end
                        end   
                    end
                end
            end
            %----------------------------------%
            DetectionTime_all(idxIniIAV(i))  = DetectionTime;

            if cc_par==1
                idxCC_TestedandInfected = idxCC_toTest(ismember(idxCC_toTest,idexS_Q));
                for ii1=1:length(idxCC_TestedandInfected) 
                    Testing_time_temp = inf;
                    if DetectionTime ~= inf
                        Testing_time_temp = round(DetectionTime-TS_ZQ(ii1)+close_contact_test_period-mod(close_contact_test_period,close_contact_test_freq),2);
                    end
                    Testing_time(idxCC_TestedandInfected(ii1)) = min([Testing_time(idxCC_TestedandInfected(ii1)) Testing_time_temp]);
                end     
                CloseContact_data = [CloseContact_data;idxIniIAV(i) {idxCC} {idexS_Q} {TS_ZQ} {idxCC_toTest} {idxCC_TestedandInfected}];
            end

            % Secondary infection   
            if cc_par==1
                if DetectionTime == inf 
                    State(idxIniIAV(i),5) = 15;
                    Z2 = length(TS_ZR);
                elseif (ismember(idxIniIAV(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) ||...
                (isnan(Testing_time(idxIniIAV(i)))==0 && DetectionTime==PosCCTestTime+turnAround_time_cc)  % QT
                    State(idxIniIAV(i),5) = 17;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniIAV(i),5) = 16;
                    Z2 = length(TS_ZQ); 
                end
            else
                if DetectionTime == inf        % not quarantined
                    State(idxIniIAV(i),5) = 15;
                    Z2 = length(TS_ZR); 
                elseif (ismember(idxIniIAV(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) % QT
                    State(idxIniIAV(i),5) = 17;
                    Z2 = length(TS_ZQ); 
                else % QS
                    State(idxIniIAV(i),5) = 16;
                    Z2 = length(TS_ZQ); 
                end
            end

            % Update next state
            State(idxIniIAV(i),3) = State(idxIniIAV(i),5);

            % Z2
            State(idxIniIAV(i),8) = Z2;

            % Update next state

            if (State(idxIniIAV(i),5) == 15) % (UV) L >> R
                idv = find(State(idexS_R,2)==11);   
                State(idexS_R,2) = 10;          % State(:,2) = current state
                State(idexS_R,3) = 2;           % State(:,3) = next state
                State(idexS_R(idv),3) = 12;
                State(idexS_R,6) = idxIniIAV(i); % State(:,6) = Spreader
                Alarm(idexS_R,2) = 0;           % Alarm(:,2) = t in current state % still confuse here
                Alarm(idexS_R,3) = TS_ZR;       % Alarm(:,3) = t until infection    % time being in S, before infected 
                tL = min(TS_ZR);                    % latent period of infector
                tI = max(TS_ZR)-min(TS_ZR);         % infection period of infector
                tRt = ceil(t+TS_ZR);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZR+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            else    % Q
                idv = find(State(idexS_Q,2)==11);
                State(idexS_Q,2) = 10;
                State(idexS_Q,3) = 2;
                State(idexS_Q(idv),3) = 12;
                State(idexS_Q,6) = idxIniIAV(i);
                Alarm(idexS_Q,2) = 0; 
                Alarm(idexS_Q,3) = TS_ZQ;
                tL = min(TS_ZQ);
                tI = max(TS_ZQ) - min(TS_ZQ);
                tRt = ceil(t+TS_ZQ);
                Rt(tRt,2) = Rt(tRt,2)+1;
                tRt2 = ceil((TS_ZQ+t)/dt);
                Rt2(tRt2,2) = Rt2(tRt2,2)+1; 
            end

            % Time to progress
            if length(tL) == 0
                Alarm(idxIniIAV(i),4) = 0;
            else
                Alarm(idxIniIAV(i),4) = tL; 
            end 

            if length(tI) == 0
                Alarm(idxIniIAV(i),5) = 0;
            else
                Alarm(idxIniIAV(i),5) = tI;
            end  

        end
        
%         idx_seeded_save=[];

    %% Time evolution

        for istep = 1:Nstep
            
            
           idS = find(State(:,2) == 1);     SS = State(idS,:);      AS = Alarm(idS,:);      NS = length(SS(:,1));
           idV = find(State(:,2) == 11);    SV = State(idV,:);      AV = Alarm(idV,:);      NV = length(SV(:,1));
           
           %----import case----%
           if  mod(t,7) < dt
               day_seed = round(t+rand(num_import_case,1)*7,2);
           end
           
           if sum(round(t,2)==day_seed)>0
               idx_seeded = randperm(NS+NV,sum(round(t,2)==day_seed));
               if idx_seeded <= NS
                   State(idx_seeded,2) = 10; % update state for check
                   State(idx_seeded,3) = 2; 
                   State(idx_seeded,6) = -4;
                   Alarm(idx_seeded,3) = day_seed(day_seed==round(t,2))-dt;
               else
                   State(idV(idx_seeded-NS),2) = 20;
                   State(idV(idx_seeded-NS),3) = 12;
                   State(idV(idx_seeded-NS),6) = -4;
                   Alarm(idV(idx_seeded-NS),3) = day_seed(day_seed==round(t,2))-dt;
               end
%                idx_seeded_save = [idx_seeded_save idx_seeded];
              
               idS = find(State(:,2) == 1);     SS = State(idS,:);      AS = Alarm(idS,:);      NS = length(SS(:,1)); %reupdate S,V
               idV = find(State(:,2) == 11);    SV = State(idV,:);      AV = Alarm(idV,:);      NV = length(SV(:,1));
           end
           %------------------%
           
           idLS = find(State(:,2) == 2);    SLS = State(idLS,:);    ALS = Alarm(idLS,:);    NLS = length(SLS(:,1));
           idISS = find(State(:,2) == 3);   SISS = State(idISS,:);  AISS = Alarm(idISS,:);  NISS = length(SISS(:,1));    
           idIAS = find(State(:,2) == 4);   SIAS = State(idIAS,:);  AIAS = Alarm(idIAS,:);  NIAS= length(SIAS(:,1));

           idLV = find(State(:,2) == 12);   SLV = State(idLV,:);    ALV = Alarm(idLV,:);    NLV = length(SLV(:,1));
           idISV = find(State(:,2) == 13);  SISV = State(idISV,:);  AISV = Alarm(idISV,:);  NISV = length(SISV(:,1));    
           idIAV = find(State(:,2) == 14);  SIAV = State(idIAV,:);  AIAV = Alarm(idIAV,:);  NIAV= length(SIAV(:,1));

           idRS = find(State(:,2) == 5);    SRS = State(idRS,:);    ARS = Alarm(idRS,:);    NRS = length(SRS(:,1));
           idRV = find(State(:,2) == 15);   SRV = State(idRV,:);    ARV = Alarm(idRV,:);    NRV = length(SRV(:,1));

           idQSS = find(State(:,2) == 6);   SQSS = State(idQSS,:);  AQSS = Alarm(idQSS,:);  NQSS = length(SQSS(:,1));
           idQTS = find(State(:,2) == 7);   SQTS = State(idQTS,:);  AQTS = Alarm(idQTS,:);  NQTS = length(SQTS(:,1));
           idQSV = find(State(:,2) == 16);  SQSV = State(idQSV,:);  AQSV = Alarm(idQSV,:);  NQSV = length(SQSV(:,1));
           idQTV = find(State(:,2) == 17);  SQTV = State(idQTV,:);  AQTV = Alarm(idQTV,:);  NQTV = length(SQTV(:,1));

           idPLS = find(State(:,2) == 10);   SPLS = State(idPLS,:);  APLS = Alarm(idPLS,:);  NPLS = length(SPLS(:,1));
           idPLV = find(State(:,2) == 20);   SPLV = State(idPLV,:);  APLV = Alarm(idPLV,:);  NPLV = length(SPLV(:,1));

           N = NS+NLS+NISS+NIAS+NV+NLV+NISV+NIAV+NRS+NRV+NQSS+NQTS+NQSV+NQTV+NPLS+NPLV;
           if N > N0
               break;
           end

           % Record data
           if num_import_case == 0
               
               if (NLS+NISS+NIAS+NLV+NISV+NIAV+NQSS+NQTS+NQSV+NQTV+NPLS+NPLV) == 0 % check extinct
                   day = day+1;
                   tsave(day:end,irepeat) = ceil(t-dt):tf-1;
                   Ssave(day:end,irepeat) = NS;         LSsave(day:end,irepeat) = NLS;
                   ISSsave(day:end,irepeat) = NISS;     IASsave(day:end,irepeat) = NIAS;

                   Vsave(day:end,irepeat) = NV;         LVsave(day:end,irepeat) = NLV;
                   ISVsave(day:end,irepeat) = NISV;     IAVsave(day:end,irepeat) = NIAV;
                   Iallsave(day:end,irepeat) = NISS+NIAS+NISV+NIAV;

                   RSsave(day:end,irepeat) = NRS;       RVsave(day:end,irepeat) = NRV; 
                   Rallsave(day:end,irepeat) = NRS+NRV;             

                   QSSsave(day:end,irepeat) = NQSS;     QTSsave(day:end,irepeat) = NQTS;
                   QSVsave(day:end,irepeat) = NQSV;     QTVsave(day:end,irepeat) = NQTV; 
                   Qallsave(day:end,irepeat) = NQSS+NQTS+NQSV+NQTV;     

                   PLSsave(day:end,irepeat) = NPLS;     PLVsave(day:end,irepeat) = NPLV;

                   CumISSsave(day:end,irepeat) = CumISS;    CumIASsave(day:end,irepeat) = CumIAS;
                   CumISVsave(day:end,irepeat) = CumISV;    CumIAVsave(day:end,irepeat) = CumIAV;
                   CumISsave(day:end,irepeat) = CumISS+CumISV; 
                   CumIAsave(day:end,irepeat) = CumIAS+CumIAV;
                   CumIallsave(day:end,irepeat) = CumISS+CumISV+CumIAS+CumIAV;
                   CumRSsave(day:end,irepeat) = CumRS;      CumRVsave(day:end,irepeat) = CumRV; 
                   CumRsave(day:end,irepeat) = CumRS+CumRV;
                   CumQSSsave(day:end,irepeat) = CumQSS;    CumQTSsave(day:end,irepeat) = CumQTS;
                   CumQSVsave(day:end,irepeat) = CumQSV;    CumQTVsave(day:end,irepeat) = CumQTV;
                   CumXsave(day:end,irepeat) = CumQSS+CumQTS+CumQSV+CumQTV;
                   display([iBatch irepeat day]) 
                   break
               end
           end

           if mod(t,dtsave) < dt
               day = day+1;
               tsave(day,irepeat) = t;
               Ssave(day,irepeat) = NS;         LSsave(day,irepeat) = NLS;
               ISSsave(day,irepeat) = NISS;     IASsave(day,irepeat) = NIAS;

               Vsave(day,irepeat) = NV;         LVsave(day,irepeat) = NLV;
               ISVsave(day,irepeat) = NISV;     IAVsave(day,irepeat) = NIAV;
               Iallsave(day:end,irepeat) = NISS+NIAS+NISV+NIAV;

               RSsave(day,irepeat) = NRS;       RVsave(day,irepeat) = NRV;  
               Rallsave(day,irepeat) = NRS+NRV;     

               QSSsave(day,irepeat) = NQSS;     QTSsave(day,irepeat) = NQTS;
               QSVsave(day,irepeat) = NQSV;     QTVsave(day,irepeat) = NQTV; 
               Qallsave(day,irepeat) = NQSS+NQTS+NQSV+NQTV;     

               PLSsave(day,irepeat) = NPLS;     PLVsave(day,irepeat) = NPLV;

               CumISSsave(day,irepeat) = CumISS;  CumIASsave(day,irepeat) = CumIAS;
               CumISVsave(day,irepeat) = CumISV;  CumIAVsave(day,irepeat) = CumIAV;
               CumISsave(day,irepeat) = CumISS+CumISV;  
               CumIAsave(day,irepeat) = CumIAS+CumIAV;
               CumIallsave(day,irepeat) = CumISS+CumIAS+CumISV+CumIAV;
               CumRSsave(day,irepeat) = CumRS;      CumRVsave(day,irepeat) = CumRV;  
               CumRsave(day,irepeat) = CumRS+CumRV;  
               CumQSSsave(day,irepeat) = CumQSS;    CumQTSsave(day,irepeat) = CumQTS;
               CumQSVsave(day,irepeat) = CumQSV;    CumQTVsave(day,irepeat) = CumQTV;
               CumXsave(day,irepeat) = CumQSS+CumQTS+CumQSV+CumQTV;
               if mod(day,10)==0
                   display([iBatch irepeat day])
               end
           end

           %% Update state
           % PLS => LS

           LS_count = 0;
           for i = 1:NPLS
               if APLS(i,2) > APLS(i,3)
                   LS_count = LS_count+1;
                   SPLS(i,2) = SPLS(i,3);     % Update state
                   tc = 5.1;
                   APLS(i,2) = APLS(i,2) - APLS(i,3);

                   % Next state
                    SorA = rand();
                    if SorA <= k_sym_S     % IS
                        SPLS(i,4) = 3;
                        R0 = R00;
                    else                % IA
                        SPLS(i,4) = 4;
                        R0 = R00*k_reduce;
                    end

                    SPLS(i,3) = SPLS(i,4);  % Update next state

                    p = (1+R0/k)^-1 ;
                    Z = nbinrnd(k,p);
                    if cc_par==1
                        idxCC = CC_matrix{idPLS(i)};
                        idxCC_toTest = CC_matrix_toTest{idPLS(i)};
                        if SPLS(i,4) == 3 % IS
                            Z = Z_all(idPLS(i),1);             
                        elseif SPLS(i,4) == 4 % IA
                            Z = Z_all(idPLS(i),2);   
                        end
                    end

                    PosRouTestTime  = inf;
                    RanTestTime     = inf;  
                    PosCCTestTime   = inf;
    %                 IsoP            = 0;

                    DetectionTime = inf;
            
                    VL_since_iso            = [];
                    num_test_till_first_neg = [];

                    if SPLS(i,4) == 3 % IS
    %                     DetectionTime = tc;  
                        if rand() < symptom_test_pop
                            DetectionTime = tc + turnAround_time_symp;  
            
                            VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                            num_test_till_first_neg = min(find(VL_since_iso < LOD_symp))-1;  % number of test untill first negative result
             
                        end         
                    end 

                    %---------routine testing--------
                    if ismember(idPLS(i),idxRouT) == 1  
                        PosRouTestTime = round(routine_test_freq*rand(1) + min(UV_all_pos_rou),2) + turnAround_time_rou;
                        DetectionTime  = min([DetectionTime PosRouTestTime]);
            
                        VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                        num_test_till_first_neg = min(find(VL_since_iso < LOD_rou))-1;  % number of test untill first negative result
            
                    end
                    %---------------------------------

                    %----------random testing---------
                    RanTT = rand(size(UV_all,1),1) < random_test_capa/N0*dt; 
                    RanTT_dt = find(RanTT==1);
                    for RanTT_ii=1:length(RanTT_dt)
                        if UV_all(RanTT_dt(RanTT_ii),2) > LOD_ran
                            RanTestTime = UV_all(RanTT_dt(RanTT_ii),1) + turnAround_time_ran;
                            DetectionTime = min([DetectionTime RanTestTime]);
            
                            VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                            num_test_till_first_neg = min(find(VL_since_iso < LOD_ran))-1;  % number of test untill first negative result
            
                        end
                    end
                    %----------------------------------

                    %------------cc testing------------
                    if cc_par==1
                        if (Testing_time(idPLS(i)) > 0) && (Testing_time(idPLS(i))~=inf)
                            if UV_all(UV_all(:,1)==Testing_time(idPLS(i)),2) > LOD_cc
                                LastCCTestTime = Testing_time(idPLS(i));
                                while UV_all(UV_all(:,1) == round(LastCCTestTime-close_contact_test_freq,2),2) > LOD_cc
                                    LastCCTestTime = LastCCTestTime - close_contact_test_freq;
                                end
                                PosCCTestTime = LastCCTestTime + turnAround_time_cc;
                                DetectionTime = min([DetectionTime PosCCTestTime]);
            
                                VL_since_iso = UV_all(find(ismember(UV_all(:,1),round([DetectionTime:isolation_test_freq:UV_all(end,1)],2))),2);
                                num_test_till_first_neg = min(find(VL_since_iso < LOD_cc))-1;  % number of test untill first negative result
            
                            end
                        end
                    end
                    %----------------------------------
                    
                    %-------------isolation------------
                    IsoP = num_test_till_first_neg*isolation_test_freq;
                    date_out_Q = DetectionTime + IsoP;
                    %----------------------------------

                    Iso_data = [Iso_data;idPLS(i) DetectionTime~=inf ...
                        DetectionTime {date_out_Q} {IsoP} {num_test_till_first_neg}];

                    Testing_data = [Testing_data; idPLS(i) ismember(idPLS(i),idxRouT) PosRouTestTime...
                       ~isempty(find(RanTT==1)) RanTestTime ...
                        PosCCTestTime DetectionTime];

                    SPLS(i,7) = Z; 

                    % Making infection
                    if cc_par==1
                        idxInfected = idxCC(randperm(length(idxCC),Z)); 
                    else
                        idxInfected = randperm(N,Z);   
                    end

                    TS_ZQ = [];     TS_ZR = [];     % ts for each infectee
                    TV_ZQ = [];     TV_ZR = [];
                    iNZQ_S = 0;     iNZR_S = 0;       % N successful infection
                    iNZQ_V = 0;     iNZR_V = 0;
                    idexS_Q = [];   idexS_R = [];   % index of infectee (successful infection)
                    idexV_Q = [];   idexV_R = [];    

                    % If I => R                    
                    for iZR = 1:Z
                        TSrand = randpdf(inf_profile_UV(:,2),inf_profile_UV(:,1),[1,1]);
                        if SS(ismember(SS(:,1),idxInfected(iZR))) == idxInfected(iZR)
                            iNZR_S = iNZR_S+1;
                            idexS_R(iNZR_S) = idxInfected(iZR);
                            TS_ZR(iNZR_S) = TSrand;
                        elseif SV(ismember(SV(:,1),idxInfected(iZR))) == idxInfected(iZR)
                            pv = rand();
                            if pv > eS
                                iNZR_V = iNZR_V+1; 
                                idexV_R(iNZR_V) = idxInfected(iZR);        
                                TV_ZR(iNZR_V) = TSrand;
                            end
                        elseif SPLS(ismember(SPLS(:,1),idxInfected(iZR))) == idxInfected(iZR)
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);
                            TS_ZR(iNZR_S) = min([TSrand  APLS(APLS(:,1)==idxInfected(iZR),3)]); 
                        elseif SPLV(ismember(SPLV(:,1),idxInfected(iZR))) == idxInfected(iZR)
                            pv = rand();
                            if pv > eS
                                iNZR_V = iNZR_V+1; 
                                idexV_R(iNZR_V) = idxInfected(iZR);
                                TV_ZR(iNZR_V) = min([TSrand  APLV(APLV(:,1)==idxInfected(iZR),3)]);
                            end
                        end
                    end

                    % If I => Q
                    for iZQ = 1:Z
                        TSrand = randpdf(inf_profile_UV(:,2),inf_profile_UV(:,1),[1,1]);
                        if SS(ismember(SS(:,1),idxInfected(iZQ))) == idxInfected(iZQ)
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S+1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = TSrand;
                            end
                        elseif SV(ismember(SV(:,1),idxInfected(iZQ))) == idxInfected(iZQ)
                            pv = rand();
                            if pv > eS
                                if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                    iNZQ_V = iNZQ_V+1;
                                    idexV_Q(iNZQ_V) = idxInfected(iZQ);
                                    TV_ZQ(iNZQ_V) = TSrand;
                                end
                            end
                        elseif SPLS(ismember(SPLS(:,1),idxInfected(iZQ))) == idxInfected(iZQ)
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S+1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = min([TSrand APLS(APLS(:,1)==idxInfected(iZQ),3)]);
                            end
                        elseif SPLV(ismember(SPLV(:,1),idxInfected(iZQ))) == idxInfected(iZQ)
                            pv = rand();
                            if pv > eS
                                if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                    iNZQ_V = iNZQ_V+1;
                                    idexV_Q(iNZQ_V) = idxInfected(iZQ);
                                    TV_ZQ(iNZQ_V) = min([TSrand APLV(APLV(:,1)==idxInfected(iZQ),3)]);
                                end
                            end   
                        end
                    end

                    DetectionTime_all(idPLS(i))  = DetectionTime;

                    TSV_ZR = [TS_ZR TV_ZR];

                    TSV_ZQ = [TS_ZQ TV_ZQ];
                    idexSV_Q = [idexS_Q idexV_Q];

                    if cc_par==1
                        idxCC_TestedandInfected = idxCC_toTest(ismember(idxCC_toTest,idexSV_Q));
                        for ii1=1:length(idxCC_TestedandInfected) 
                            Testing_time_temp = inf;
                            if DetectionTime ~= inf
                                Testing_time_temp = round(DetectionTime-TSV_ZQ(ii1)+close_contact_test_period-mod(close_contact_test_period,close_contact_test_freq),2);
                            end
                            Testing_time(idxCC_TestedandInfected(ii1)) = min([Testing_time(idxCC_TestedandInfected(ii1)) Testing_time_temp]);
                        end
                        CloseContact_data = [CloseContact_data;idPLS(i) {idxCC} {idexSV_Q} {TSV_ZQ} {idxCC_toTest} {idxCC_TestedandInfected}];
                    end

                    % Secondary infection
                    if cc_par==1
                        if DetectionTime == inf        % not quarantined
                            SPLS(i,5) = 5;
                            Z2 = length(TSV_ZR); 
                        elseif (ismember(idPLS(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                        (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) ||...
                        (isnan(Testing_time(idPLS(i)))==0 && DetectionTime==PosCCTestTime+turnAround_time_cc)  % QT
                            SPLS(i,5) = 7;
                            Z2 = length(TSV_ZQ); 
                        else % QS
                            SPLS(i,5) = 6;
                            Z2 = length(TSV_ZQ); 
                        end
                    else
                        if DetectionTime == inf        % not quarantined
                            SPLS(i,5) = 5;
                            Z2 = length(TSV_ZR);  
                        elseif (ismember(idPLS(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                        (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) % QT
                            SPLS(i,5) = 7;
                            Z2 = length(TSV_ZQ); 
                        else % QS
                            SPLS(i,5) = 6;
                            Z2 = length(TSV_ZQ); 
                        end
                    end

                    SPLS(i,8) = Z2;

                    if (SPLS(i,5) == 5) % R
                        for id_run=1:length(idexS_R)
                            if length(find(SS(:,1)==idexS_R(id_run))) == 1                                
                                SS(SS(:,1)==idexS_R(id_run),2) = 10;
                                SS(SS(:,1)==idexS_R(id_run),3) = 2;
                                SS(SS(:,1)==idexS_R(id_run),6) = SPLS(i,1);
                                AS(AS(:,1)==idexS_R(id_run),2) = 0; 
                                AS(AS(:,1)==idexS_R(id_run),3) = TS_ZR(id_run);
                            elseif length(find(SPLS(:,1)==idexS_R(id_run))) == 1                              
                                SPLS(SPLS(:,1)==idexS_R(id_run),2) = 10;
                                SPLS(SPLS(:,1)==idexS_R(id_run),3) = 2;
                                SPLS(SPLS(:,1)==idexS_R(id_run),6) = SPLS(i,1);
                                APLS(APLS(:,1)==idexS_R(id_run),2) = 0; 
                                APLS(APLS(:,1)==idexS_R(id_run),3) = TS_ZR(id_run);
                            end
                        end

                        for id_run=1:length(idexV_R)
                            if length(find(SV(:,1)==idexV_R(id_run))) == 1                                
                                SV(SV(:,1)==idexV_R(id_run),2) = 10;
                                SV(SV(:,1)==idexV_R(id_run),3) = 12;
                                SV(SV(:,1)==idexV_R(id_run),6) = SPLS(i,1);
                                AV(AV(:,1)==idexV_R(id_run),2) = 0; 
                                AV(AV(:,1)==idexV_R(id_run),3) = TV_ZR(id_run);
                            elseif length(find(SPLV(:,1)==idexV_R(id_run))) == 1                              
                                SPLV(SPLV(:,1)==idexV_R(id_run),2) = 10;
                                SPLV(SPLV(:,1)==idexV_R(id_run),3) = 2;
                                SPLV(SPLV(:,1)==idexV_R(id_run),6) = SPLS(i,1);
                                APLV(APLV(:,1)==idexV_R(id_run),2) = 0; 
                                APLV(APLV(:,1)==idexV_R(id_run),3) = TV_ZR(id_run);
                            end
                        end

                        if isempty(TSV_ZR)
                            tL = 0;
                            tI = 0;
                        else
                            TSV_ZR = sort(TSV_ZR);
                            tL = min(TSV_ZR);
                            tI = TSV_ZR(end)-TSV_ZR(1);
                        end

                        tRt = ceil(t+TSV_ZR);
                        tRt(tRt>Nsave) = [];
                        Rt(tRt,2) = Rt(tRt,2)+1;

                        tRt2 = ceil((TSV_ZR+t)/dt);
                        tRt2(tRt2>Nstep) = [];
                        Rt2(tRt2,2) = Rt2(tRt2,2)+1; 

                    else    % Q
                        for id_run=1:length(idexS_Q)
                            if length(find(SS(:,1)==idexS_Q(id_run))) == 1                                
                                SS(SS(:,1)==idexS_Q(id_run),2) = 10;
                                SS(SS(:,1)==idexS_Q(id_run),3) = 2;
                                SS(SS(:,1)==idexS_Q(id_run),6) = SPLS(i,1);
                                AS(AS(:,1)==idexS_Q(id_run),2) = 0; 
                                AS(AS(:,1)==idexS_Q(id_run),3) = TS_ZQ(id_run);
                            elseif length(find(SPLS(:,1)==idexS_Q(id_run))) == 1                              
                                SPLS(SPLS(:,1)==idexS_Q(id_run),2) = 10;
                                SPLS(SPLS(:,1)==idexS_Q(id_run),3) = 2;
                                SPLS(SPLS(:,1)==idexS_Q(id_run),6) = SPLS(i,1);
                                APLS(APLS(:,1)==idexS_Q(id_run),2) = 0; 
                                APLS(APLS(:,1)==idexS_Q(id_run),3) = TS_ZQ(id_run);
                            end
                        end

                        for id_run=1:length(idexV_Q)
                            if length(find(SV(:,1)==idexV_Q(id_run))) == 1                                
                                SV(SV(:,1)==idexV_Q(id_run),2) = 10;
                                SV(SV(:,1)==idexV_Q(id_run),3) = 12;
                                SV(SV(:,1)==idexV_Q(id_run),6) = SPLS(i,1);
                                AV(AV(:,1)==idexV_Q(id_run),2) = 0; 
                                AV(AV(:,1)==idexV_Q(id_run),3) = TV_ZQ(id_run);
                            elseif length(find(SPLV(:,1)==idexV_Q(id_run))) == 1                              
                                SPLV(SPLV(:,1)==idexV_Q(id_run),2) = 10;
                                SPLV(SPLV(:,1)==idexV_Q(id_run),3) = 2;
                                SPLV(SPLV(:,1)==idexV_Q(id_run),6) = SPLS(i,1);
                                APLV(APLV(:,1)==idexV_Q(id_run),2) = 0; 
                                APLV(APLV(:,1)==idexV_Q(id_run),3) = TV_ZQ(id_run);
                            end
                        end

                        if isempty(TSV_ZQ) 
                            tL = 0;
                            tI = 0;
                        else
                            TSV_ZQ = sort(TSV_ZQ);
                            tL = min(TSV_ZQ);
                            tI = TSV_ZQ(end)-TSV_ZQ(1);
                        end

                        tRt = ceil(t+TSV_ZQ);
                        tRt(tRt>Nsave) = [];
                        Rt(tRt,2) = Rt(tRt,2)+1;

                        tRt2 = ceil((TSV_ZQ+t)/dt);
                        tRt2(tRt2>Nstep) = [];
                        Rt2(tRt2,2) = Rt2(tRt2,2)+1; 

                    end

                    APLS(i,4) = tL;
                    APLS(i,5)  = tI;
               end
           end 

           %NPLV
           LV_count = 0;
           for i = 1:NPLV
               if APLV(i,2) > APLV(i,3)
                   LV_count = LV_count+1;
                   SPLV(i,2) = SPLV(i,3);     % Update state
                   tc = 5.1;
                   APLV(i,2) = APLV(i,2) - APLV(i,3);

                   % Next state
                    SorA = rand();
                    if SorA <= k_sym_V     % IS
                        SPLV(i,4) = 13;
                        R0 = R00*(1-eI);
                    else                % IA
                        SPLV(i,4) = 14;
                        R0 = R00*k_reduce*(1-eI);
                    end

                    SPLV(i,3) = SPLV(i,4); % Update next state

                    p = (1+R0/k)^-1 ;
                    Z = nbinrnd(k,p);
                    if cc_par==1
                        idxCC = CC_matrix{idPLV(i)};
                        idxCC_toTest = CC_matrix_toTest{idPLV(i)};
                        if SPLV(i,4) == 13 % IS
                            Z = Z_all(idPLV(i),3);             
                        elseif SPLV(i,4) == 14 % IA
                            Z = Z_all(idPLV(i),4);   
                        end
                    end

                    PosRouTestTime  = inf;
                    RanTestTime     = inf;  
                    PosCCTestTime   = inf;
    %                 IsoP            = 0;

                    DetectionTime = inf;
            
                    VL_since_iso            = [];
                    num_test_till_first_neg = [];

                    if SPLV(i,4) == 13 % IS
    %                     DetectionTime = tc; 
                        if rand() < symptom_test_pop
                            DetectionTime = tc + turnAround_time_symp; 
            
                            VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                            num_test_till_first_neg = min(find(VL_since_iso < LOD_symp))-1;  % number of test untill first negative result
              
                        end          
                    end 
                    %---------routine testing--------
                    if ismember(idPLV(i),idxRouT) == 1  
                        PosRouTestTime = round(routine_test_freq*rand(1) + min(V_all_pos_rou),2) + turnAround_time_rou;
                        DetectionTime  = min([DetectionTime PosRouTestTime]);
            
                        VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                        num_test_till_first_neg = min(find(VL_since_iso < LOD_rou))-1;  % number of test untill first negative result
            
                    end
                    %---------------------------------

                    %----------random testing---------
                    RanTT = rand(size(V_all,1),1) < random_test_capa/N0*dt; 
                    RanTT_dt = find(RanTT==1);
                    for RanTT_ii=1:length(RanTT_dt)
                        if V_all(RanTT_dt(RanTT_ii),2) > LOD_ran
                            RanTestTime = V_all(RanTT_dt(RanTT_ii),1) + turnAround_time_ran;
                            DetectionTime = min([DetectionTime RanTestTime]);
            
                            VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                            num_test_till_first_neg = min(find(VL_since_iso < LOD_ran))-1;  % number of test untill first negative result
            
                        end
                    end
                    %----------------------------------

                    %------------cc testing------------
                    if cc_par==1
                        if (Testing_time(idPLV(i)) > 0) && (Testing_time(idPLV(i))~=inf)
                            if V_all(V_all(:,1)==Testing_time(idPLV(i)),2) > LOD_cc
                                LastCCTestTime = Testing_time(idPLV(i));
                                while V_all(V_all(:,1) == round(LastCCTestTime-close_contact_test_freq,2),2) > LOD_cc
                                    LastCCTestTime = LastCCTestTime - close_contact_test_freq;
                                end
                                PosCCTestTime = LastCCTestTime + turnAround_time_cc;
                                DetectionTime = min([DetectionTime PosCCTestTime]);
            
                                VL_since_iso = V_all(find(ismember(V_all(:,1),round([DetectionTime:isolation_test_freq:V_all(end,1)],2))),2);
                                num_test_till_first_neg = min(find(VL_since_iso < LOD_cc))-1;  % number of test untill first negative result
            
                            end
                        end
                    end
                    %----------------------------------

                    %-------------isolation------------
                    IsoP = num_test_till_first_neg*isolation_test_freq;
                    date_out_Q = DetectionTime + IsoP;
                    %----------------------------------

                    Iso_data = [Iso_data;idPLV(i) DetectionTime~=inf ...
                        DetectionTime {date_out_Q} {IsoP} {num_test_till_first_neg}];

                    Testing_data = [Testing_data; idPLV(i) ismember(idPLV(i),idxRouT) PosRouTestTime...
                           ~isempty(find(RanTT==1)) RanTestTime ...
                           PosCCTestTime DetectionTime];

                    SPLV(i,7) = Z; 

                    % Making infection
                    if cc_par==1
                        idxInfected = idxCC(randperm(length(idxCC),Z)); 
                    else
                        idxInfected = randperm(N,Z);   % infectee
                    end

                    TS_ZQ = [];     TS_ZR = [];     % ts for each infectee
                    TV_ZQ = [];     TV_ZR = [];
                    iNZQ_S = 0;     iNZR_S = 0;       % N successful infection
                    iNZQ_V = 0;     iNZR_V = 0;
                    idexS_Q = [];   idexS_R = [];   % index of infectee (successful infection)
                    idexV_Q = [];   idexV_R = []; 

                    % If I => R                    
                    for iZR = 1:Z
                        TSrand = randpdf(inf_profile_V(:,2),inf_profile_V(:,1),[1,1]);
                        if SS(ismember(SS(:,1),idxInfected(iZR))) == idxInfected(iZR)
                            iNZR_S = iNZR_S+1;
                            idexS_R(iNZR_S) = idxInfected(iZR);
                            TS_ZR(iNZR_S) = TSrand;
                        elseif SV(ismember(SV(:,1),idxInfected(iZR))) == idxInfected(iZR)
                            pv = rand();
                            if pv > eS
                                iNZR_V = iNZR_V+1; 
                                idexV_R(iNZR_V) = idxInfected(iZR);        
                                TV_ZR(iNZR_V) = TSrand;
                            end
                        elseif SPLS(ismember(SPLS(:,1),idxInfected(iZR))) == idxInfected(iZR)
                            iNZR_S = iNZR_S+1; 
                            idexS_R(iNZR_S) = idxInfected(iZR);
                            TS_ZR(iNZR_S) = min([TSrand  APLS(APLS(:,1)==idxInfected(iZR),3)]); 
                        elseif SPLV(ismember(SPLV(:,1),idxInfected(iZR))) == idxInfected(iZR)
                            pv = rand();
                            if pv > eS
                                iNZR_V = iNZR_V+1; 
                                idexV_R(iNZR_V) = idxInfected(iZR);
                                TV_ZR(iNZR_V) = min([TSrand  APLV(APLV(:,1)==idxInfected(iZR),3)]);
                            end
                        end
                    end

                    % If I => Q
                    for iZQ = 1:Z
                        TSrand = randpdf(inf_profile_V(:,2),inf_profile_V(:,1),[1,1]);
                        if SS(ismember(SS(:,1),idxInfected(iZQ))) == idxInfected(iZQ)
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S+1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = TSrand;
                            end
                        elseif SV(ismember(SV(:,1),idxInfected(iZQ))) == idxInfected(iZQ)
                            pv = rand();
                            if pv > eS
                                if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                    iNZQ_V = iNZQ_V+1;
                                    idexV_Q(iNZQ_V) = idxInfected(iZQ);
                                    TV_ZQ(iNZQ_V) = TSrand;
                                end
                            end
                        elseif SPLS(ismember(SPLS(:,1),idxInfected(iZQ))) == idxInfected(iZQ)
                            if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                iNZQ_S = iNZQ_S+1;
                                idexS_Q(iNZQ_S) = idxInfected(iZQ);
                                TS_ZQ(iNZQ_S) = min([TSrand APLS(APLS(:,1)==idxInfected(iZQ),3)]);
                            end
                        elseif SPLV(ismember(SPLV(:,1),idxInfected(iZQ))) == idxInfected(iZQ)
                            pv = rand();
                            if pv > eS
                                if TSrand < DetectionTime || TSrand > (DetectionTime+IsoP)
                                    iNZQ_V = iNZQ_V+1;
                                    idexV_Q(iNZQ_V) = idxInfected(iZQ);
                                    TV_ZQ(iNZQ_V) = min([TSrand APLV(APLV(:,1)==idxInfected(iZQ),3)]);
                                end
                            end   
                        end
                    end

                    DetectionTime_all(idPLV(i))  = DetectionTime;

                    TSV_ZR = [TS_ZR TV_ZR];

                    TSV_ZQ = [TS_ZQ TV_ZQ];
                    idexSV_Q = [idexS_Q idexV_Q];

                    if cc_par==1
                        idxCC_TestedandInfected = idxCC_toTest(ismember(idxCC_toTest,idexSV_Q));
                        for ii1=1:length(idxCC_TestedandInfected) 
                            Testing_time_temp = inf;
                            if DetectionTime ~= inf
                                Testing_time_temp = round(DetectionTime-TSV_ZQ(ii1)+close_contact_test_period-mod(close_contact_test_period,close_contact_test_freq),2);
                            end
                            Testing_time(idxCC_TestedandInfected(ii1)) = min([Testing_time(idxCC_TestedandInfected(ii1)) Testing_time_temp]);
                        end
                        CloseContact_data = [CloseContact_data;idPLV(i) {idxCC} {idexSV_Q} {TSV_ZQ} {idxCC_toTest} {idxCC_TestedandInfected}];
                    end

                    % Secondary infection   
                    if cc_par==1
                        if DetectionTime == inf        % not quarantined
                            SPLV(i,5) = 15;
                            Z2 = length(TSV_ZR);  
                        elseif (ismember(idPLV(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                        (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) ||...
                        (isnan(Testing_time(idPLV(i)))==0 && DetectionTime==PosCCTestTime+turnAround_time_cc)  % QT
                            SPLV(i,5) = 17;
                            Z2 = length(TSV_ZQ); 
                        else % QS
                            SPLV(i,5) = 16;
                            Z2 = length(TSV_ZQ); 
                        end
                    else
                        if DetectionTime == inf        % not quarantined
                            SPLV(i,5) = 15;
                            Z2 = length(TSV_ZR);  
                        elseif (ismember(idPLV(i),idxRouT)==1 && DetectionTime==PosRouTestTime+turnAround_time_rou) ||...
                        (isempty(find(RanTT==1,1))==0 && DetectionTime==RanTestTime+turnAround_time_ran) % QT
                            SPLV(i,5) = 17;
                            Z2 = length(TSV_ZQ); 
                        else % QS
                            SPLV(i,5) = 16;
                            Z2 = length(TSV_ZQ); 
                        end
                    end

                    SPLV(i,8) = Z2;

                    if (SPLV(i,5) == 15) % R
                        for id_run=1:length(idexS_R)
                            if length(find(SS(:,1)==idexS_R(id_run))) == 1                                
                                SS(SS(:,1)==idexS_R(id_run),2) = 10;
                                SS(SS(:,1)==idexS_R(id_run),3) = 2;
                                SS(SS(:,1)==idexS_R(id_run),6) = SPLV(i,1);
                                AS(AS(:,1)==idexS_R(id_run),2) = 0; 
                                AS(AS(:,1)==idexS_R(id_run),3) = TS_ZR(id_run);
                            elseif length(find(SPLS(:,1)==idexS_R(id_run))) == 1                              
                                SPLS(SPLS(:,1)==idexS_R(id_run),2) = 10;
                                SPLS(SPLS(:,1)==idexS_R(id_run),3) = 2;
                                SPLS(SPLS(:,1)==idexS_R(id_run),6) = SPLV(i,1);
                                APLS(APLS(:,1)==idexS_R(id_run),2) = 0; 
                                APLS(APLS(:,1)==idexS_R(id_run),3) = TS_ZR(id_run);
                            end
                        end

                        for id_run=1:length(idexV_R)
                            if length(find(SV(:,1)==idexV_R(id_run))) == 1                                
                                SV(SV(:,1)==idexV_R(id_run),2) = 10;
                                SV(SV(:,1)==idexV_R(id_run),3) = 12;
                                SV(SV(:,1)==idexV_R(id_run),6) = SPLV(i,1);
                                AV(AV(:,1)==idexV_R(id_run),2) = 0; 
                                AV(AV(:,1)==idexV_R(id_run),3) = TV_ZR(id_run);
                            elseif length(find(SPLV(:,1)==idexV_R(id_run))) == 1                              
                                SPLV(SPLV(:,1)==idexV_R(id_run),2) = 10;
                                SPLV(SPLV(:,1)==idexV_R(id_run),3) = 2;
                                SPLV(SPLV(:,1)==idexV_R(id_run),6) = SPLV(i,1);
                                APLV(APLV(:,1)==idexV_R(id_run),2) = 0; 
                                APLV(APLV(:,1)==idexV_R(id_run),3) = TV_ZR(id_run);
                            end
                        end

                        if isempty(TSV_ZR)
                            tL = 0;
                            tI = 0;
                        else
                            TSV_ZR = sort(TSV_ZR);
                            tL = min(TSV_ZR);
                            tI = TSV_ZR(end)-TSV_ZR(1);
                        end

                        tRt = ceil(t+TSV_ZR); 
                        tRt(tRt>Nsave) = [];
                        Rt(tRt,2) = Rt(tRt,2)+1;

                        tRt2 = ceil((TSV_ZR+t)/dt);
                        tRt2(tRt2>Nstep) = [];
                        Rt2(tRt2,2) = Rt2(tRt2,2)+1; 

                    else    % Q
                        for id_run=1:length(idexS_Q)
                            if length(find(SS(:,1)==idexS_Q(id_run))) == 1                                
                                SS(SS(:,1)==idexS_Q(id_run),2) = 10;
                                SS(SS(:,1)==idexS_Q(id_run),3) = 2;
                                SS(SS(:,1)==idexS_Q(id_run),6) = SPLV(i,1);
                                AS(AS(:,1)==idexS_Q(id_run),2) = 0; 
                                AS(AS(:,1)==idexS_Q(id_run),3) = TS_ZQ(id_run);
                            elseif length(find(SPLS(:,1)==idexS_Q(id_run))) == 1                              
                                SPLS(SPLS(:,1)==idexS_Q(id_run),2) = 10;
                                SPLS(SPLS(:,1)==idexS_Q(id_run),3) = 2;
                                SPLS(SPLS(:,1)==idexS_Q(id_run),6) = SPLV(i,1);
                                APLS(APLS(:,1)==idexS_Q(id_run),2) = 0; 
                                APLS(APLS(:,1)==idexS_Q(id_run),3) = TS_ZQ(id_run);
                            end
                        end

                        for id_run=1:length(idexV_Q)
                            if length(find(SV(:,1)==idexV_Q(id_run))) == 1                                
                                SV(SV(:,1)==idexV_Q(id_run),2) = 10;
                                SV(SV(:,1)==idexV_Q(id_run),3) = 12;
                                SV(SV(:,1)==idexV_Q(id_run),6) = SPLV(i,1);
                                AV(AV(:,1)==idexV_Q(id_run),2) = 0; 
                                AV(AV(:,1)==idexV_Q(id_run),3) = TV_ZQ(id_run);
                            elseif length(find(SPLV(:,1)==idexV_Q(id_run))) == 1                              
                                SPLV(SPLV(:,1)==idexV_Q(id_run),2) = 10;
                                SPLV(SPLV(:,1)==idexV_Q(id_run),3) = 2;
                                SPLV(SPLV(:,1)==idexV_Q(id_run),6) = SPLV(i,1);
                                APLV(APLV(:,1)==idexV_Q(id_run),2) = 0; 
                                APLV(APLV(:,1)==idexV_Q(id_run),3) = TV_ZQ(id_run);
                            end
                        end

                        if isempty(TSV_ZQ) 
                            tL = 0;
                            tI = 0;
                        else
                            TSV_ZQ = sort(TSV_ZQ);
                            tL = min(TSV_ZQ);
                            tI = TSV_ZQ(end)-TSV_ZQ(1);
                        end

                        tRt = ceil(t+TSV_ZQ);
                        tRt(tRt>Nsave) = [];
                        Rt(tRt,2) = Rt(tRt,2)+1;

                        tRt2 = ceil((TSV_ZQ+t)/dt);
                        tRt2(tRt2>Nstep) = [];
                        Rt2(tRt2,2) = Rt2(tRt2,2)+1; 

                    end
                    APLV(i,4) = tL;
                    APLV(i,5)  = tI;
               end
           end 

           % LS -> I
           for i = 1:NLS
               if ALS(i,2) >= ALS(i,4)
                   SLS(i,2) = SLS(i,3);   % Update state
                   if SLS(i,3) == 3
                       CumISS = CumISS+1;
                   elseif SLS(i,3) == 4
                       CumIAS = CumIAS+1;
                   end
                   SLS(i,3) = SLS(i,5);   % Update next state
                   ALS(i,2) = ALS(i,2) - ALS(i,4);     
               end 
           end

           % LV -> I
           for i = 1:NLV
               if ALV(i,2) >= ALV(i,4)
                   SLV(i,2) = SLV(i,3);   % Update state
                   if SLV(i,3) == 13
                       CumISV = CumISV+1;
                   elseif SLV(i,3) == 14
                       CumIAV = CumIAV+1;
                   end
                   SLV(i,3) = SLV(i,5);   % Update next state
                   ALV(i,2) = ALV(i,2) - ALV(i,4);     
               end 
           end

           % ISS -> QTS/QSS/RS
           for i =1:NISS
               if AISS(i,2) >= AISS(i,5)
                   SISS(i,2) = SISS(i,3);
                   if SISS(i,2) == 5
                       CumRS = CumRS + 1;
                   elseif SISS(i,2) == 6
                       CumQSS = CumQSS + 1;
                   else
                       CumQTS = CumQTS + 1;
                   end
               SISS(i,3) = SISS(i,5);    
               AISS(i,2) = AISS(i,2) - AISS(i,5);          
               end
           end

           % ISV -> QTV/QSV/RV
           for i =1:NISV
               if AISV(i,2) >= AISV(i,5)
                   SISV(i,2) = SISV(i,3);
                   if SISV(i,2) == 15
                       CumRV = CumRV + 1;
                   elseif SISV(i,2) == 16
                       CumQSV = CumQSV + 1;
                   else 
                       CumQTV = CumQTV + 1;
                   end
               SISV(i,3) = SISV(i,5); 
               AISV(i,2) = AISV(i,2) - AISV(i,5);          
               end
           end

           % IAS -> QTS/QSS/RS
           for i = 1:NIAS
               if AIAS(i,2) >= AIAS(i,5)
                   SIAS(i,2) = SIAS(i,3);
                   if SIAS(i,2) == 5
                       CumRS = CumRS + 1;
                   elseif SIAS(i,2) == 6
                       CumQSS = CumQSS + 1;
                   else
                       CumQTS = CumQTS + 1;
                   end
               SIAS(i,3) = SIAS(i,5);  
               AIAS(i,2) = AIAS(i,2) - AIAS(i,5);          
               end
           end

           % IAV -> QTV/QSV/RV
           for i = 1:NIAV
               if AIAV(i,2) >= AIAV(i,5)
                   SIAV(i,2) = SIAV(i,3);
                   if SIAV(i,2) == 15
                       CumRV = CumRV + 1;
                   elseif SIAV(i,2) == 16
                       CumQSV = CumQSV + 1;
                   else 
                       CumQTV = CumQTV + 1;
                   end
               SIAV(i,3) = SIAV(i,5); 
               AIAV(i,2) = AIAV(i,2) - AIAV(i,5);          
               end
           end

%            State = [SS; SLS; SISS; SIAS; SV; SLV; SISV; SIAV; SRS; SRV; SQSS; SQTS; SQSV; SQTV; SPLS; SPLV];
%            Alarm = [AS; ALS; AISS; AIAS; AV; ALV; AISV; AIAV; ARS; ARV; AQSS; AQTS; AQSV; AQTV; APLS; APLV];
           %----------------------------------------------------------------------------------%
           State(SS(:,1),2:end)     = SS(:,2:end);      State(SLS(:,1),2:end)   = SLS(:,2:end); 
           State(SISS(:,1),2:end)   = SISS(:,2:end);    State(SIAS(:,1),2:end)  = SIAS(:,2:end); 
           State(SV(:,1),2:end)     = SV(:,2:end);      State(SLV(:,1),2:end)   = SLV(:,2:end);
           State(SISV(:,1),2:end)   = SISV(:,2:end);    State(SIAV(:,1),2:end)  = SIAV(:,2:end); 
           State(SRS(:,1),2:end)    = SRS(:,2:end);     State(SRV(:,1),2:end)   = SRV(:,2:end); 
           State(SQSS(:,1),2:end)   = SQSS(:,2:end);    State(SQTS(:,1),2:end)  = SQTS(:,2:end);
           State(SQSV(:,1),2:end)   = SQSV(:,2:end);    State(SQTV(:,1),2:end)  = SQTV(:,2:end); 
           State(SPLS(:,1),2:end)   = SPLS(:,2:end);    State(SPLV(:,1),2:end)  = SPLV(:,2:end);
           
           Alarm(AS(:,1),2:end)     = AS(:,2:end);      Alarm(ALS(:,1),2:end)   = ALS(:,2:end); 
           Alarm(AISS(:,1),2:end)   = AISS(:,2:end);    Alarm(AIAS(:,1),2:end)  = AIAS(:,2:end); 
           Alarm(AV(:,1),2:end)     = AV(:,2:end);      Alarm(ALV(:,1),2:end)   = ALV(:,2:end);
           Alarm(AISV(:,1),2:end)   = AISV(:,2:end);    Alarm(AIAV(:,1),2:end)  = AIAV(:,2:end); 
           Alarm(ARS(:,1),2:end)    = ARS(:,2:end);     Alarm(ARV(:,1),2:end)   = ARV(:,2:end); 
           Alarm(AQSS(:,1),2:end)   = AQSS(:,2:end);    Alarm(AQTS(:,1),2:end)  = AQTS(:,2:end);
           Alarm(AQSV(:,1),2:end)   = AQSV(:,2:end);    Alarm(AQTV(:,1),2:end)  = AQTV(:,2:end); 
           Alarm(APLS(:,1),2:end)   = APLS(:,2:end);    Alarm(APLV(:,1),2:end)  = APLV(:,2:end);
           %-----------------------------------------------------------------------------------%
           
           SS = []; SLS = []; SISS = []; SIAS = []; SV = []; SLV = []; SISV = []; SIAV =[]; SRS = []; SRV = []; SQSS=[]; SQTS=[]; SQSV=[]; SQTV=[]; SPLS = []; SPLV = [];
           AS = []; ALS = []; AISS = []; AIAS = []; AV = []; ALV = []; AISV = []; AIAV =[]; ARS = []; ARV = []; AQSS=[]; AQTS=[]; AQSV=[]; AQTV=[]; APLS = []; APLV = [];

           Lall_count = LS_count+LV_count;
           LS_count_save(istep,irepeat) = LS_count;
           LV_count_save(istep,irepeat) = LV_count;
           Lall_count_save(istep,irepeat) = Lall_count;

           %% Update time    
           t = t + dt;
           Alarm(:,2) = Alarm(:,2) + dt;  

        end

        Z_save(:,irepeat) = State(:,7);
        Z2_save(:,irepeat) = State(:,8); 
        State_save(:,irepeat) = State(:,2); 

        meanDetectionTime(irepeat)      = mean(DetectionTime_all(DetectionTime_all~=inf));

        Testing_time(isnan(Testing_time)) = [];
        meanTesting_time_save(irepeat)  = mean(Testing_time(Testing_time~=Inf));
        Testing_time_save(irepeat,:)    = [length(find(Testing_time == inf)) length(find(Testing_time(Testing_time~=Inf) < 0)) length(find(Testing_time(Testing_time~=Inf) >= 0))];
        
        CCall_save=0; CCtest_save=0; CCinfectedtested_save=0;
        
        if cc_par==1
            for cc_temp3=1:length(CloseContact_data(:,2))
                CCall_save              = CCall_save+length(CloseContact_data{cc_temp3,2});
                CCtest_save             = CCtest_save+length(CloseContact_data{cc_temp3,5});
                CCinfectedtested_save   = CCinfectedtested_save+length(CloseContact_data{cc_temp3,6});
            end
            CCtest_matrix_save(irepeat,:) = [CCall_save CCtest_save CCinfectedtested_save];
        end
        
        Test_count(irepeat) = sum(Testing_data(:,7) ~= inf);
        
    end

    %% Time to elimination
    Echeck = PLSsave +PLVsave + LSsave + LVsave + ISSsave + IASsave + ISVsave + IAVsave;
    CumCheck = Rallsave + Qallsave;

    [~,ind_ex] = min(Echeck,[],1);

    %% New and cumulative cases
    CumIall = CumISsave+CumIAsave;
    Newcase = zeros(Nsave,Nrepeat);
    NewX = zeros(Nsave,Nrepeat);
    Newcase(2:Nsave,:) = CumIall(2:Nsave,:) - CumIall(1:(Nsave-1),:);
    NewX(2:Nsave,:) = CumXsave(2:Nsave,:) - CumXsave(1:(Nsave-1),:);

    filename = ['cc_pop' num2str(N0) '_' num2str(close_contact_test_pc*100) 'ccf' num2str(close_contact_test_freq) 'p' num2str(close_contact_test_period) '_'...
        num2str(symptom_test_pop*100) 'qf' num2str(isolation_test_freq) '_' ...
        num2str(routine_test_pc*100) 'rouf' num2str(routine_test_freq) '_' ...
        num2str(random_test_capa/N0*100) 'ran_' ...
        'ini' num2str(NPLS0+NPLV0+NLS0+NLV0+NISS0+NIAS0+NISV0+NIAV0) '_LOD' num2str(LOD_atk) '_ttsymp' num2str(turnAround_time_symp) ...
         '_ttrou' num2str(turnAround_time_rou) '_R0' num2str(R00*100) '_ttcc' num2str(turnAround_time_cc) ...
         '_ic' num2str(num_import_case) '_nrs0' num2str(NRS0) '_fix_lamp_batch' num2str(iBatch)];

    save(filename);

    close all
end