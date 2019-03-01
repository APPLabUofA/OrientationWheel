
% Code is based on loading saved ERP data

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
%% Select electrodes for analysis

% List electrodes to get ERPs
% elect_erp = [2 3 4 5 6 7 8 9 10 15 16 19 20 21 22 25 26 29 30];
% el_erp_names = {'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P3';'P4';'CP1';'CP2';'C3';'C4';'FC1';'FC2';'F3';'F4'};

elect_erp = [2 3 4 5 6 9 10 19 20 25 26];
el_erp_names = {'Oz';'Pz';'Cz';'FCz';'Fz';'PO3';'PO4';'CP1';'CP2';'FC1';'FC2'};

% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% ------------------ Load Previously Saved Data --------------------------
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Load data from Anal_ERP_Orient_byTarget.m
load('erp_out_byTarget.mat')

% Load data from Anal_ERP_Orient_byCatchTrial.m
erp_catchtrials = cell2mat(struct2cell(load('erp_out_catch_trials.mat')));

% load behavior data
load('ALLEEG_filt_byTargets_v3.mat');

% load settings
load('filt_byTargets_v3_Settings.mat');



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% ::::::::::::::::::  Plot the ERPs by electrode  ::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% ------------------------------------------------------------------------- 
% Average across subjects by errors
erp_out_byerr(:,:,1) = squeeze(mean(erp_out_x(:,:,:),1)); %small errors
erp_out_byerr(:,:,2) = squeeze(mean(erp_out_n(:,:,:),1)); %large errors
% erp_out_byerr(:,:,3) = squeeze(mean((erp_out_x(:,:,:)-erp_out_n(:,:,:)),1)); %difference
erp_out_byerr(:,:,3) = squeeze(mean(erp_catchtrials(:,:,:),1)); %catch trials (no target)
% Average across subjects - catch trials
erp_out_bycatch(:,:) = squeeze(mean(erp_catchtrials(:,:,:),1));
% ------------------------------------------------------------------------- 

% Plot ERPs with error bars
for ii = 1:length(elect_erp) 
    i_elect = elect_erp(ii); %for doing only a selection of electrodes
%     i_elect = ii; %selection done when making ERPs
    % get axes limits
    ymin = -6; ymax = 10;
    xmin = -200; xmax = 800;
    
    figure('Color',[1 1 1]); 
    boundedline(EEG.times,erp_out_bycatch(i_elect,:),squeeze(std(erp_catchtrials(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'g',...
            EEG.times,erp_out_byerr(i_elect,:,1),squeeze(std(erp_out_x(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'c',...
            EEG.times,erp_out_byerr(i_elect,:,2),squeeze(std(erp_out_n(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'m');
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    set(gca,'ydir','reverse'); xlim([xmin xmax]); ylim([ymin ymax])
    
    title([el_erp_names{ii} ': ERPs by Response Error']); 
    xlabel('Time (ms)'); ylabel('Voltage (uV)'); 
    xticks(xmin:100:xmax); yticks(ymin:2:ymax)    
        
    legend({'No Target','Accurate','Guess'},'Location','best');
    
%     savefig(['M:\Personal_Folders\Sarah\Manuscripts\Orientation_Wheel\Figures\ERP_' el_erp_names{ii}])
end
clear ii xmax xmin    


% /////////////////////////////////////////////////////////////////////////
%% ::::::::::::::::::::  Statistics: t-Test  ::::::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% This is for an approximation and should NOT be used in reports (need to
% use a parametric or non-parametric ANOVA & corrected for multiple
% comparisons)

%%% Pick your time window 350-550 for P3 %%%
% time1 = 350;
% time2 = 550;

time1 = 200;
time2 = 300;

% time1 = 150;
% time2 = 250;

% time1 = 0;
% time2 = 50;

% time1 = 500; %response screen at 567 ms
% time2 = 640;

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-1;

for ii = 1:length(elect_erp)    
    i_elect = elect_erp(ii); %for doing only a selection of electrodes

    % Accurate vs Guess
    [h,p,ci,stat] = ttest(mean(erp_out_x(:,i_elect,time_window),3),mean(erp_out_n(:,i_elect,time_window),3),.05,'both',1);
    ttest_elect(ii,1).AvG = h;
    ttest_elect(ii,2).AvG = p;
%     ttest_elect_ci(ii,1).AvG = ci(1);
%     ttest_elect_ci(ii,2).AvG = ci(2);
%     ttest_elect_stats(ii,1).AvG = stat;
    ttest_elect_id{ii,1} = el_erp_names{ii};
    clear h p ci stat
    
    % Accurate vs Catch
    [h,p,ci,stat] = ttest(mean(erp_out_x(:,i_elect,time_window),3),mean(erp_catchtrials(:,i_elect,time_window),3),.05,'both',1);
    ttest_elect(ii,1).AvC = h;
    ttest_elect(ii,2).AvC = p;
%     ttest_elect_ci(ii,1).AvC = ci(1);
%     ttest_elect_ci(ii,2).AvC = ci(2);
%     ttest_elect_stats(ii,1).AvC = stat;
    ttest_elect_id{ii,1} = el_erp_names{ii};
    clear h p ci stat
    
    % Guesses vs Catch
    [h,p,ci,stat] = ttest(mean(erp_out_n(:,i_elect,time_window),3),mean(erp_catchtrials(:,i_elect,time_window),3),.05,'both',1);
    ttest_elect(ii,1).GvC = h;
    ttest_elect(ii,2).GvC = p;
%     ttest_elect_ci(ii,1).GvC = ci(1);
%     ttest_elect_ci(ii,2).GvC = ci(2);
%     ttest_elect_stats(ii,1).GvC = stat;
    ttest_elect_id{ii,1} = el_erp_names{ii};
    clear i_elect h p ci stat
end
clear ii


clear time_window time1 time2




% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     '''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};


% Set the range of time to consider
% tWin{1} = [50 150];
% tWin{2} = [150 250];
% tWin{3} = [250 350];
% tWin{4} = [350 450];
% tWin{5} = [450 550];

% tWin{1} = [100 200];
% tWin{2} = [200 300];
% tWin{3} = [300 400];
% tWin{4} = [400 500];
% tWin{5} = [500 600];

tWin{1} = [0 200];
tWin{2} = [200 400];
tWin{3} = [400 600];

CLims1 = [-8 8]; %range in microvolts
nconds = 3; %number of plots
conds = {'Accurate';'Guess';'No Target'}; %labels for plots

for tw_i = 1:length(tWin) %loop through several time windows 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %this code finds the times you want from the timess variable
    time_window = find(EEG.times>= itWin(1),1):find(EEG.times>= itWin(2),1)-1;
    
    figure('Color',[1 1 1],'Position',[1 1 941 349]);

    for i_cond = 1:nconds %loop through conditions to make plot of each        
        subtightplot(1,3,i_cond,[0.02,0.02],[0.05,0.07],[0.05,0.05]);
        set(gca,'Color',[1 1 1]);
        
        temp = mean(erp_out_byerr(:,time_window,i_cond),2)'; %ERP within time window
        temp(1) = NaN; %so M2 is not included
        
        if i_cond == 4 %for making topography from conditon differences
            CLims = [-4 4]; %need smaller scale
            topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
                'plotchans',elect_erp,'emarker',{'.','k',11,1})
    %         topoplot(temp,EEG.chanlocs,'whitebk','on',0.6,'maplimits',...
    %             'plotchans',elect_erp,'emarker',{'.','k',11,1})
        else
            topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims1,...
            'plotchans',elect_erp,'emarker',{'.','k',11,1})
        end
        title(conds{i_cond});
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Voltage (uV)');
        clear temp
    end
    
    % Overall subplot title
    supertitle([num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'],...
        'FontSize',10.5)
    
    savefig(['M:\Personal_Folders\Sarah\Manuscripts\Orientation_Wheel\Figures\ERP_' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'])
    
    clear itWin time_window i_cond

end
clear tw_i nconds conds tWin CLims CLims1 t


% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Topographys     '''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};


% Set the range of time to consider
% tWin{1} = [50 150];
% tWin{2} = [150 250];
% tWin{3} = [250 350];
% tWin{4} = [350 450];
% tWin{5} = [450 550];

% tWin{1} = [100 200];
% tWin{2} = [200 300];
% tWin{3} = [300 400];
% tWin{4} = [400 500];
% tWin{5} = [500 600];

tWin{1} = [0 200];
tWin{2} = [200 400];
tWin{3} = [400 600];

CLims1 = [-8 8]; %range in microvolts
nconds = 3; %number of plots
conds = {'Accurate';'Guess';'No Target'}; %labels for plots

for tw_i = 1:length(tWin) %loop through several time windows 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %this code finds the times you want from the timess variable
    time_window = find(EEG.times>= itWin(1),1):find(EEG.times>= itWin(2),1)-1;
    
    figure('Color',[1 1 1],'Position',[1 1 941 349]);

    for i_cond = 1:nconds %loop through conditions to make plot of each        
        subtightplot(1,3,i_cond,[0.02,0.02],[0.05,0.07],[0.05,0.05]);
        set(gca,'Color',[1 1 1]);
        
        temp = mean(erp_out_byerr(:,time_window,i_cond),2)'; %ERP within time window
        temp(1) = NaN; %so M2 is not included
        
        if i_cond == 4 %for making topography from conditon differences
            CLims = [-4 4]; %need smaller scale
            topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
                'plotchans',elect_erp,'emarker',{'.','k',11,1})
    %         topoplot(temp,EEG.chanlocs,'whitebk','on',0.6,'maplimits',...
    %             'plotchans',elect_erp,'emarker',{'.','k',11,1})
        else
            topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims1,...
            'plotchans',elect_erp,'emarker',{'.','k',11,1})
        end
        title(conds{i_cond});
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Voltage (uV)');
        clear temp
    end
    
    % Overall subplot title
    supertitle([num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'],...
        'FontSize',10.5)
    
    savefig(['M:\Personal_Folders\Sarah\Manuscripts\Orientation_Wheel\Figures\ERP_' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'])
    
    clear itWin time_window i_cond

end
clear tw_i nconds conds tWin CLims CLims1 t

nconds = 3; %number of plots
% conds = {'Accurate';'Guess';'No Target'}; %labels for plots
conds = {'Accurate-Guesses';'Accurate-No Target';'Guesses-No Target'}; %labels for plots



% /////////////////////////////////////////////////////////////////////////
%% ::::::::::::::::::  Statistics: Correlation  :::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Put all the response errors across subjects into vector
resp_errdeg_cat = cat(2,resp_errdeg{1:end});
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%% Pick your time window 300-500 for P3 %%%
% time1 = 300;
% time2 = 550;

time1 = 200;
time2 = 300;

% time1 = 150;
% time2 = 250;

% time1 = 0;
% time2 = 50;

% time1 = 500; %response screen at 567 ms
% time2 = 640;

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-1;

amp_out_timewin = cell(length(exp.participants),length(elect_erp)); %pre-allocate
% get single trial amplitudes during time window
for i_part = 1:length(exp.participants) % --------------------------------- 
    % Calculate ERP
    for ii = 1:length(elect_erp)
        i_elect = elect_erp(ii); %for doing only a selection of electrodes
        amp_out_timewin{i_part,i_elect} = squeeze(mean(ALLEEG(i_part).data(i_elect,time_window,:),2));
    end
    clear ii i_elect
end
clear i_part


%% Run Permutation test

nperms = 10000; %number of permutations used to estimate null distribution
resp_errrad_cat = circ_ang2rad(resp_errdeg_cat); %convert to radians for circular corr
n_resp = length(resp_errrad_cat); %number of observations to permute

% re-set values
permtest.rho_obs = []; 
permtest.p_z = []; 
permtest.p_n = []; 

% Permutation test at each electrode
for ii = 1:5 %test central electrodes only
    i_elect = elect_erp(ii); %for doing only a selection of electrodes
    
    % Put all the amplitudes across subjects into vector
    tmp_amp = cat(1,amp_out_timewin{1:end,i_elect});
    
    % Make distribution of null-hypothesis test statistic
    rho_perm = zeros(1,nperms); %pre-allocate
    for i_perm = 1:nperms
        order_resp = randperm(n_resp); %randomly set order of data
        [rho(i_perm),pval] = circ_corrcl(resp_errrad_cat(order_resp),tmp_amp');
        clear order_resp pval
    end
    clear i_perm pval i_elect

    % Get observed rho value
    [permtest.rho_obs(ii), pval] = circ_corrcl(resp_errrad_cat,tmp_amp');

    % Get p-value based on Z distribution
    %   **null distribution needs to be at least approximately Gaussian
    %   *can use 2-tail only when null distribution is Gaussian, else 1-tail
%     Z_val = (permtest.rho_obs(ii) - mean(rho_perm))/std(rho_perm);
    % p_z = normcdf(Z_val); %lower-tail
%     permtest.p_z(ii) = normcdf(Z_val,'upper'); %upper-tailed
    [h,permtest.p_z(ii)] = ztest(rho_obs, mean(rho_perm), std(rho_perm)); %two-tailed

    % Get p-value based on count
%     permtest.p_n(ii) = sum(permtest.rho_obs(ii) < rho_perm)/nperms; %upper-tailed
    permtest.p_n(ii) = sum(abs(rho_obs) < abs(rho_perm))/nperms; %two-tailed
    
    clear Z_val h pval rho_perm obs_power
end
clear n_resp nperms ii












