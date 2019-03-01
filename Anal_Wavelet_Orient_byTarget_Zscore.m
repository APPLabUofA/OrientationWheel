% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% SPECTOGRAM
% A spectogram is a 3d figure that plots time on the x-axis, frequency on the 
% y-axis, and shows you the power or phase-locking value for each point. 
% We compute spectograms if we have power and phase information, averaged 
% across trials, for at least one electrode. 
% This can help us understand the changes of power and phase throughout the 
% trial.

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% Variables working with:
% ersp(i_sub,i_cond,i_perm,i_chan,:,:)
% itc(i_sub,i_cond,i_perm,i_chan,:,:)
% powbase,times,freqs

% The variables ersp and itc will be a 6D variable: 
% (participants x conditions x events x electrodes x frequencies x timepoints)
% (participants x sets x events x electrodes x frequencies x timepoints)

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% period = 1/EEG.srate; 
% time (in s) = [EEG.event.latency]*period
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

eeglab redraw


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load previously processed target-aligned epoch data
% data has been converted to log scale
all_ersp_Z = struct2cell(load('all_ersp_Z.mat'));  %gets loaded as a struct
all_ersp_Z = all_ersp_Z{1};

% load behavior data
load('ALLEEG_filt_byTargets_v3.mat');

% load settings
load('filt_byTargets_v3_Settings.mat');


% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%%                  $$$$$$$ BEH Data $$$$$$$

% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants) % --------------------
    [n,m] = size(ALLEEG(i_part).rejtrial);
    % Get list of rejected trials
    pip = 1;
    for ni = 1:n %for when there are more than 1 column
        for mi = 1:m
            if ~isempty(ALLEEG(i_part).rejtrial(ni,mi).ids)
                rejlist{pip} = ALLEEG(i_part).rejtrial(ni,mi).ids;
                pip = 1 + pip;
            end
        end
        clear mi
    end
    if pip > 1 %if trials were rejected
        err_deg_tmp = ALLEEG(i_part).error_deg; %start with all the errors
        % each set of rejected trials needs to be removed in order
        % sequentially
        for mi = 1:length(rejlist)
            tmplist = [rejlist{mi}];
            err_deg_tmp(tmplist) = []; %removes the trials
            clear tmplist
        end
        clear mi
    elseif pip == 1 %if no trials were rejected, rejlist variable not created
        err_deg_tmp = ALLEEG(i_part).error_deg;
    end
    % create variable with selected BEH 
    resp_errdeg{i_part} = err_deg_tmp;
    
    clear rejlist n m err_deg_tmp pip ni
end
clear i_part
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////

% Fit errors to mixed model
model_out = cell(1,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
    error_deg{ii} = resp_errdeg{ii};
%     error_deg{ii} = ALLEEG(ii).error_deg; %comment out if loaded data above
%     model_out{ii} = MemFit(error_deg{ii});
    model_out{ii} = MLE(error_deg{ii}); %fits without plotting
end
clear ii error_deg
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////






% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%% >>>>>>>>>>>>>>>>>>>>  POWER ANALYSES  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%% Baseline correction of power by trial and frequency band
% /////////////////////////////////////////////////////////////////////////
% #########################################################################
% Right now it is set-up to subtract mean power across the epoch at each
% frequency from each time point at that same frequency
% note: does not include the catch trial data
tic %see how long this takes
all_erspN = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % --
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        tmp_ersp = abs(all_ersp{i_part,i_elect});
        for i_trial = 1:size(tmp_ersp,3)
            for fq = 1:length(freqs)
%                 bl_freq = mean(tmp_ersp(fq,:,i_trial),2); %average each trial
                bl_freq = mean(mean(tmp_ersp(fq,:,:),2),3); %average all trials
                all_erspN{i_part,i_elect}.trials(fq,:,i_trial) = tmp_ersp(fq,:,i_trial)-bl_freq;
            end
        clear fq bl_freq
        end
    end
    clear ii i_elect tmp_ersp
end
clear i_part
toc %see how long this takes

% #########################################################################
% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%% OR use raw ERS values (log scaled)
% /////////////////////////////////////////////////////////////////////////

% --For data with targets--
all_erspN = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % --
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        tmp_ersp = abs(all_ersp{i_part,i_elect});
        for i_trial = 1:size(tmp_ersp,3)
            all_erspN{i_part,i_elect}.trials(:,:,i_trial) = 10*log10(tmp_ersp(:,:,i_trial)); %dB converted
        end
        clear i_trial
    end
    clear ii i_elect tmp_ersp
end
clear i_part


% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%%                      Standardize Power
% /////////////////////////////////////////////////////////////////////////
% #########################################################################

all_ersp_Z = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
% Change power to z-score values per person
for i_part = 1:length(exp.participants)
    % Get power across trials
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_erspN{i_part,i_elect}.trials; %get single subject's baseline corrected power
%         all_ersp_Z{i_part,i_elect}.trials = normalize(part_ersp,3,'zscore','robust');
        all_ersp_Z{i_part,i_elect}.trials = (part_ersp - mean(part_ersp(:))) / std(part_ersp(:));
        clear part_ersp i_elect
    end
    clear ii
end
clear i_part




% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%% ERS: Power by errors
% /////////////////////////////////////////////////////////////////////////
% #########################################################################

% Create ERS by errors
x_errdeg_m = cell(1,length(exp.participants)); %pre-allocate
n_errdeg_m = cell(1,length(exp.participants)); %pre-allocate
x_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
n_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
errlims = cell(1,length(exp.participants));    %pre-allocate
for i_part = 1:length(exp.participants) % ----------------------
    
    % Get upper and lower limits based on model fit
%     errlims{i_part}(1) = -(model_out{1,i_part}.maxPosterior(2)); %negative value
%     errlims{i_part}(2) = model_out{1,i_part}.maxPosterior(2);
    errlims{i_part}(1) = -(model_out{1,i_part}(2)); %negative value
    errlims{i_part}(2) = model_out{1,i_part}(2);
    
    % Get errors values
    x_errdeg_m{i_part} = resp_errdeg{i_part}(resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)); %small errors
    n_errdeg_m{i_part} = resp_errdeg{i_part}([find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))]);

    % Calculate power
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp_Z{i_part,i_elect}.trials; %get single subject's baseline corrected power
        
        % Get trials with small errors
        x_pwr{1,i_elect}(i_part,:,:) = squeeze(mean(part_ersp(:,:,[...
            find((resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)))] ),3));

        % Get trials with large errors
        n_pwr{1,i_elect}(i_part,:,:) = squeeze(mean(part_ersp(:,:,[...
            find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))] ),3));
        
        clear part_ersp i_elect
    end
end
clear ii i_part

% /////////////////////////////////////////////////////////////////////////


% #########################################################################

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%% Plot spectogram across subjects &&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% Raw ERS plots
% for ii = 1:length(exp.singletrialselecs)
for ii = 1:5
    
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    %mean across subjects
    plot_ers_x = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %small errors
    plot_ers_n = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %large errors
    
    CLim = [-1.5 1.5]; %set power scale of plot
    
    % Plot Small Errors
    figure('Position', [1 1 1685 405]); colormap('jet') %open a new figure
    subplot(1,2,1)
    imagesc(times,freqs,plot_ers_x,CLim);
    title(['Accurate: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    % Plot Large Errors
    subplot(1,2,2)
    imagesc(times,freqs,plot_ers_n,CLim);
    title(['Guesses: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
   
    clear plot_ers_x plot_ers_n CLim
end
clear ii i_elect

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Difference ERS plots
% for ii = 1:length(exp.singletrialselecs)
for ii = 1:5
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    %mean across subjects
    plot_ers_x = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %small errors
    plot_ers_n = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %large errors
    
    CLim = [-0.2 0.2]; %set power scale of plot
    
    % Plot Accurate-Guesses
    figure; colormap('jet') %open a new figure
    imagesc(times,freqs,plot_ers_x-plot_ers_n,CLim);
    title(['Accurate-Guesses: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    
    clear plot_ers_x plot_ers_n CLim
end
clear ii i_elect



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%% Plot spectogram for each subject &
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
for i_part = 1:length(exp.participants)

    figure('Position', [1 1 624 1016]); colormap('jet') %open a new figure
    
    for ii = 1:5 %just central electrodes
%     for ii = 1:length(exp.singletrialselecs)    
        
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        plot_ers_x = squeeze(x_pwr{1,i_elect}(i_part,:,:)); %small errors data for each subject
        plot_ers_n = squeeze(n_pwr{1,i_elect}(i_part,:,:)); %large errors data for each subject

        CLim = [0 2000]; %set power scale of plot

        figure('Position', [1 1 1685 405]); colormap('jet') %open a new figure
        
        % Plot Small Errors
        subplot(1,2,1)
        imagesc(times,freqs,plot_ers_x,CLim);
        title('Small Errors'); set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freq (Hz)'); xlabel('Time (ms)');
        colorbar
        
        % Plot Large Errors
        subplot(1,2,2)
        imagesc(times,freqs,plot_ers_n,CLim);
        title('Large Errors'); set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freq (Hz)'); xlabel('Time (ms)');
        colorbar

        clear plot_ers_x plot_ers_n plot_ers_c
    end
    % Overall subplot title
    supertitle(['Subj ' num2str(exp.participants{i_part}) ': ' exp.singtrlelec_name{ii}],...
        'FontSize',10.5)
    
    clear ii i_elect
end
clear i_part CLim
% -------------------------------------------------------------------------







% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%%    Compute power in time and frequency windows for errors
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------

clear x_pwr_win n_pwr_win c_pwr_win

%finds the frequencies you want (gamma (35–90 Hz))
freqband = [15 35]; %beta
% freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [10 14]; %high alpha
% freqband = [3 8]; %theta
 
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-600 -400];
% timewin = [-400 -200];
% timewin = [-200 -100];
% timewin = [-200 0];
% timewin = [0 200];
% timewin = [200 400];
% timewin = [400 600];
% timewin = [300 600];
% timewin = [100 200];
% timewin = [200 300];
% timewin = [300 500];
% timewin = [400 500];
% timewin = [500 600];
timelim = find(times>=timewin(1) & times<=timewin(2));

x_pwr_win = cell(length(exp.singletrialselecs),1); %pre-allocate
n_pwr_win = cell(length(exp.singletrialselecs),1); %pre-allocate
for i_part = 1:length(exp.participants)
%     for ii = 1:5 %only central electrodes
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        x_pwr_win{i_elect}(i_part) = mean(mean(x_pwr{1,i_elect}(i_part,freqlim,timelim),2),3); %small errors
        n_pwr_win{i_elect}(i_part) = mean(mean(n_pwr{1,i_elect}(i_part,freqlim,timelim),2),3); %large errors
    end
    clear ii i_elect
end
clear i_part

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
%% Run nonparametric statistics
% for ii = 1:length(exp.singletrialselecs)
for ii = 1:5 %test central electrodes only
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    sgnrank_pwr_win_elect{ii} = exp.singtrlelec_name{ii}; %save name
    % Sign rank test
    % accurate v guesses
    [p,h,stat] = signrank(x_pwr_win{i_elect}(:),n_pwr_win{i_elect}(:));
    sgnrank_pwr_win(ii,1) = p;
    sgnrank_pwr_win(ii,2) = h;
    sgnrank_pwr_win_stats(ii,1) = stat;
    clear h p stat i_elect
end
clear ii i_elect

% Correction w/FDR
[h,crit_p,adj_ci,adj_p] = fdr_bh(sgnrank_pwr_win(:,1),0.05);
sgnrank_pwr_win(:,3) = adj_p;
sgnrank_pwr_win(:,4) = h;
clear h crit_p adj_ci adj_p

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Run parametric statistics
% for ii = 1:length(exp.singletrialselecs)
for ii = 1:5 %test central electrodes only
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    ttest_pwr.elect{ii} = exp.singtrlelec_name{ii}; %save name
    % t-test
    [h,p,ci,stat] = ttest(x_pwr_win{i_elect}(:),n_pwr_win{i_elect}(:));
    ttest_pwr_win(ii,1) = h;
    ttest_pwr_win(ii,2) = p;
    ttest_pwr_ci_win(ii,1) = ci(1);
    ttest_pwr_ci_win(ii,2) = ci(2);
    ttest_pwr_stats_win(ii,1).electrode = stat;
    clear h p ci stat
end
clear ii i_elect

% Correction w/FDR
[h,crit_p,adj_ci,adj_p] = fdr_bh(ttest_pwr_win(:,2),0.05);
ttest_pwr_win(:,4) = adj_p;
ttest_pwr_win(:,3) = h;
clear h crit_p adj_ci adj_p






% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%% Run Permutation test

nperms = 10000; %number of permutations used to estimate null distribution
% re-set values
permtest.zval_obs = []; 
permtest.p_z = []; 
permtest.p_n = []; 

% Permutation test at each electrode
for ii = 1:5 %test central electrodes only
    
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    permtest.elect{ii,1} = exp.singtrlelec_name{i_elect};
    
    obs_pwr = n_pwr_win{i_elect}(:); %large error
    n_resp = length(obs_pwr);

    % Make distribution of null-hypothesis test statistic
    zval_perm = zeros(1,nperms); %pre-allocate
    for i_perm = 1:nperms
        order_resp = randperm(n_resp); %randomly set order of data
        [p,h,stat] = signrank(x_pwr_win{i_elect}(:),n_pwr_win{i_elect}(order_resp));
        zval_perm(i_perm) = stat.zval;
        clear order_resp p h stat
    end
    clear i_perm pval

    % Get observed zval value
    [p,h,stat] = signrank(x_pwr_win{i_elect}(:),obs_pwr);
    permtest.obs_zval(ii) = stat.zval;

    % Plot null distribution
    % figure; histogram(zval_perm)

    % Get p-value based on Z distribution
    %   **null distribution needs to be at least approximately Gaussian
    %   *can use 2-tail only when null distribution is Gaussian, else 1-tail
%     Z_val = (permtest.obs_zval(ii) - mean(zval_perm))/std(zval_perm);
    % p_z = normcdf(Z_val); %lower-tail
%     permtest.p_z(ii) = normcdf(Z_val,'upper'); %upper-tailed
    [h,permtest.p_z(ii)] = ztest(permtest.obs_zval(ii), mean(zval_perm), std(zval_perm)); %two-tailed

    % Get p-value based on count
%     permtest.p_n(ii) = sum(permtest.obs_zval(ii) < zval_perm)/nperms; %upper-tailed
    permtest.p_n(ii) = sum(abs(permtest.obs_zval(ii)) < abs(zval_perm))/nperms; %two-tailed
    
    clear Z_val h pval zval_perm obs_pwer p stat i_elect
end
clear n_resp nperms ii


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Descriptive statistics
ii = 3;
i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes

% nanmean(x_pwr_win{i_elect}(:))
% nanmean(n_pwr_win{i_elect}(:))
% 
% nanstd(x_pwr_win{i_elect}(:))
% nanstd(n_pwr_win{i_elect}(:))

bar_vals = [nanmean(x_pwr_win{i_elect}(:)) nanmean(n_pwr_win{i_elect}(:))];
bar_errs = [(nanstd(x_pwr_win{i_elect}(:))/sqrt(length(exp.participants)))...
            (nanstd(n_pwr_win{i_elect}(:))/sqrt(length(exp.participants)))];

% Bar graph
figure;
barweb(bar_vals, bar_errs); 
% ylim([0.1 0.2]);
legend('Accurate','Guesses');
title([exp.singtrlelec_name{ii} ': ' num2str(timewin(1)) ' to ' num2str(timewin(2)) ' ms; ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz'])


clear bar_vals bar_errs




% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%%              Correlate power with errors
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%           Compute power in time and frequency windows
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

clear current_power

%finds the frequencies you want
% freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [11 14]; %high alpha
freqband = [3 8]; %theta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-600 -400];
% timewin = [-400 -200];
% timewin = [-200 0];
% timewin = [0 200];
% timewin = [200 400];
% timewin = [400 600];
timelim = find(times>=timewin(1) & times<=timewin(2));

current_power = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
%     for ii = 1:5 %only central electrodes
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
%         part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        part_ersp = all_ersp_Z{i_part,i_elect}.trials; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
%             current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(part_ersp(freqlim,timelim,i_trial),1),2));
        end
        clear part_ersp i_trial

        % plot trial power on a histogram
%         figure; hist(current_power{i_part,i_elect},30)
%         ylabel('Count');
%         title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{ii}])
    end
end
clear i_elect i_part timelim 


% /////////////////////////////////////////////////////////////////////////
%%           Correlate power with degrees error each subject
% /////////////////////////////////////////////////////////////////////////
% Loop through each Participant & plot correlation
for i_part = 1:length(exp.participants) % --------------
    for ii = 1 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect});
        % correlation betwen power and errors and plot
        figure; polarscatter(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect})
        hold on
        pog = convhull(polyshape(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect}));
        polarplot(pog.Vertices(:,1),pog.Vertices(:,2))
        title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{ii}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval pog pog1
    end
    clear i_elect ii
end
clear i_part


% /////////////////////////////////////////////////////////////////////////
%%           Correlate power with degrees error
% /////////////////////////////////////////////////////////////////////////
% Create one big array of power
all_currentpwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    all_currentpwr{1,i_elect} = cat(2,current_power{1:end,i_elect});
end
clear ii i_elect

% Put all the response errors across subjects into vector
resp_errdeg_cat = cat(2,resp_errdeg{1:end});
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Plot correlation overall all subjects and trials
for ii = 1 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg_cat),all_currentpwr{1,i_elect});
    % correlation betwen power and errors and plot
    figure; polarscatter(circ_ang2rad(resp_errdeg_cat),all_currentpwr{1,i_elect},'filled','MarkerFaceAlpha',.5)
    title([exp.singtrlelec_name{ii} ': ' num2str(timewin(1)) ' to ' num2str(timewin(2)) ' ms; ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz'...
        ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
    clear x y rho pval
end
clear i_elect ii
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


% /////////////////////////////////////////////////////////////////////////
%%    Correlate power with degrees error each subject in one plot
% /////////////////////////////////////////////////////////////////////////

% Loop through each Participant & plot correlation
figure('Position', [1 1 1893 402])
for i_part = 1:length(exp.participants) % --------------
    for ii = 2 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect});
        % correlation betwen power and errors and plot
        subtightplot(2,13,i_part,[0.01 0.03],[0.001 0.001],[0.05 0.05])
        polarscatter(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect},16)
        hold on
        pog = convhull(polyshape(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect}));
        polarplot(pog.Vertices(:,1),pog.Vertices(:,2))
        
%         title(['Subj ' num2str(exp.participants{i_part}) ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        title(['rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))],'FontSize',8)
        clear x y rho pval pog
    end
    % Overall subplot title
    supertitle([exp.singtrlelec_name{ii} ': ' num2str(timewin(1)) ' to ' num2str(timewin(2)) ' ms; '...
        num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz'])
    
    clear i_elect ii
end
clear i_part

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% -------------------------------------------------------------------------
%                     --- Permutation test ---
% -------------------------------------------------------------------------
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


nperms = 10000; %number of permutations used to estimate null distribution
resp_errrad_cat = circ_ang2rad(resp_errdeg_cat); %convert to radians for circular corr
n_resp = length(resp_errrad_cat); %number of observations to permute

% re-set values
permtest.rho_obs = []; 
permtest.p_z = []; 
permtest.p_n = []; 

% Permutation test at each electrode
for ii = 1:5 %test central electrodes only
    
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    permtest.elect{ii,1} = exp.singtrlelec_name{i_elect};
    obs_power = all_currentpwr{1,i_elect}; %get power data from electrode

    % Make distribution of null-hypothesis test statistic
    rho_perm = zeros(1,nperms); %pre-allocate
    for i_perm = 1:nperms
        order_resp = randperm(n_resp); %randomly set order of data
        [rho_perm(i_perm), pval] = circ_corrcl(resp_errrad_cat(order_resp),obs_power);
        clear order_resp
    end
    clear i_perm pval i_elect

    % Get observed rho value
    [permtest.rho_obs(ii), pval] = circ_corrcl(resp_errrad_cat,obs_power);

    % Plot null distribution
    % figure; histogram(rho_perm)

    % Get p-value based on Z distribution
    %   **null distribution needs to be at least approximately Gaussian
    %   *can use 2-tail only when null distribution is Gaussian, else 1-tail
    Z_val = (permtest.rho_obs(ii) - mean(rho_perm))/std(rho_perm);
    % p_z = normcdf(Z_val); %lower-tail
    permtest.p_z(ii) = normcdf(Z_val,'upper'); %upper-tailed
    % [h,p_z] = ztest(rho_obs, mean(rho_perm), std(rho_perm)); %two-tailed

    % Get p-value based on count
    permtest.p_n(ii) = sum(permtest.rho_obs(ii) < rho_perm)/nperms; %upper-tailed
    % p_n = sum(abs(rho_obs) < abs(rho_perm))/nperms; %two-tailed
    
    clear Z_val h pval rho_perm obs_power
end

clear n_resp nperms ii


% Correction w/FDR
[h,crit_p,adj_ci,permtest.adj_pz] = fdr_bh(permtest.p_z,0.05);
[h,crit_p,adj_ci,permtest.adj_pn] = fdr_bh(permtest.p_n,0.05);
clear h crit_p adj_ci adj_p




% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     ''''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

% Set the range of time to consider
tWin{1} = [-600 -400];
tWin{2} = [-400 -200];
tWin{3} = [-200 0];
tWin{4} = [0 200];
tWin{5} = [200 400];
tWin{6} = [400 600];

% tWin{1} = [-600 -300];
% tWin{2} = [-300 0];
% tWin{3} = [0 300];
% tWin{4} = [300 600];

%finds the frequencies you want
freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [10 14]; %high alpha
% freqband = [3 8]; %theta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

% ERSP averaged across subjects
% ersp(participants x conditions x events x electrodes x frequencies x timepoints)
out_ersp_elect = squeeze(mean(ersp(:,1,1,:,:,:),1));

CLim = [40 70]; %set power scale of plot
colormap('jet')

for tw_i = 1:length(tWin) %loop through several time windows 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %finds the times you want from the times variable
    time_window = find(times>= itWin(1),1):find(times>= itWin(2),1)-1;
    
%     figure('Color',[1 1 1],'Position',[1 1 941 349]);
    figure('Color',[1 1 1]);
    set(gca,'Color',[1 1 1]);
    
    temp = mean(mean(out_ersp_elect(:,freqlim,time_window),2),3)';
    temp(1) = NaN; %not M2 electrode
       
    topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
        'plotchans',elect_erp,'emarker',{'.','k',11,1})

    title([num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz: ' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms']);
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Power (uV^2)');

    clear itWin time_window temp
end
clear tw_i

clear freqlim CLim



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Correlate power in time and frequency windows with response error
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

clear current_power

%finds the frequencies you want
freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [10 14]; %high alpha
% freqband = [3 8]; %theta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
% timewin = [-600 -400];
% timewin = [-400 -200];
timewin = [-200 0];
% timewin = [0 200];
% timewin = [200 400];
% timewin = [400 600];
timelim = find(times>=timewin(1) & times<=timewin(2));

current_power = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
%     for ii = 1:5 %only central electrodes
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
%         part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        part_ersp = all_ersp_Z{i_part,i_elect}.trials; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
%             current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(part_ersp(freqlim,timelim,i_trial),1),2));
        end
        clear part_ersp i_trial

        % plot trial power on a histogram
%         figure; hist(current_power{i_part,i_elect},30)
%         ylabel('Count');
%         title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{ii}])
    end
end
clear i_elect i_part timelim 


% /////////////////////////////////////////////////////////////////////////
%%           Correlate power with degrees error each subject
% /////////////////////////////////////////////////////////////////////////
% Loop through each participant & plot correlation
for i_part = 1:length(exp.participants) % --------------
    for ii = 1 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect});
        % correlation betwen power and errors and plot
        figure; polarscatter(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect})
        title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{ii}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
    end
    clear i_elect ii
end
clear i_part


% /////////////////////////////////////////////////////////////////////////
%%           Correlate power with degrees error
% /////////////////////////////////////////////////////////////////////////
% Create one big array of power
all_currentpwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    all_currentpwr{1,i_elect} = cat(2,current_power{1:end,i_elect});
end
clear ii i_elect

% Put all the response errors across subjects into vector
resp_errdeg_cat = cat(2,resp_errdeg{1:end});
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Plot correlation overall all subjects and trials
for ii = 1:5 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg_cat),all_currentpwr{1,i_elect});
    % correlation betwen power and errors and plot
    figure; polarscatter(circ_ang2rad(resp_errdeg_cat),all_currentpwr{1,i_elect},'filled','MarkerFaceAlpha',.5)
    title([exp.singtrlelec_name{ii} ': ' num2str(timewin(1)) ' to ' num2str(timewin(2)) ' ms; ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz'...
        ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
    clear x y rho pval
end
clear i_elect ii
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% -------------------------------------------------------------------------
%                     --- Permutation test ---
% -------------------------------------------------------------------------
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


nperms = 10000; %number of permutations used to estimate null distribution
resp_errrad_cat = circ_ang2rad(resp_errdeg_cat); %convert to radians for circular corr
n_resp = length(resp_errrad_cat); %number of observations to permute

% re-set values
permtest.rho_obs = []; 
permtest.p_z = []; 
permtest.p_n = []; 

% Permutation test at each electrode
for ii = 1:5 %test central electrodes only
    
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    permtest.elect{ii,1} = exp.singtrlelec_name{i_elect};
    obs_power = all_currentpwr{1,i_elect}; %get power data from electrode

    % Make distribution of null-hypothesis test statistic
    rho_perm = zeros(1,nperms); %pre-allocate
    for i_perm = 1:nperms
        order_resp = randperm(n_resp); %randomly set order of data
        [rho_perm(i_perm), pval] = circ_corrcl(resp_errrad_cat(order_resp),obs_power);
        clear order_resp
    end
    clear i_perm pval i_elect

    % Get observed rho value
    [permtest.rho_obs(ii), pval] = circ_corrcl(resp_errrad_cat,obs_power);

    % Plot null distribution
    % figure; histogram(rho_perm)

    % Get p-value based on Z distribution
    %   **null distribution needs to be at least approximately Gaussian
    %   *can use 2-tail only when null distribution is Gaussian, else 1-tail
    Z_val = (permtest.rho_obs(ii) - mean(rho_perm))/std(rho_perm);
    % p_z = normcdf(Z_val); %lower-tail
    permtest.p_z(ii) = normcdf(Z_val,'upper'); %upper-tailed
    % [h,p_z] = ztest(rho_obs, mean(rho_perm), std(rho_perm)); %two-tailed

    % Get p-value based on count
    permtest.p_n(ii) = sum(permtest.rho_obs(ii) < rho_perm)/nperms; %upper-tailed
    % p_n = sum(abs(rho_obs) < abs(rho_perm))/nperms; %two-tailed
    
    clear Z_val h pval rho_perm obs_power
end

clear n_resp nperms ii






% #########################################################################
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%%            Separate response errors by power
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
% #########################################################################

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%           Compute power in time and frequency windows
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

clear current_power 

%finds the frequencies you want
% freqband = [15 35]; %beta
freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [10 14]; %high alpha
% freqband = [3 8]; %theta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

%finds the times you want from the times variable
% timewin = [-600 -500];
% timewin = [-500 -400];
% timewin = [-400 -300];
% timewin = [-300 -200];
% timewin = [-200 -100];
% timewin = [-100 0];
% timewin = [0 100];
% timewin = [100 200];
% timewin = [200 300];
% timewin = [300 400];
% timewin = [400 500];
% timewin = [500 600];

% timewin = [-600 -400];
% timewin = [-400 -200];
% timewin = [-200 0];
timewin = [0 200];
% timewin = [200 400];
% timewin = [400 600];
timelim = find(times>=timewin(1) & times<=timewin(2));
% -------------------------------------------------------------------------
% Calculate power in time window at freq band
current_power = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
%     for ii = 1:5 %only central electrodes
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
%         part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        part_ersp = all_ersp_Z{i_part,i_elect}.trials; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
%             current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(part_ersp(freqlim,timelim,i_trial),1),2));
        end
        clear part_ersp i_trial

        % plot trial power on a histogram
%         figure; hist(current_power{i_part,i_elect},30)
%         ylabel('Count');
%         title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{ii}])
    end
end
clear i_elect i_part timelim 

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Separate trials based on median band power
errdeg_Hpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
errdeg_Lpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
pwr_Hpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
pwr_Lpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
%     for ii = 1:5 %only central electrodes
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % Get trials by median power split
        mdn_pwr = nanmedian(current_power{i_part,i_elect});
        errdeg_Hpwr{i_part,i_elect} = resp_errdeg{i_part}(current_power{i_part,i_elect}>=mdn_pwr);
        errdeg_Lpwr{i_part,i_elect} = resp_errdeg{i_part}(current_power{i_part,i_elect}<mdn_pwr);
        % get high/low power on trials
        pwr_Hpwr{i_part,i_elect} = current_power{i_part,i_elect}(current_power{i_part,i_elect}>=mdn_pwr);
        pwr_Lpwr{i_part,i_elect} = current_power{i_part,i_elect}(current_power{i_part,i_elect}<mdn_pwr);
        clear i_elect mdn_pwr
    end
    clear ii
end
clear i_part

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Fit errors to model by electrode and subject
model_out_Hpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
model_out_Lpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
mean_errdeg_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
mean_errdeg_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
std_errdeg_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
std_errdeg_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
g_out_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
g_out_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
sd_out_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
sd_out_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
mu_out_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
mu_out_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
% Fit errors to mixed model
model = StandardMixtureModel(); %standard 2 parameter model
% model = WithBias(StandardMixtureModel); %model with mu
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
%     for ii = 1:5 %only central electrodes
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % Get mean error
        mean_errdeg_Hpwr(i_part,i_elect) = mean(abs(errdeg_Hpwr{i_part,i_elect}));
        mean_errdeg_Lpwr(i_part,i_elect) = mean(abs(errdeg_Lpwr{i_part,i_elect}));
        % Get SD error
        std_errdeg_Hpwr(i_part,i_elect) = std(abs(errdeg_Hpwr{i_part,i_elect}));
        std_errdeg_Lpwr(i_part,i_elect) = std(abs(errdeg_Lpwr{i_part,i_elect}));
        
%        model_out_Hpwr{i_part,i_elect} = MemFit(errdeg_Hpwr{i_part,i_elect},model);
%        model_out_Lpwr{i_part,i_elect} = MemFit(errdeg_Lpwr{i_part,i_elect},model);
        model_out_Hpwr{i_part,i_elect} = MLE(errdeg_Hpwr{i_part,i_elect},model); %fits without plotting
        model_out_Lpwr{i_part,i_elect} = MLE(errdeg_Lpwr{i_part,i_elect},model); %fits without plotting
        
        % Save specific outputs
        if strcmpi(model.name,'Standard mixture model with bias') %model w/mu
            mu_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(1);
            mu_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(1);
            g_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(2);
            g_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(2);
            sd_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(3);
            sd_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(3);
        else
            g_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(1);
            g_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(1);
            sd_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(2);
            sd_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(2);
        end
        
        clear i_elect
    end
    clear ii
end
clear i_part

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Run nonparametric statistics
for ii = 1:length(exp.singletrialselecs)
% for ii = 1:5 %only central electrodes
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    
    [p,h,stat] = signrank(mean_errdeg_Hpwr(:,i_elect),mean_errdeg_Lpwr(:,i_elect));
    sgnrank_pwr.mean(ii,1) = p;
    sgnrank_pwr.mean(ii,2) = h;
    sgnrank_pwr_stats.mean(ii,1).electrode = stat;
    clear h p stat
    
    [p,h,stat] = signrank(std_errdeg_Hpwr(:,i_elect),std_errdeg_Lpwr(:,i_elect));
    sgnrank_pwr.std(ii,1) = p;
    sgnrank_pwr.std(ii,2) = h;
    sgnrank_pwr_stats.std(ii,1).electrode = stat;
    clear h p stat
    
    
    [p,h,stat] = signrank(g_out_Hpwr(:,i_elect),g_out_Lpwr(:,i_elect));
    sgnrank_pwr.g(ii,1) = p;
    sgnrank_pwr.g(ii,2) = h;
    sgnrank_pwr_stats.g(ii,1).electrode = stat;
    clear h p stat

    [p,h,stat] = signrank(sd_out_Hpwr(:,i_elect),sd_out_Lpwr(:,i_elect));
    sgnrank_pwr.sd(ii,1) = p;
    sgnrank_pwr.sd(ii,2) = h;
    sgnrank_pwr_stats.sd(ii,1).electrode = stat;
    clear h p stat
    
    if strcmpi(model.name,'Standard mixture model with bias') %model w/mu
        [p,h,stat] = signrank(mu_out_Hpwr(:,i_elect),mu_out_Lpwr(:,i_elect));
        sgnrank_pwr.mu(ii,1) = p;
        sgnrank_pwr.mu(ii,2) = h;
        sgnrank_pwr_stats.mu(ii,1).electrode = stat;
        clear h p stat
    end
    clear i_elect
end
clear ii

% Turn into table
T = struct2table(sgnrank_pwr);

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Run parametric statistics
for ii = 1:length(exp.singletrialselecs)
% for ii = 1:5 %only central electrodes
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    ttest_pwr.elect{ii} = exp.singtrlelec_name{ii}; %save name
%     nanmean(g_out_Hpwr(:,i_elect))
%     nanmean(g_out_Lpwr(:,i_elect))
%     
%     nanmean(sd_out_Hpwr(:,i_elect))
%     nanmean(sd_out_Lpwr(:,i_elect))

    [h,p,ci,stat] = ttest(mean_errdeg_Hpwr(:,i_elect),mean_errdeg_Lpwr(:,i_elect));
    ttest_pwr.mean(ii,1) = h;
    ttest_pwr.mean(ii,2) = p;
    ttest_pwr_ci.mean(ii,1) = ci(1);
    ttest_pwr_ci.mean(ii,2) = ci(2);
    ttest_pwr_stats.mean(ii,1).electrode = stat;
    clear h p ci stat
    
    [h,p,ci,stat] = ttest(std_errdeg_Hpwr(:,i_elect),std_errdeg_Lpwr(:,i_elect));
    ttest_pwr.std(ii,1) = h;
    ttest_pwr.std(ii,2) = p;
    ttest_pwr_ci.std(ii,1) = ci(1);
    ttest_pwr_ci.std(ii,2) = ci(2);
    ttest_pwr_stats.std(ii,1).electrode = stat;
    clear h p ci stat
    
    
    [h,p,ci,stat] = ttest(g_out_Hpwr(:,i_elect),g_out_Lpwr(:,i_elect));
    ttest_pwr.g(ii,1) = h;
    ttest_pwr.g(ii,2) = p;
    ttest_pwr_ci.g(ii,1) = ci(1);
    ttest_pwr_ci.g(ii,2) = ci(2);
    ttest_pwr_stats.g(ii,1).electrode = stat;
    clear h p ci stat

    [h,p,ci,stat] = ttest(sd_out_Hpwr(:,i_elect),sd_out_Lpwr(:,i_elect));
    ttest_pwr.sd(ii,1) = h;
    ttest_pwr.sd(ii,2) = p;
    ttest_pwr_ci.sd(ii,1) = ci(1);
    ttest_pwr_ci.sd(ii,2) = ci(2);
    ttest_pwr_stats.sd(ii,1).electrode = stat;
    clear h p ci stat
    
    if strcmpi(model.name,'Standard mixture model with bias') %model w/mu
        [h,p,ci,stat] = ttest(mu_out_Hpwr(:,i_elect),mu_out_Lpwr(:,i_elect));
        ttest_pwr.mu(ii,1) = h;
        ttest_pwr.mu(ii,2) = p;
        ttest_pwr_ci.mu(ii,1) = ci(1);
        ttest_pwr_ci.mu(ii,2) = ci(2);
        ttest_pwr_stats.mu(ii,1).electrode = stat;
        clear h p ci stat
    end
    clear i_elect
end
clear ii
% -------------------------------------------------------------------------
%% Add p-values to table
if strcmpi(model.name,'Standard mixture model with bias') %model w/mu
    T.mean(:,2) = ttest_pwr.mean(:,2);
    T.std(:,2) = ttest_pwr.std(:,2);
    T.g(:,2) = ttest_pwr.g(:,2);
    T.sd(:,2) = ttest_pwr.sd(:,2);
    T.mu(:,2) = ttest_pwr.mu(:,2);
    T.elect = ttest_pwr.elect'; % electrode labels
else
    T.mean(:,2) = ttest_pwr.mean(:,2);
    T.std(:,2) = ttest_pwr.std(:,2);
    T.g(:,2) = ttest_pwr.g(:,2);
    T.sd(:,2) = ttest_pwr.sd(:,2);
    T.elect = ttest_pwr.elect'; % electrode labels
end
    

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

% clears variables that end with...
clear -regexp _Hpwr\> _Lpwr\>
clear sgnrank_pwr ttest_pwr ttest_pwr_stats ttest_pwr_ci sgnrank_pwr_stats sgnrank_pwr_ci...
    T


% Correction w/FDR
[h,crit_p,adj_ci,adj_pz] = fdr_bh(T.g(:,2),0.05);
[h,crit_p,adj_ci,adj_pn] = fdr_bh(T.g(:,1),0.05);
clear h crit_p adj_ci adj_p adj_pz adj_pn

i_elect=exp.singletrialselecs(ii);

mean(g_out_Hpwr(:,i_elect))
mean(g_out_Lpwr(:,i_elect))

std(g_out_Hpwr(:,i_elect))
std(g_out_Lpwr(:,i_elect))


mean(sd_out_Hpwr(:,i_elect))
mean(sd_out_Lpwr(:,i_elect))

std(sd_out_Hpwr(:,i_elect))
std(sd_out_Lpwr(:,i_elect))



mean(mean_errdeg_Hpwr(:,i_elect))
mean(mean_errdeg_Lpwr(:,i_elect))





% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     '''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};


clear current_power 

%finds the frequencies you want
% freqband = [15 35]; %beta
freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [10 14]; %high alpha
% freqband = [3 8]; %theta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

%finds the times you want from the times variable
% timewin = [-600 -500];
% timewin = [-500 -400];
% timewin = [-400 -300];
% timewin = [-300 -200];
% timewin = [-200 -100];
% timewin = [-100 0];
% timewin = [0 100];
% timewin = [100 200];
% timewin = [200 300];
% timewin = [300 400];
% timewin = [400 500];
% timewin = [500 600];

% timewin = [-600 -400];
% timewin = [-400 -200];
timewin = [-200 0];
% timewin = [0 200];
% timewin = [200 400];
% timewin = [400 600];
timelim = find(times>=timewin(1) & times<=timewin(2));
% -------------------------------------------------------------------------
% Calculate power in time window at freq band
current_power = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
%     for ii = 1:5 %only central electrodes
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
%         part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        part_ersp = all_ersp_Z{i_part,i_elect}.trials; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
%             current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(part_ersp(freqlim,timelim,i_trial),1),2));
        end
        clear part_ersp i_trial

        % plot trial power on a histogram
%         figure; hist(current_power{i_part,i_elect},30)
%         ylabel('Count');
%         title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{ii}])
    end
end
clear i_elect i_part timelim 

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Separate trials based on median band power
errdeg_Hpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
errdeg_Lpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
pwr_Hpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
pwr_Lpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
%     for ii = 1:5 %only central electrodes
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % Get trials by median power split
        mdn_pwr = nanmedian(current_power{i_part,i_elect});
        errdeg_Hpwr{i_part,i_elect} = resp_errdeg{i_part}(current_power{i_part,i_elect}>=mdn_pwr);
        errdeg_Lpwr{i_part,i_elect} = resp_errdeg{i_part}(current_power{i_part,i_elect}<mdn_pwr);
        % get high/low power on trials
        pwr_Hpwr{i_part,i_elect} = current_power{i_part,i_elect}(current_power{i_part,i_elect}>=mdn_pwr);
        pwr_Lpwr{i_part,i_elect} = current_power{i_part,i_elect}(current_power{i_part,i_elect}<mdn_pwr);
        clear i_elect mdn_pwr
    end
    clear ii
end
clear i_part

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Fit errors to model by electrode and subject
model_out_Hpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
model_out_Lpwr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
mean_errdeg_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
mean_errdeg_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
std_errdeg_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
std_errdeg_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
g_out_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
g_out_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
sd_out_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
sd_out_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
mu_out_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
mu_out_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
% Fit errors to mixed model
model = StandardMixtureModel(); %standard 2 parameter model
% model = WithBias(StandardMixtureModel); %model with mu
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
%     for ii = 1:5 %only central electrodes
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % Get mean error
        mean_errdeg_Hpwr(i_part,i_elect) = mean(abs(errdeg_Hpwr{i_part,i_elect}));
        mean_errdeg_Lpwr(i_part,i_elect) = mean(abs(errdeg_Lpwr{i_part,i_elect}));
        % Get SD error
        std_errdeg_Hpwr(i_part,i_elect) = std(abs(errdeg_Hpwr{i_part,i_elect}));
        std_errdeg_Lpwr(i_part,i_elect) = std(abs(errdeg_Lpwr{i_part,i_elect}));
        
%        model_out_Hpwr{i_part,i_elect} = MemFit(errdeg_Hpwr{i_part,i_elect},model);
%        model_out_Lpwr{i_part,i_elect} = MemFit(errdeg_Lpwr{i_part,i_elect},model);
        model_out_Hpwr{i_part,i_elect} = MLE(errdeg_Hpwr{i_part,i_elect},model); %fits without plotting
        model_out_Lpwr{i_part,i_elect} = MLE(errdeg_Lpwr{i_part,i_elect},model); %fits without plotting
        
        % Save specific outputs
        if strcmpi(model.name,'Standard mixture model with bias') %model w/mu
            mu_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(1);
            mu_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(1);
            g_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(2);
            g_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(2);
            sd_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(3);
            sd_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(3);
        else
            g_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(1);
            g_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(1);
            sd_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(2);
            sd_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(2);
        end
        
        clear i_elect
    end
    clear ii
end
clear i_part

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Create topoplots
% Guess rate topoplots
CLim = [0.15 0.25]; %set power scale of plot
colormap('warm')
    
figure('Color',[1 1 1],'Position',[1 1 941 349]);
temp = mean(g_out_Hpwr(:,:),1)';
temp(1) = NaN; %not M2 electrode
subplot(1,2,1);
topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
    'plotchans',elect_erp,'emarker',{'.','k',11,1})
title('High Power');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Guess Rate (g)');
clear temp

temp = mean(g_out_Lpwr(:,:),1)';
temp(1) = NaN; %not M2 electrode
subplot(1,2,2);
topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
    'plotchans',elect_erp,'emarker',{'.','k',11,1})
title('Low Power');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Guess Rate (g)');
clear temp
    
% Overall subplot title
supertitle(['Guess Rate: ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz: ' num2str(timewin(1)) ' to ' num2str(timewin(2)) ' ms'],...
    'FontSize',10.5)

clear CLim

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SD topoplots
CLim = [9 13]; %set power scale of plot
colormap('jet')
    
figure('Color',[1 1 1],'Position',[1 1 941 349]);
temp = mean(sd_out_Hpwr(:,:),1)';
temp(1) = NaN; %not M2 electrode
subplot(1,2,1);
topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
    'plotchans',elect_erp,'emarker',{'.','k',11,1})
title('High Power');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Quality (SD)');
clear temp

temp = mean(sd_out_Lpwr(:,:),1)';
temp(1) = NaN; %not M2 electrode
subplot(1,2,2);
topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
    'plotchans',elect_erp,'emarker',{'.','k',11,1})
title('Low Power');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Quality (SD)');
clear temp
    
% Overall subplot title
supertitle(['SD: ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz: ' num2str(timewin(1)) ' to ' num2str(timewin(2)) ' ms'],...
    'FontSize',10.5)

clear CLim

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Mu topoplots
CLim = [-1.5 1.5]; %set power scale of plot
colormap('jet')
    
figure('Color',[1 1 1],'Position',[1 1 941 349]);
temp = mean(mu_out_Hpwr(:,:),1)';
temp(1) = NaN; %not M2 electrode
subplot(1,2,1);
topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
    'plotchans',elect_erp,'emarker',{'.','k',11,1})
title('High Power');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Bias (mu)');
clear temp

temp = mean(mu_out_Lpwr(:,:),1)';
temp(1) = NaN; %not M2 electrode
subplot(1,2,2);
topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
    'plotchans',elect_erp,'emarker',{'.','k',11,1})
title('Low Power');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Bias (mu)');
clear temp
    
% Overall subplot title
supertitle(['Mu (3 params): ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz: ' num2str(timewin(1)) ' to ' num2str(timewin(2)) ' ms'],...
    'FontSize',10.5)

clear CLim

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Difference topoplots
% Guess rate topoplots
CLim = [-0.02 0.02]; %set power scale of plot

temp = mean((g_out_Hpwr(:,:)-g_out_Lpwr(:,:)),1)';
temp(1) = NaN; %not M2 electrode
topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
    'plotchans',elect_erp,'emarker',{'.','k',11,1})
colormap('summer')
title('High Power - Low Power');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Guess Rate (g)');
clear temp




% clears variables that end with...
clear -regexp _Hpwr\> _Lpwr\>
clear sgnrank_pwr ttest_pwr ttest_pwr_stats ttest_pwr_ci sgnrank_pwr_stats sgnrank_pwr_ci...
    T







