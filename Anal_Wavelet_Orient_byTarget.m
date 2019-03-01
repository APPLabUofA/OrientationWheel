
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
all_erspN = struct2cell(load('all_erspN.mat')); %gets loaded as a struct
all_erspN = all_erspN{1}; 

% load behavior data
load('ALLEEG_filt_byTargets_v3.mat');

% load settings
load('filt_byTargets_v3_Settings.mat');


% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Load previously processed catch trial data
% assumes that the catch trials have been processed with the same
% parameters as the all_ersp data (check their settings in exp to make sure) 
catch_trials_v1_wav = load("all_ersp_catch_trials.mat");
all_ersp_byC = catch_trials_v1_wav.all_ersp;

% Check to make sure their frequency and time scales are the same
% output should be 0 or they are not the same!
sum(catch_trials_v1_wav.freqs ~= freqs)
sum(catch_trials_v1_wav.times ~= times)

clear catch_trials_v1_wav


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
model1 = StandardMixtureModel();
model2 = WithBias(StandardMixtureModel);
model_out = cell(1,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
    error_deg{ii} = resp_errdeg{ii};
%     error_deg{ii} = ALLEEG(ii).error_deg; %comment out if loaded data above
%     model_out{ii} = MemFit(error_deg{ii},model1);
%     MemFit(error_deg{ii},{model1,model2}); %comparing models
    model_out{ii} = MLE(error_deg{ii},model1); %fits without plotting
%     model_out{ii} = MLE(error_deg{ii},model2); %fits without plotting
%     model_out{ii} = model_out{ii}(2:3);
end
clear ii error_deg model1 model2


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


% --For catch trial data--
all_erspC = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % --
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        tmp_ersp = abs(all_ersp_byC{i_part,i_elect});
        for i_trial = 1:size(tmp_ersp,3)
            all_erspC{i_part,i_elect}.trials(:,:,i_trial) = 10*log10(tmp_ersp(:,:,i_trial)); %dB converted
        end
        clear i_trial
    end
    clear ii i_elect tmp_ersp
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
        part_ersp = all_erspN{i_part,i_elect}.trials; %get single subject's baseline corrected power
        
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

% Create ERS for catch trials
c_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % ----------------------
    % Calculate power
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_erspC{i_part,i_elect}.trials; %get single subject's baseline corrected power
        c_pwr{1,i_elect}(i_part,:,:) = squeeze(mean(part_ersp(:,:,:),3));
    end
    clear part_ersp i_elect
end
clear ii i_part

% /////////////////////////////////////////////////////////////////////////
% #########################################################################
%% Gets a count of trials
err_trl_count(:,1) = cellfun(@numel,x_errdeg_m); %small errors
err_trl_count(:,2) = cellfun(@numel,n_errdeg_m); %large errors
err_trl_count(:,3) = cell2mat({ALLEEG(1:end).trials}); %total trial count
err_trl_count(:,4) = cell2mat({ALLEEG(1:end).catch_trials}); %catch trials
% #########################################################################
% #########################################################################


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% /////////////////////////////////////////////////////////////////////////
%% >>> Plot spectogram across subjects 
% /////////////////////////////////////////////////////////////////////////

% Raw ERS plots
% for ii = 1:length(exp.singletrialselecs)
for ii = 1:5
    
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    %mean across subjects
    plot_ers_x = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %small errors
    plot_ers_n = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %large errors
    plot_ers_c = squeeze(mean(c_pwr{1,i_elect}(:,:,:),1)); %catch trials
    
    CLim = [20 35]; %set power scale of plot
    
    % Plot Small Errors
    figure('Position', [1 1 1685 405]); colormap('jet') %open a new figure
    subplot(1,3,1)
    imagesc(times,freqs,plot_ers_x,CLim);
    title(['Accurate: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    % Plot Large Errors
    subplot(1,3,2)
    imagesc(times,freqs,plot_ers_n,CLim);
    title(['Guesses: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    % Plot Catch Trials
    subplot(1,3,3)
    imagesc(times,freqs,plot_ers_c,CLim);
    title(['No Target: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar

    clear plot_ers_x plot_ers_n plot_ers_c CLim
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
    plot_ers_c = squeeze(mean(c_pwr{1,i_elect}(:,:,:),1)); %catch trials
    
    CLim = [-0.8 0.8]; %set power scale of plot
    
    % Plot Accurate-Guesses
    figure('Position', [1 1 1685 405]); colormap('jet') %open a new figure
    subplot(1,3,1)
    imagesc(times,freqs,plot_ers_x-plot_ers_n,CLim);
    title(['Accurate-Guesses: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    % Plot Guesses-No Target
    subplot(1,3,2)
    imagesc(times,freqs,plot_ers_x-plot_ers_c,CLim);
    title(['Accurate-No Target: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    % Plot Accurate-No Target
    subplot(1,3,3)
    imagesc(times,freqs,plot_ers_n-plot_ers_c,CLim);
    title(['Guesses-No Target: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([3 35]); yticks(5:5:35)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    
    clear plot_ers_x plot_ers_n plot_ers_c CLim
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
        plot_ers_c = squeeze(c_pwr{1,i_elect}(i_part,:,:)); %large errors data for each subject

        CLim = [0 2000]; %set power scale of plot

        figure('Position', [1 1 1685 405]); colormap('jet') %open a new figure
        
        % Plot Small Errors
        subplot(1,3,1)
        imagesc(times,freqs,plot_ers_x,CLim);
        title('Small Errors'); set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freq (Hz)'); xlabel('Time (ms)');
        colorbar
        
        % Plot Large Errors
        subplot(1,3,2)
        imagesc(times,freqs,plot_ers_n,CLim);
        title('Large Errors'); set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freq (Hz)'); xlabel('Time (ms)');
        colorbar
    
        % Plot Catch Trials
        subplot(1,3,3)
        imagesc(times,freqs,plot_ers_c);
        title('No Target'); set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line for target onset
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freqency (Hz)'); xlabel('Time (ms)');
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

%finds the frequencies you want
freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [11 14]; %high alpha
% freqband = [3 8]; %theta
% freqband = [15 35]; %beta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-600 -400];
% timewin = [-400 -200];
% timewin = [-200 -100];
% timewin = [-200 0];
% timewin = [-100 0];
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
c_pwr_win = cell(length(exp.singletrialselecs),1); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:5 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        x_pwr_win{i_elect}(i_part) = mean(mean(x_pwr{1,i_elect}(i_part,freqlim,timelim),2),3); %small errors
        n_pwr_win{i_elect}(i_part) = mean(mean(n_pwr{1,i_elect}(i_part,freqlim,timelim),2),3); %large errors
        c_pwr_win{i_elect}(i_part) = mean(mean(c_pwr{1,i_elect}(i_part,freqlim,timelim),2),3); %no target
    end
    clear ii i_elect
end
clear i_part

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
%% Run nonparametric statistics
% for ii = 1:length(exp.singletrialselecs)
for ii = 1:5 %only central electrodes
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    sgnrank_pwr_win_elect{ii} = exp.singtrlelec_name{ii}; %save name
    % Sign rank test
    % accurate v guesses
    [p,h,stat] = signrank(x_pwr_win{i_elect}(:),n_pwr_win{i_elect}(:));
    sgnrank_pwr_win.AvG(ii,1) = p;
    sgnrank_pwr_win.AvG(ii,2) = h;
    sgnrank_pwr_win_stats.AvG(ii,1) = stat;
    clear h p stat
    % accurate v catch trials
    [p,h,stat] = signrank(x_pwr_win{i_elect}(:),c_pwr_win{i_elect}(:));
    sgnrank_pwr_win.AvC(ii,1) = p;
    sgnrank_pwr_win.AvC(ii,2) = h;
    sgnrank_pwr_win_stats.AvC(ii,1) = stat;
    clear h p stat
    % accurate v catch trials
    [p,h,stat] = signrank(n_pwr_win{i_elect}(:),c_pwr_win{i_elect}(:));
    sgnrank_pwr_win.GvC(ii,1) = p;
    sgnrank_pwr_win.GvC(ii,2) = h;
    sgnrank_pwr_win_stats.GvC(ii,1) = stat;
    clear h p stat i_elect
end
clear ii i_elect

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Run parametric statistics
% for ii = 1:length(exp.singletrialselecs)
%     i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
%     ttest_pwr.elect{ii} = exp.singtrlelec_name{ii}; %save name
%     % t-test
%     [h,p,ci,stat] = ttest(x_pwr_win{i_elect}(:),n_pwr_win{i_elect}(:));
%     ttest_pwr_win(ii,1) = h;
%     ttest_pwr_win(ii,2) = p;
%     ttest_pwr_ci_win(ii,1) = ci(1);
%     ttest_pwr_ci_win(ii,2) = ci(2);
%     ttest_pwr_stats_win(ii,1).electrode = stat;
%     clear h p ci stat
% end
% clear ii i_elect


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
% freqband = [7 13]; %alpha
freqband = [8 14]; %alpha
% freqband = [7 10]; %low alpha
% freqband = [10 14]; %high alpha
% freqband = [4 7]; %theta
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
        part_ersp = all_erspN{i_part,i_elect}.trials; %get single subject's ersp
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
for i_part = 1:length(exp.participants)
    for ii = 1 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect});
        % correlation betwen power and errors and plot
        figure; polarscatter(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect})
        title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{ii}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval pog
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



% /////////////////////////////////////////////////////////////////////////
%%    Correlate power with degrees error each subject in one plot
% /////////////////////////////////////////////////////////////////////////

% Loop through each Participant & plot correlation
figure('Position', [1 1 1893 402])
for i_part = 1:length(exp.participants) % --------------
    for ii = 1 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect});
        % correlation betwen power and errors and plot
        subtightplot(2,13,i_part,[0.01 0.03],[0.001 0.001],[0.05 0.05])
        polarscatter(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect},16)
%         title(['Subj ' num2str(exp.participants{i_part}) ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
    end
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





% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     ''''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% List electrodes to get topograph plots (need all of them) 
elect_topplot = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
% el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
%     'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};


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

% Finds the frequencies you want from the freqs variable
freqband = [8 14]; %alpha
% freqband = [3 8]; %theta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

% Get mean power at frequency band for each electrode
pwr_top = NaN([length(times),length(elect_topplot),3]); %pre-allocate
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes

    pwr_top(:,i_elect,1) = squeeze(mean(mean(x_pwr{1,i_elect}(:,freqlim,:),2),1)); %small errors
    pwr_top(:,i_elect,2) = squeeze(mean(mean(n_pwr{1,i_elect}(:,freqlim,:),2),1)); %large errors
    pwr_top(:,i_elect,3) = squeeze(mean(mean(c_pwr{1,i_elect}(:,freqlim,:),2),1)); %no target
end
clear ii i_elect

% Get difference in mean power at frequency band for each electrode
pwr_top_diff = NaN([length(times),length(elect_topplot),3]); %pre-allocate
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes

    pwr_top_diff(:,i_elect,1) = pwr_top(:,i_elect,1)-pwr_top(:,i_elect,2); %small-large
    pwr_top_diff(:,i_elect,2) = pwr_top(:,i_elect,1)-pwr_top(:,i_elect,3); %small-no target
    pwr_top_diff(:,i_elect,3) = pwr_top(:,i_elect,2)-pwr_top(:,i_elect,3); %large-no target
end
clear ii i_elect


CLims = [-0.8 0.8]; %set power scale of plot
colormap('jet')

nconds = 3; %number of plots
% conds = {'Accurate';'Guess';'No Target'}; %labels for plots
conds = {'Accurate-Guesses';'Accurate-No Target';'Guesses-No Target'}; %labels for plots

for tw_i = 1:length(tWin) %loop through several time windows 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %finds the times you want from the times variable
    time_window = find(times>= itWin(1),1):find(times>= itWin(2),1)-1;
    
    figure('Color',[1 1 1],'Position',[1 1 941 349]); %need extra big figure

    for i_cond = 1:nconds %loop through conditions to make plot of each        
        subtightplot(1,3,i_cond,[0.02,0.02],[0.05,0.07],[0.05,0.05]);
        set(gca,'Color',[1 1 1]);
        
        temp = mean(pwr_top_diff(time_window,:,i_cond),1)'; %mean power within time window
        temp(1) = NaN; %so M2 is not included

        topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
        'plotchans',elect_topplot,'emarker',{'.','k',11,1})
        
        title(conds{i_cond});
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Power (dB)');
        clear temp
    end
    
    % Overall subplot title
    supertitle([num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz: ' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'],...
        'FontSize',10.5)

    clear itWin time_window i_cond 
end
clear tw_i nconds conds

clear freqlim CLims pwr_top pwr_top_diff





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
% freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [10 14]; %high alpha
freqband = [3 8]; %theta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
% timewin = [-500 -300];
timewin = [-600 -400];
% timewin = [-400 -200];
% timewin = [-200 0];
% timewin = [-100 0];

% timewin = [0 200];
% timewin = [200 400];
% timewin = [400 600];
% timewin = [300 600];
timelim = find(times>=timewin(1) & times<=timewin(2));
% -------------------------------------------------------------------------
% Calculate power in time window at freq band
current_power = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:5 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
%         part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        part_ersp = all_erspN{i_part,i_elect}.trials; %get single subject's ersp
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
    for ii = 1:5 %only central electrodes
%     for ii = 1:length(exp.singletrialselecs)
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
for i_part = 1:length(exp.participants)
%     for ii = 1:length(exp.singletrialselecs)
    for ii = 1:5 %only central electrodes
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % Get mean error
        mean_errdeg_Hpwr(i_part,i_elect) = mean(abs(errdeg_Hpwr{i_part,i_elect}));
        mean_errdeg_Lpwr(i_part,i_elect) = mean(abs(errdeg_Lpwr{i_part,i_elect}));
        % Get SD error
        std_errdeg_Hpwr(i_part,i_elect) = std(abs(errdeg_Hpwr{i_part,i_elect}));
        std_errdeg_Lpwr(i_part,i_elect) = std(abs(errdeg_Lpwr{i_part,i_elect}));
        
%        model_out_Hpwr{i_part,i_elect} = MemFit(errdeg_Hpwr{i_part,i_elect});
%        model_out_Lpwr{i_part,i_elect} = MemFit(errdeg_Lpwr{i_part,i_elect});
        model_out_Hpwr{i_part,i_elect} = MLE(errdeg_Hpwr{i_part,i_elect}); %fits without plotting
        model_out_Lpwr{i_part,i_elect} = MLE(errdeg_Lpwr{i_part,i_elect}); %fits without plotting
        
        % Save specific outputs
        g_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(1);
        g_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(1);
        sd_out_Hpwr(i_part,i_elect) = model_out_Hpwr{i_part,i_elect}(2);
        sd_out_Lpwr(i_part,i_elect) = model_out_Lpwr{i_part,i_elect}(2);
        
        clear i_elect
    end
    clear ii
end
clear i_part

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Run nonparametric statistics
% for ii = 1:length(exp.singletrialselecs)
for ii = 1:5 %only central electrodes
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    sgnrank_pwr.elect{1,ii} = exp.singtrlelec_name{ii}; %save name
    % nanmedian(g_out_Hpwr(:,i_elect))
    % nanmedian(g_out_Lpwr(:,i_elect))
    % 
    % nanmedian(sd_out_Hpwr(:,i_elect))
    % nanmedian(sd_out_Lpwr(:,i_elect))

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
    clear h p stat i_elect
end
clear ii

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Run parametric statistics
% for ii = 1:length(exp.singletrialselecs)
for ii = 1:5 %only central electrodes
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    ttest_pwr.elect{1,ii} = exp.singtrlelec_name{ii}; %save name
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
    clear h p ci stat i_elect
end
clear ii

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------









%##########################################################################
%##########################################################################
%--------------------------------------------------------------------------
%% <<<<<<<<<<<<<<<<<<<<<<<<<< PHASE ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%--------------------------------------------------------------------------
%##########################################################################
%##########################################################################


% exp is a function that might get used below
exp2 = exp;
clear exp

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Separate phase trials by response errors

current_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
x_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
n_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
errlims = cell(1,length(exp2.participants));    %pre-allocate
for i_part = 1:length(exp2.participants)
   
    % Get upper and lower limits based on model fit
    errlims{i_part}(1) = -(model_out{1,i_part}(2)); %negative value
    errlims{i_part}(2) = model_out{1,i_part}(2);
    
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp(:,:,:) = all_ersp{i_part,i_elect}; %get single subject's ersp
        current_phase{i_part,i_elect}(:,:,:) = squeeze(angle(part_ersp(:,:,:)));

        
        % Get phase values separated by errors in trial
        % small errors
        x_phase{i_part,i_elect}(:,:,:) = current_phase{i_part,i_elect}(:,:,...
            [find((resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)))] );
        % larger errors
        n_phase{i_part,i_elect}(:,:,:) = current_phase{i_part,i_elect}(:,:,...
            [find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))] );
        
        clear part_ersp i_trial
    end
    clear ii
end
clear i_elect i_part


% .............
% Phase for catch trials
c_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp2.participants)
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp(:,:,:) = all_ersp_byC{i_part,i_elect}; %get single subject's ersp

        % Get phase values from catch trials
        c_phase{i_part,i_elect}(:,:,:) = squeeze(angle(part_ersp(:,:,:)));
        
        clear part_ersp
    end
    clear ii
end
clear i_elect i_part


% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------

% Plot ITCPz by subject
for i_part = 1:length(exp2.participants)
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % Avg phase across trials (freq x time x trial)
        tmp_x = (size(x_phase{i_part,i_elect},3)).*(squeeze(abs(mean(exp(1i.*x_phase{i_part,i_elect}),3))).^2);
        tmp_n = (size(n_phase{i_part,i_elect},3)).*(squeeze(abs(mean(exp(1i.*n_phase{i_part,i_elect}),3))).^2);
        tmp_c = (size(c_phase{i_part,i_elect},3)).*(squeeze(abs(mean(exp(1i.*c_phase{i_part,i_elect}),3))).^2);
        
        figure('Position', [1 1 1685 405]); %open a new figure
        
        % Plot Small Errors
        subplot(1,3,1)
        contourf(times,freqs,tmp_x,40,'linecolor','none')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freqency (Hz)'); xlabel('Time (ms)');
        colorbar
%         colormap jet
        title('Accurate');
        
        % Plot Large Errors
        subplot(1,3,2)
        contourf(times,freqs,tmp_n,40,'linecolor','none')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freqency (Hz)'); xlabel('Time (ms)');
        colorbar
        colormap jet
        title('Guesses');
        
        % Plot Large Errors
        subplot(1,3,3)
        contourf(times,freqs,tmp_c,40,'linecolor','none')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freqency (Hz)'); xlabel('Time (ms)');
        colorbar
        colormap jet
        title('No Target');
        
        % Overall subplot title
        supertitle(['ITPCz: Subj ' num2str(exp2.participants{i_part}) ': ' exp2.singtrlelec_name{ii}],'FontSize',10.5)
      
       clear i_elect tmp_x tmp_n tmp_c  
    end
    clear ii
end
clear i_part


% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------

% Plot POS p-value by subject
for i_part = 1:length(exp2.participants)
    for ii = 1
%     for ii = 1:length(exp2.singletrialselecs)    
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % Get POS p-vals at each freq x time across trial
        [p_xn.circWW{i_part,i_elect}, p_xn.POS{i_part,i_elect}, p_xn.zPOS{i_part,i_elect}] =...
            PhaseOpposition(x_phase{i_part,i_elect},n_phase{i_part,i_elect});
        [p_cn.circWW{i_part,i_elect}, p_cn.POS{i_part,i_elect}, p_cn.zPOS{i_part,i_elect}] =...
            PhaseOpposition(c_phase{i_part,i_elect},n_phase{i_part,i_elect});
        [p_cx.circWW{i_part,i_elect}, p_cx.POS{i_part,i_elect}, p_cx.zPOS{i_part,i_elect}] =...
            PhaseOpposition(c_phase{i_part,i_elect},x_phase{i_part,i_elect});
        
        figure('Position', [1 1 1685 405]); %open a new figure
        CLim = [0 0.055];
        
        % Plot Accurate vs Guesses
        subplot(1,3,1)
        imagesc(times,freqs,p_xn.circWW{i_part,i_elect},CLim)
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freqency (Hz)'); xlabel('Time (ms)');
        colorbar
        title('Accurate vs Guesses');
        % Plot No Target vs Guesses
        subplot(1,3,2)
        imagesc(times,freqs,p_cn.circWW{i_part,i_elect},CLim)
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freqency (Hz)'); xlabel('Time (ms)');
        colorbar
        title('No Target vs Guesses');
        % Plot No Target vs Accuracy
        subplot(1,3,3)
        imagesc(times,freqs,p_cx.circWW{i_part,i_elect},CLim)
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
        ylim([3 35]); yticks(5:5:35)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freqency (Hz)'); xlabel('Time (ms)');
        colorbar
        title('No Target vs Accuracy');
        
        % Overall subplot title
        supertitle(['POS p-value: Subj ' num2str(exp2.participants{i_part}) ': ' exp2.singtrlelec_name{ii}],'FontSize',10.5)
      
        
       clear i_elect
    end
    clear ii
end
clear i_part









% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%% >>>>>>>>>>>>>>>>>>>>  PHASE ANALYSES  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@






















