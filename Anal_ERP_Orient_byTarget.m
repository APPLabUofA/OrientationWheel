


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
%%                  $$$$$$$ BEH Data $$$$$$$

% Get BEH data excluding trials that were rejected when creating epochs
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



% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%%                  ERPs by Errors - Mean
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

errlims = cell(1,length(exp.participants));    %pre-allocate
errs_x = cell(1,length(exp.participants));    %pre-allocate
errs_n = cell(1,length(exp.participants));    %pre-allocate
% create separate ERPs for large and small errors
for i_part = 1:length(exp.participants) % --------------------------------- 
    
    % Get 25th and 75th quantile values (can use because data is normally
    % distributed)
%     errlims{i_part} = quantile(resp_errdeg{i_part},[0.25 0.75]);
    
    % Get upper and lower limits based on model fit
%     errlims{i_part}(1) = -(model_out{1,i_part}.maxPosterior(2)); %negative value
%     errlims{i_part}(2) = model_out{1,i_part}.maxPosterior(2);
    errlims{i_part}(1) = -(model_out{1,i_part}(2)); %negative value
    errlims{i_part}(2) = model_out{1,i_part}(2);
    
    % Calculate ERP
    for ii = 1:length(elect_erp)
        i_elect = elect_erp(ii); %for doing only a selection of electrodes
        
        % Get trials with small errors
%         erp_out_x(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,...
%             [find((resp_errdeg{i_part}<errlims{i_part}(2) & resp_errdeg{i_part}>errlims{i_part}(1)))] ),3));
        erp_out_x(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,...
            [find((resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)))] ),3));
        %save small errors
%         errs_x{i_part} = resp_errdeg{i_part}(resp_errdeg{i_part}<errlims{i_part}(2) & resp_errdeg{i_part}>errlims{i_part}(1));
        errs_x{i_part} = resp_errdeg{i_part}(resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75));
        
        % Get trials with large errors
%         erp_out_n(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,...
%             [find(resp_errdeg{i_part}>=errlims{i_part}(2)) find(resp_errdeg{i_part}<=errlims{i_part}(1))] ),3));
        erp_out_n(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,...
            [find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))] ),3));
        %save large errors
%         errs_n{i_part} = resp_errdeg{i_part}([find(resp_errdeg{i_part}>=errlims{i_part}(2)) find(resp_errdeg{i_part}<=errlims{i_part}(1))]);
        errs_n{i_part} = resp_errdeg{i_part}([find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))]);
    end
    clear ii i_elect
end
clear i_part

% /////////////////////////////////////////////////////////////////////////



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% ::::::::::::::::  Plot single subj ERPs by electrode  :::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% for i_part = 1:length(exp.participants)
for i_part = 22
    % Get data by subject into a temporary variable for plotting
    erp_out_temp(:,:,1) = squeeze(erp_out_x(i_part,:,:,:)); %small errors
    erp_out_temp(:,:,2) = squeeze(erp_out_n(i_part,:,:,:)); %large errors
    erp_out_temp(:,:,3) = squeeze(erp_out_x(i_part,:,:,:) - erp_out_n(i_part,:,:,:)); %difference
    
    figure('Position', [1 1 624 1016]) %open a new figure
    
    % Plot by small and large error trials & difference ERP
%     for ii = 1:length(elect_erp)
    for ii = 1:5 %only the central electrodes   
        
    %     i_elect = elect_erp(ii); %for doing only a selection of electrodes
        i_elect = ii; %selection done when making ERPs
        
        % Get axes limits
        ymin = floor(min([erp_out_temp(i_elect,:,1) erp_out_temp(i_elect,:,2)]));
        ymax = ceil(max([erp_out_temp(i_elect,:,1) erp_out_temp(i_elect,:,2)]));
%         ymin = -3; ymax = 12;
        xmin = -300; xmax = 600;

        subplot(5,1,ii) %for plotting only the central electrodes
        plot(EEG.times,erp_out_temp(i_elect,:,1),'color',[0.14 0.93 0.9],'LineWidth',1.5)
        hold on
        plot(EEG.times,erp_out_temp(i_elect,:,2),'-.m','LineWidth',1.5)
        hold on
%         plot(EEG.times,erp_out_temp(i_elect,:,3),'color',[0 1 0.5],'LineWidth',1)
%         hold on
        line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
        line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
        line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
        line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
        set(gca,'ydir','reverse'); xlim([xmin xmax]); ylim([ymin ymax])

        title(el_erp_names{ii}); 
        xlabel('Time (ms)'); ylabel('Voltage (uV)'); xticks(xmin:200:xmax);
        if ii == 3 %add a legend only to the bottom plot
%             legend({'Small','Large','Diff: S-L'},'Location','best');
            legend({'Small','Large'},'Location','best');
        end
        
        clear ymin ymax
    end
    % Overall subplot title
    supertitle(['Subj ' num2str(exp.participants{i_part}) ': ERPs by Response Error'],...
        'FontSize',10.5)
    
    clear ii xmax xmin erp_out_temp
end
clear i_part

% /////////////////////////////////////////////////////////////////////////



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% ::::::::::::::::::  Plot the ERPs by electrode  :::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% ------------------------------------------------------------------------- 
% average across subjects
erp_out_byerr(:,:,1) = squeeze(mean(erp_out_x(:,:,:),1)); %small errors
erp_out_byerr(:,:,2) = squeeze(mean(erp_out_n(:,:,:),1)); %large errors
erp_out_byerr(:,:,3) = squeeze(mean((erp_out_x(:,:,:)-erp_out_n(:,:,:)),1)); %difference
% ------------------------------------------------------------------------- 


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Plot by small and large error trials & difference ERP
% for ii = 1:8
for ii = 1:length(elect_erp)    
    i_elect = elect_erp(ii); %for doing only a selection of electrodes
%     i_elect = ii; %selection done when making ERPs
    % get axes limits
%     ymin = floor(min([erp_out_byerr(i_elect,:,1) erp_out_byerr(i_elect,:,2)]));
%     ymax = ceil(max([erp_out_byerr(i_elect,:,1) erp_out_byerr(i_elect,:,2)]));
    ymin = -5; ymax = 9;
    xmin = -400; xmax = 800;
    
    figure
    %small errors
    plot(EEG.times,erp_out_byerr(i_elect,:,1),'color',[0.14 0.93 0.9],'LineWidth',1.5)
    hold on
    %large errors
    plot(EEG.times,erp_out_byerr(i_elect,:,2),'-.m','LineWidth',1.5)
    hold on
    %difference
    plot(EEG.times,erp_out_byerr(i_elect,:,3),'color',[0 1 0.5],'LineWidth',1)
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
    line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    set(gca,'ydir','reverse'); xlim([xmin xmax]); ylim([ymin ymax])
    
    title([el_erp_names{ii} ': ERPs by Response Error']); 
    xlabel('Time (ms)'); ylabel('Voltage (uV)'); xticks(xmin:200:xmax);
    legend({'Small','Large','Diff: S-L'},'Location','best');
    clear ymin ymax
end
clear ii xmax xmin


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Plot ERPs with error bars
% for ii = 1:length(elect_erp) 
for ii = 1:5 
    i_elect = elect_erp(ii); %for doing only a selection of electrodes
%     i_elect = ii; %selection done when making ERPs
    % get axes limits
    ymin = -6; ymax = 10;
    xmin = -200; xmax = 800;
    
    figure('Color',[1 1 1]); 
    boundedline(EEG.times,erp_out_byerr(i_elect,:,1),squeeze(std(erp_out_x(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'c',...
            EEG.times,erp_out_byerr(i_elect,:,2),squeeze(std(erp_out_n(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'m');
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    set(gca,'ydir','reverse'); xlim([xmin xmax]); ylim([ymin ymax])
    
    title([el_erp_names{ii} ': ERPs by Response Error']); 
    xlabel('Time (ms)'); ylabel('Voltage (uV)'); 
    xticks(xmin:200:xmax); yticks(ymin:2:ymax)    
        
    legend({'Accurate','Guess'},'Location','best');
end
clear ii xmax xmin    


% /////////////////////////////////////////////////////////////////////////
%% ::::::::::::::::::::  Statistics: t-Test  ::::::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%%% Pick your time window 300-500 for P3, 150-250 for MMN %%%
% time1 = 300;0
% time2 = 500;

% time1 = 200;
% time2 = 300;

% time1 = 150;
% time2 = 250;

time1 = 0;
time2 = 50;

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-1;

for ii = 1:length(elect_erp)    
    i_elect = elect_erp(ii); %for doing only a selection of electrodes

    [h,p,ci,stat] = ttest(mean(erp_out_x(:,i_elect,time_window),3),mean(erp_out_n(:,i_elect,time_window),3),.05,'both',1);

    ttest_electrodes(ii,1) = h;
    ttest_electrodes(ii,2) = p;
    ttest_electrodes_ci(ii,1) = ci(1);
    ttest_electrodes_ci(ii,2) = ci(2);
    ttest_electrodes_stats(ii,1).electrode = stat;
    ttest_electrodes_id{ii,1} = el_erp_names{ii};

    clear i_elect h p ci stat
end
clear ii


% /////////////////////////////////////////////////////////////////////////
%% ::::::::::::::::::  Statistics: Correlation  :::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%%% Pick your time window 350-550 for P3 %%%
% time1 = 350;
% time2 = 550;

time1 = 200;
time2 = 300;

% time1 = 100;
% time2 = 200;

% time1 = 150;
% time2 = 250;


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


% /////////////////////////////////////////////////////////////////////////
% Correlation across all subjects
% /////////////////////////////////////////////////////////////////////////
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Put all the response errors across subjects into vector
resp_errdeg_cat = cat(2,resp_errdeg{1:end});
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Test correlations and polar scatter plot
% for ii = 1:length(elect_erp)  
for ii = 2:6
    i_elect = elect_erp(ii); %for doing only a selection of electrodes
    
    % Put all the amplitudes across subjects into vector
    tmp_amp = cat(1,amp_out_timewin{1:end,i_elect});
    % Actual correlation test
    [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg_cat),tmp_amp');
    
    % correlation betwen power and errors and plot
    figure; polarscatter(circ_ang2rad(resp_errdeg_cat),tmp_amp','filled','MarkerFaceAlpha',.5)
    rlim([0 60])
%     figure; scatter(resp_errdeg_cat,tmp_amp')
    title([el_erp_names{ii} ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,4)) '; ' num2str(time1) '-' num2str(time2) ' ms'])
    
    erp_corr.p(i_elect) = pval;
    
    clear rho pval tmp_amp i_elect
end
clear ii



% /////////////////////////////////////////////////////////////////////////
% Correlation for each subject
% /////////////////////////////////////////////////////////////////////////

% Test correlations and polar scatter plot
for i_part = 1:length(exp.participants)
    % for ii = 1:length(elect_erp)  
    for ii = 2
        i_elect = elect_erp(ii); %for doing only a selection of electrodes

        % Actual correlation test
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),amp_out_timewin{i_part,i_elect}');

        % Polar scatterplot with transparent circles
        figure; polarscatter(circ_ang2rad(resp_errdeg{i_part}),amp_out_timewin{i_part,i_elect}','filled','MarkerFaceAlpha',.5)
        rlim([0 60])
    %     figure; scatter(resp_errdeg_cat,tmp_amp')
        title(['Subj ' num2str(exp.participants{i_part}) ' ' el_erp_names{ii} ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,4)) '; ' num2str(time1) '-' num2str(time2) ' ms'])

        erp_corr.p(i_elect) = pval;

        clear rho pval tmp_amp i_elect
    end
    clear ii
end
clear i_part


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     '''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

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

tWin{1} = [300 500];

CLims = [-8 8];
nconds = 3; %plotting the large, small, and difference ERPs
conds = {'Accurate';'Guess';'Accurate - Guess'}; %labels for plots

for tw_i = 1:length(tWin) %loop through several time windows 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %this code finds the times you want from the timess variable
    time_window = find(EEG.times>= itWin(1),1):find(EEG.times>= itWin(2),1)-1;
    
    figure('Color',[1 1 1],'Position',[1 1 941 349]);

    for i_cond = 1:nconds
        subtightplot(1,3,i_cond,[0.02,0.02],[0.05,0.07],[0.05,0.05]);
        set(gca,'Color',[1 1 1]);
        temp = mean(erp_out_byerr(:,time_window,i_cond),2)'; %ERP within time window
        temp(1) = NaN;
        if i_cond == 3
            CLims = [-4 4];
            topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
                'plotchans',elect_erp,'emarker',{'.','k',11,1})
    %         topoplot(temp,EEG.chanlocs,'whitebk','on',0.6,'maplimits',...
    %             'plotchans',elect_erp,'emarker',{'.','k',11,1})
        else
            CLims = [-8 8];
            topoplot(temp,ALLEEG(1).chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
            'plotchans',elect_erp,'emarker',{'.','k',11,1})
%         topoplot(temp,EEG.chanlocs,'whitebk','on',0.6,'maplimits',...
%             'plotchans',elect_erp,'emarker',{'.','k',11,1})
        end
        title(conds{i_cond});
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Voltage (uV)');
        clear temp
    end
    
    % Overall subplot title
    supertitle([num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'],...
        'FontSize',10.5)
    
    clear itWin time_window i_cond

end
clear tw_i nconds conds tWin






































% /////////////////////////////////////////////////////////////////////////

clear erp_out_byerr erp_out_x erp_out_n

















