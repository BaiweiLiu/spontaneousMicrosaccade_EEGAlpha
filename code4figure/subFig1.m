%% Analysis script for search task in mutil task project

%% It's always good to start with a clean sheet
clc, clear, close all, warning('off','all')

%% set project name
projectname = 'SF figure 1';

%% prepare Toolbox
parent_folder = '/Users/baiweil/AnalysisDesk/Spontaneous_alpha/EEG_eye_delay_NN/';
package_folder = '/Users/baiweil/Documents/MATLAB_pacakges/';

cfg= [];
cfg.package_folder = package_folder;
cfg.futureTB_path = [package_folder 'futureToolBox'];
cd(cfg.futureTB_path);
futureTB_setup(cfg);

%% Set directions
write_dir_eeg = [parent_folder 'eeg_results'];
read_dir_reepoch = [write_dir_eeg filesep 're_epoch/corr' ];

%% get GA data 
sublist = get_subFiles(read_dir_reepoch);

for subjInd = 1 :length(sublist)
    load(sublist{subjInd})
    GA.saccEvent{subjInd} = eeg_data.saccEvent;
    GA.gazeVel{subjInd} =eeg_data.gazeVel;
    GA.gazePos{subjInd} =eeg_data.gazePos;
    GA.time4gazeVel{subjInd} =eeg_data.time4gazeVel;
    GA.OldtrialNumber{subjInd} =eeg_data.OldtrialNumber;
    GA.timetoEncoding{subjInd} =eeg_data.timetoEncoding;
    GA.medianGaze{subjInd} =eeg_data.medianGaze;
    GA.shiftStartEnd{subjInd} =eeg_data.shiftStartEnd;
    GA.toward{subjInd} =eeg_data.toward;
    GA.dwell_time{subjInd} =eeg_data.dwell_time;
    GA.start_return{subjInd} =eeg_data.start_return;
    GA.thre_gazeVol{subjInd} =eeg_data.thre_gazeVol;
    GA.displace{subjInd} =eeg_data.displace;
    GA.time = eeg_data.time;
end

%% 
GA_struct = GA;
save([write_dir_eeg filesep 'group_results' filesep 'GA_supFigure1.mat'], 'GA_struct')

%% 
load([write_dir_eeg filesep 'group_results' filesep 'GA_supFigure1.mat'])

%% sub_fig1a

cmap = flipud(brewermap([],'*PRGn'));
c1 = cmap(end-40,:);
c2 = cmap(1+40,:);
mFig = true;
data2plot_sum = [];
num_left_sum = [];
num_right_sum = [];

figure('position', [100, 100, 1200, 800], 'color', [1 1 1]) 
for subjInd = 1:length(GA.saccEvent)
    subplot( 5, 5, subjInd)
    data2plot = GA_struct.displace{subjInd};
    data2plot = data2plot * (5.7/100);
    
    data2plot_usable = data2plot(abs(data2plot) > 0.057 & abs(data2plot) < 1);
    num_left = sum(data2plot_usable<0);
    num_right = sum(data2plot_usable>0);
    h1 = histogram(data2plot,[-1: 0.05: -0.1],'Normalization','probability'); hold on
    h2 = histogram(data2plot,[0.1: 0.05: 1],'Normalization','probability'); hold on
    
    h3 = histogram(data2plot,[-3: 0.05: -1],'Normalization','probability'); hold on
    h4 = histogram(data2plot,[1: 0.05: 3],'Normalization','probability'); hold on
    
    plot([1 1], [0 0.15], '--k');
    plot([-1 -1], [0 0.15], '--k');
    h1.FaceColor = c1;
    h1.EdgeColor = c1;
    
    
    h2.FaceColor = c2;
    h2.EdgeColor = c2;
    
    h3.FaceColor = [0.6 0.6 0.6] + 0.1;
    h3.EdgeColor = [0.6 0.6 0.6]+ 0.1;

    h4.FaceColor = [0.6 0.6 0.6] + 0.1;
    h4.EdgeColor = [0.6 0.6 0.6] + 0.1;
    
    set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2,'TickDir','out');
    box off
    xlim([-3 3])
    ax = gca; ax.TickLength = [0.015 0.03]; ax.XTick = [-3:1:3]; ax.XTickLabel = [-3:1:3];
    ax.YTick= [0 0.1]; ylim([0 0.12])

    if mFig; xticklabels({}); yticklabels({}); end
    data2plot_sum = [data2plot_sum; data2plot];
    
    title([num2str(num_left) ' ' num2str(num_right)],'FontSize', 20)
    num_left_sum  = [num_left_sum num_left];
    num_right_sum = [num_right_sum num_right];
end

subplot( 5, 5, 25)

h1 = histogram(data2plot_sum,[-1: 0.05: -0.1],'Normalization','probability'); hold on
h2 = histogram(data2plot_sum,[0.1: 0.05: 1],'Normalization','probability'); hold on

h3 = histogram(data2plot_sum,[-3: 0.05: -1],'Normalization','probability'); hold on
h4 = histogram(data2plot_sum,[1: 0.05: 3],'Normalization','probability'); hold on
    plot([1 1], [0 0.08], '--k');
    plot([-1 -1], [0 0.08], '--k');
h1.FaceColor = c1;
h1.EdgeColor = c1;

h2.FaceColor = c2;
h2.EdgeColor = c2;

h3.FaceColor = [0.6 0.6 0.6] + 0.1;
h3.EdgeColor = [0.6 0.6 0.6]+ 0.1;

h4.FaceColor = [0.6 0.6 0.6] + 0.1;
h4.EdgeColor = [0.6 0.6 0.6] + 0.1;

set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2,'TickDir','out');
box off
xlim([-3 3]);
ax = gca; ax.TickLength = [0.015 0.03]; ax.XTick = [-3:1:3]; ax.XTickLabel = [-3:1:3];
ax.YTick= [0 0.05]; ylim([0 0.06]);

if mFig; xticklabels({}); yticklabels({}); end
num_left_sum = sum(num_left_sum);
num_right_sum = sum(num_right_sum);
title([num2str(num_left_sum) ' ' num2str(num_right_sum)],'FontSize', 20)
set(gcf, 'renderer', 'Painters'); 
saveas(gcf, ['/Users/baiweil/Library/Mobile Documents/com~apple~CloudDocs/Projects/Current Projects/ReAnalysis/EEG_eye_delay_NN/figs_ms/sup_fig1a.eps'], 'epsc');

%% sub_fig1b
cmap = flipud(brewermap([],'*PRGn'));
c1 = cmap(end-40,:);
c2 = cmap(1+40,:);
mFig = true;
data2plot_sum = [];
num_trial_sum = [];
figure('position', [100, 100, 1200, 800], 'color', [1 1 1]) 
for subjInd = 1:length(GA.saccEvent)
    subplot( 5, 5, subjInd)
    
    data_size = GA_struct.displace{subjInd};
    data_size = data_size * (5.7/100);
    
    sel_usable = abs(data_size) > 0.057 & abs(data_size) < 1;
    
    data2plot = GA_struct.timetoEncoding{subjInd}(sel_usable);
    h1 = histogram(data2plot,[0: 50: 2000],'Normalization','probability'); hold on
     
     plot([500 500], [0 0.15], '--k');
     plot([1500 1500], [0 0.15], '--k');

     h1.FaceColor = [0.5 0.5 0.5] ;
     h1.EdgeColor = [0.5 0.5 0.5];

    set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2,'TickDir','out');
    box off
     xlim([0 2000])
     ax = gca; ax.TickLength = [0.015 0.03]; ax.XTick = [0 :500 :2000]; ax.XTickLabel = [0:500:2000];
%     ax.YTick= [0 0.1]; 
    ylim([0 0.1])

    if mFig; xticklabels({}); yticklabels({}); end
    data2plot_sum = [data2plot_sum; data2plot];
    
    title([num2str(length(data2plot))],'FontSize', 20)
    
    num_trial_sum = [num_trial_sum length(data2plot)];
end

subplot( 5, 5, 25)

data2plot = GA_struct.timetoEncoding{subjInd};
h1 = histogram(data2plot_sum,[0: 50: 2000],'Normalization','probability'); hold on

plot([500 500], [0 0.15], '--k');
plot([1500 1500], [0 0.15], '--k');

h1.FaceColor = [0.5 0.5 0.5];
h1.EdgeColor = [0.5 0.5 0.5];

set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2,'TickDir','out');
box off
xlim([0 2000])
ax = gca; ax.TickLength = [0.015 0.03]; ax.XTick = [0 :500 :2000]; ax.XTickLabel = [0:500:2000];
%     ax.YTick= [0 0.1]; 
ylim([0 0.1])
if mFig; xticklabels({}); yticklabels({}); end

num_trial_sum = sum(num_trial_sum);

title([num2str(num_trial_sum)],'FontSize', 20)
set(gcf, 'renderer', 'Painters'); 
saveas(gcf, ['/Users/baiweil/Library/Mobile Documents/com~apple~CloudDocs/Projects/Current Projects/ReAnalysis/EEG_eye_delay_NN/figs_ms/sup_fig1b.eps'], 'epsc');