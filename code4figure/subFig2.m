%% sup figure 2 plot 
load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_cvsi_TF_O12.mat')
GA_struct_cvsi= GA_struct;

load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_TF.mat')
GA_struct_tf = GA_struct;

load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_cvsi_TF.mat')
GA_struct_cvsi_figure1= GA_struct;

load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_TF.mat')
GA_struct_tf_figure1 = GA_struct;
%% 
mFig = true;

time_i_lamda = [0.05 0.2];
time_i_alpha = [0.05 0.2];
figure('position', [100, 100, 1200, 1050], 'color', [1 1 1]) 
subplot(3,4,[1:2] +4)
GA_struct = GA_struct_cvsi;
timeLim = [-0.5 0.5];
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(1:40);
cfg.tf_data = squeeze(GA_struct.power_cvsi(:,1,1,1,1:40,:));
cfg.zli = [-5 5];
cfg.cluster_mask = false;
cfg.npermutations = 100;
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
plot(timeLim, [8 8], ':k');
plot(timeLim, [12 12], ':k');

plot([time_i_lamda(1) time_i_lamda(2)], [3 3], 'k','LineWidth',3);
plot([time_i_lamda(1) time_i_lamda(2)], [7 7], 'k','LineWidth',3);
plot([time_i_lamda(1) time_i_lamda(1)], [3 7], 'k','LineWidth',3);
plot([time_i_lamda(2) time_i_lamda(2)], [3 7], 'k','LineWidth',3);

set(gca,'LineWidth',3)
    set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
    set(gca,'TickDir','out');

if mFig; xticklabels({}); yticklabels({}) 
else 
    colorbar('eastoutside')	
    ylabel('Freq')    
    xlabel('time from onset of saccade')
end
ax = gca;
ax.XTick = [-0.5 -0.25 0 0.25 0.5];

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.2);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on

subplot(3,4,[3:4] +4)
GA_struct = GA_struct_cvsi_figure1;
timeLim = [-0.5 0.5];
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(1:40);
cfg.tf_data = squeeze(GA_struct.power_cvsi(:,1,1,1,1:40,:));
cfg.zli = [-8 8];
cfg.npermutations = 100;
cfg.cluster_mask = false;
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
plot(timeLim, [8 8], ':k');
plot(timeLim, [12 12], ':k');

plot([time_i_alpha(1) time_i_alpha(2)], [8 8], 'k','LineWidth',3);
plot([time_i_alpha(1) time_i_alpha(2)], [12 12], 'k','LineWidth',3);
plot([time_i_alpha(1) time_i_alpha(1)], [8 12], 'k','LineWidth',3);
plot([time_i_alpha(2) time_i_alpha(2)], [8 12], 'k','LineWidth',3);

set(gca,'LineWidth',3)
    set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
    set(gca,'TickDir','out');

if mFig; xticklabels({}); yticklabels({}) 
else 
    colorbar('eastoutside')	
    ylabel('Freq')    
    xlabel('time from onset of saccade')
end
ax = gca;
ax.XTick = [-0.5 -0.25 0 0.25 0.5];

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.2);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on

GA_struct = GA_struct_tf;

cond = 'MStype' % MStype or toward 

subplot(3,4,6+4)
freq_i = [3:7];

for condInd = 1:1
    %data2plot1 = squeeze(mean(GA_struct.diff_pow(:,1,2,1,:,:,:)));
    %data2plot2 = squeeze(mean(GA_struct.diff_pow(:,1,3,1,:,:,:)));

    dummy = [];
    dummy.time = GA_struct.time;
    dummy.freq = GA_struct.freq;
    dummy.label = GA_struct.eleLabel;
    dummy.dimord = 'chan_freq_time';
    if strcmp(cond, 'MStype')
        dummy.data2plot1 = squeeze(mean(GA_struct.diff_pow(:,1,condInd,1,:,:,:))); % mean over trials (reduce to 3 dimensions: electrodes, frequency points, time points (on which we centred our sliding time windows)
    elseif strcmp(cond, 'toward')
        dummy.data2plot1 = squeeze(mean(GA_struct.diff_pow(:,1,1,condInd,:,:,:)));
    end
    % now lets plot it

    cmap = brewermap([],'*RdBu');
    colormap(cmap)
    cfg = [];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.marker = 'off';;
     cfg.xlim = time_i_lamda;
    cfg.figure = 'gcf';
    cfg.colormap = cmap;
    
    cfg.parameter = 'data2plot1';
    cfg.ylim = freq_i;
    cfg.zlim =  [-10,10];
    ft_topoplotTFR(cfg, dummy); %title([num2str(timeStart(timeInd)) ' to ' num2str(timeStart(timeInd)+ timewin)])
end

subplot(3,4,8+4)
freq_i = [8:15];

for condInd = 1:1
    %data2plot1 = squeeze(mean(GA_struct.diff_pow(:,1,2,1,:,:,:)));
    %data2plot2 = squeeze(mean(GA_struct.diff_pow(:,1,3,1,:,:,:)));

    dummy = [];
    dummy.time = GA_struct.time;
    dummy.freq = GA_struct.freq;
    dummy.label = GA_struct.eleLabel;
    dummy.dimord = 'chan_freq_time';
    if strcmp(cond, 'MStype')
        dummy.data2plot1 = squeeze(mean(GA_struct.diff_pow(:,1,condInd,1,:,:,:))); % mean over trials (reduce to 3 dimensions: electrodes, frequency points, time points (on which we centred our sliding time windows)
    elseif strcmp(cond, 'toward')
        dummy.data2plot1 = squeeze(mean(GA_struct.diff_pow(:,1,1,condInd,:,:,:)));
    end
    % now lets plot it

    cmap = brewermap([],'*RdBu');
    colormap(cmap)
    cfg = [];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.marker = 'off';;
     cfg.xlim = time_i_alpha;
    cfg.figure = 'gcf';
    cfg.colormap = cmap;
    
    cfg.parameter = 'data2plot1';
    cfg.ylim = freq_i;
    cfg.zlim =  [-20 20];
    ft_topoplotTFR(cfg, dummy); %title([num2str(timeStart(timeInd)) ' to ' num2str(timeStart(timeInd)+ timewin)])
end

%
choiL = ismember(GA_struct.eleLabel, {'O1' 'O2' }); % choi = channel-of-interest
choiR = ismember(GA_struct.eleLabel, {'PO7' 'PO8'}); % choi = channel-of-interest

% plot electrodes on an empty topography
cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.style = 'blank';
cfg.comment = 'no';
cfg.highlight = 'electrode';
cfg.parameter = 'data2plot1';

subplot(3,4,2); cfg.highlightchannel = find(choiL); ft_topoplotER(cfg, dummy);

subplot(3,4,4); cfg.highlightchannel = find(choiR); ft_topoplotER(cfg, dummy);
%%
set(gcf, 'renderer', 'Painters'); 
saveas(gcf, ['/Users/baiweil/Library/Mobile Documents/com~apple~CloudDocs/Projects/Current Projects/ReAnalysis/EEG_eye_delay_NN/figs_ms/sup_fig2.eps'], 'epsc');
