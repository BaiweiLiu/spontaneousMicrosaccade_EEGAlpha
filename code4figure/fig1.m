%% figure 1 plot 
load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_cvsi_TF.mat')
GA_struct_cvsi= GA_struct;
load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_TF.mat')
GA_struct_tf = GA_struct;
%%
timeLim = [-0.5 0.5];
cmap = flipud(brewermap([],'*RdBu'));
mFig = true;

GA_struct = GA_struct_cvsi;

figure('position', [100, 100, 600, 800], 'color', [1 1 1])
subplot(11,4,1:4)
data2plot= squeeze(GA_struct.power_gazeVel(:,1,1,1,:)) * (5.7/100);
lineErr_plot(GA_struct.time4gazeVel,data2plot, 'k', 'ci')
%plot([0 0], [0 2], '--k')
xlim(timeLim)
ylim([0 0.1])
set(gca,'LineWidth',3)
set(gca, 'FontSize', 15,'Fontname','Arial', 'LineWidth', 1.2);
set(gca,'TickDir','out');

ax = gca;
ax.YTick = [0 0.1]; 

if mFig; xticklabels({}); yticklabels({}) 
else 
    colorbar('eastoutside')	
    ylabel('ALpha lateralisation')    
    xlabel('Time from onset of saccade')
end


subplot(11,4,5:20)
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.power_cvsi(:,1,1,1,3:40,:));
cfg.zli = [-8 8];
%cfg.cluster_mask = false;
cfg.npermutations =10000
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
plot(timeLim, [8 8], ':k');
plot(timeLim, [12 12], ':k');
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
% ax.XTick = [0 0.4 0.5 0.8 1];
% ax.TickLength = [0.025 0.03]

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.2);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on

% sub fig 2
subplot(11,4,21:36)
data2plot_alpha = squeeze(mean(GA_struct.power_cvsi(:,1,1,1,[8:12],:),5));
%data2plot_beta = squeeze(mean(GA_struct.power_cvsi(:,1,1,1,[15:25],:),5));
lineErr_plot(GA_struct.time, data2plot_alpha, cmap(end-50,:), 'ci');
%lineErr_plot(GA_struct.time, data2plot_beta, cmap(1+50,:), 'ci');
plot([GA_struct.time(1) GA_struct.time(end)], [0 0], '--k')
plot([0 0], [-12 3], '--k')
xlim(timeLim)
ylim([-12 3])
set(gca,'LineWidth',3)
set(gca, 'FontSize', 15,'Fontname','Arial', 'LineWidth', 1.2);
set(gca,'TickDir','out');

%
time2plot = GA_struct.time >= -0.5 & GA_struct.time <= 0.5;
compData = zeros(size(data2plot_alpha(:,time2plot)));
cfg = [];
cfg.xax = GA_struct.time(time2plot);
cfg.npermutations = 10000;
cfg.clusterStatEvalaluationAlpha= 0.05;
cfg.nsub=23;
cfg.statMethod = 'montecarlo';
cfg.time_i = [-0.5,0.5];
state_t = cluster_perm_1D(cfg,data2plot_alpha(:,time2plot),compData);

mask_xxx = double(state_t.mask); mask_xxx(mask_xxx==0) = nan;
plot(GA_struct.time(time2plot), mask_xxx * -11.5, 'k', 'LineWidth', 3);
ax = gca; ax.XTick = [-0.5 : 0.25: 0.5]; ax.YTick = [-12 : 3: 3]; 

if mFig; xticklabels({}); yticklabels({}) 
else 
    colorbar('eastoutside')	
    ylabel('ALpha lateralisation')    
    xlabel('Time from onset of saccade')
end

% fig 3 
GA_struct = GA_struct_tf;

cond = 'MStype' % MStype or toward 
time_i = [0 0.2];
freq_i = [8:12];

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
%     cfg.xlim = [0 0.2];
    cfg.figure = 'gcf';
    cfg.colormap = cmap;
    
     timewin = 0.25; timeWin = [-0.5 0.5]; timeStep = 0.25;
     timeStart = [timeWin(1):timeStep:(timeWin(2)-timewin)];
     for timeInd = 1:length(timeStart)
        
        subplot(11,4,[36 + timeInd, 40 + timeInd])
        cfg.parameter = 'data2plot1';

            cfg.xlim = [timeStart(timeInd),timeStart(timeInd)+ timewin];
            %cfg.xlim = [0,0.25];
            cfg.ylim = freq_i;
            cfg.zlim =  [-20 20];
            ft_topoplotTFR(cfg, dummy); %title([num2str(timeStart(timeInd)) ' to ' num2str(timeStart(timeInd)+ timewin)])
     end
end

%% 
set(gcf, 'renderer', 'Painters'); 
saveas(gcf, ['/Users/baiweil/Library/Mobile Documents/com~apple~CloudDocs/Projects/Current Projects/ReAnalysis/EEG_eye_delay_NN/figs_ms/fig1.eps'], 'epsc');

