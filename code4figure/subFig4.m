%% figure sup itpc for return  

load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_cvsi_TF-PO78.mat')
GA_struct_PO78= GA_struct;

% itpc PO78
GA_struct = GA_struct_PO78;
mFig = true
timeLim = [-0.5 0.5] ; 
zl = [0 0.25];
figure('position', [100 100 1200 500], 'color', [1 1 1])

subplot(2,6,[1:2])
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.phaseLock_co(:,3,3:40,:));
cfg.zli = zl;
cfg.cluster_mask = false;
cfg.npermutations =10
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.5);
set(gca,'TickDir','out');
xlim([-0.5 0.5])

if mFig;  yticklabels({}); xlabel({}); xticklabels({}) 
else 
    title('con')
    colorbar('eastoutside')	
    ylabel('ALpha lateralisation')    
    xlabel('Time from onset of saccade')
end
ax = gca;
ax.XTick = [-0.5 : 0.25: 0.5];
ax.YTick = [3 10 20 30 40];
ax.TickLength = [0.015 0.015];

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.5);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
    
    
subplot(2,6,[3:4]) 
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.phaseLock_ip(:,3,3:40,:));
cfg.zli = zl ;
cfg.cluster_mask = false;
cfg.npermutations =10
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.5);
set(gca,'TickDir','out');

xlim([-0.5 0.5])

if mFig;  yticklabels({}); xlabel({}); xticklabels({}) 
else 
    title('ipi')
    colorbar('eastoutside')	
    ylabel('ALpha lateralisation')    
    xlabel('Time from onset of saccade')
end
ax = gca;
ax.XTick = [-0.5 : 0.25: 0.5];
ax.YTick = [3 10 20 30 40];
ax.TickLength = [0.015 0.015];

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.5);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on


subplot(2,6,[5:6]) 
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.phaseLock_ip(:,3,3:40,:) - GA_struct.phaseLock_co(:,3,3:40,:));
cfg.zli = [-0.15 0.15];
%cfg.cluster_mask = false;
cfg.npermutations =1000
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.5);
set(gca,'TickDir','out');

xlim([-0.5 0.5])

if mFig;  yticklabels({}); xlabel({}); xticklabels({}) 
else 
    title('ipsi vs. contra')
    colorbar('eastoutside')	
    ylabel('ALpha lateralisation')    
    xlabel('Time from onset of saccade')
end
ax = gca;
ax.XTick = [-0.5 : 0.25: 0.5];
ax.YTick = [3 10 20 30 40];
ax.TickLength = [0.015 0.015];

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.5);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on

%
load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_itpc_topo.mat')

time_i = [0 0.25];
freq_i = [8:12];

for freqInd = 2:2
    if freqInd == 1
        freq_i = [3:7];
        zlim = [-0 0.4];
    elseif freqInd == 2
        freq_i = [8:12];
        zlim = [-0 0.25];
    end
for directInd = 1:3
    subplot(2,6,8+ directInd)
    %data2plot1 = squeeze(mean(GA_struct.diff_pow(:,1,2,1,:,:,:)));
    %data2plot2 = squeeze(mean(GA_struct.diff_pow(:,1,3,1,:,:,:)));

    dummy = [];
    dummy.time = GA_struct.time;
    dummy.freq = GA_struct.freq;
    dummy.label = GA_struct.eleLabel;
    dummy.dimord = 'chan_freq_time';

        dummy.data2plot1 = squeeze(mean(GA_struct.phaseLock_left(:,3,:,:,:))); % mean over trials (reduce to 3 dimensions: electrodes, frequency points, time points (on which we centred our sliding time windows)

        dummy.data2plot2 = squeeze(mean(GA_struct.phaseLock_right(:,3,:,:,:)));
        
        dummy.data2plot3 = squeeze(mean(GA_struct.phaseLock_left(:,3,:,:,:))) -squeeze(mean(GA_struct.phaseLock_right(:,3,:,:,:)));
    % now lets plot it

    cmap = brewermap([],'*RdBu');
    colormap(cmap)
    cfg = [];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.marker = 'off';;
     cfg.xlim = time_i;
    cfg.figure = 'gcf';
    cfg.colormap = cmap;
    
    if directInd == 1
        cfg.parameter = 'data2plot1';
    elseif directInd == 2
        cfg.parameter = 'data2plot2';
    elseif directInd == 3
        cfg.parameter = 'data2plot3';
        zlim = [-0.12 0.12] ;
    end
    
     
    cfg.ylim = freq_i;
    cfg.zlim =  zlim;
    ft_topoplotTFR(cfg, dummy); %title([num2str(timeStart(timeInd)) ' to ' num2str(timeStart(timeInd)+ timewin)])
    if ~mFig
        colorbar
    end
end
end
%% 
set(gcf, 'renderer', 'Painters'); 
saveas(gcf, ['/Users/baiweil/Library/Mobile Documents/com~apple~CloudDocs/Projects/Current Projects/ReAnalysis/EEG_eye_delay_NN/figs_ms/fig_sup_itpc_coVSip_return.eps'], 'epsc');

