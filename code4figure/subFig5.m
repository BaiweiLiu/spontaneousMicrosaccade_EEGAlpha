%% figure 4 plot 
load('/Users/baiweil/AnalysisDesk/Spontaneous_alphas/data2plot/GA_saccSize_cvsi_TF-O12.mat')
GA_struct_O12= GA_struct;

%% itpc PO78
GA_struct = GA_struct_O12;
mFig = true
timeLim = [-0.5 0.5] ; 
zl = [0 0.4];
figure('position', [100 100 1500 250], 'color', [1 1 1])

subplot(1,3,1)
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.phaseLock_co(:,2,3:40,:));
cfg.zli = zl;
cfg.cluster_mask = false;
cfg.npermutations =10
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
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
ax.TickLength = [0.015 0.015];

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.2);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
    
    
subplot(1,3,2) 
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.phaseLock_ip(:,2,3:40,:));
cfg.zli = zl ;
cfg.cluster_mask = false;
cfg.npermutations =10
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
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
ax.TickLength = [0.015 0.015];

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.2);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on


subplot(1,3,3) 
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.phaseLock_ip(:,2,3:40,:) - GA_struct.phaseLock_co(:,2,3:40,:));
cfg.zli = [-0.15 0.15];
%cfg.cluster_mask = false;
cfg.npermutations =1000
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
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
ax.TickLength = [0.015 0.015];

box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','LineWidth', 1.2);

set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on

%% 
set(gcf, 'renderer', 'Painters'); 
saveas(gcf, ['/Users/baiweil/Library/Mobile Documents/com~apple~CloudDocs/Projects/Current Projects/ReAnalysis/EEG_eye_delay_NN/figs_ms/fig5_sup_itpc_O12.eps'], 'epsc');