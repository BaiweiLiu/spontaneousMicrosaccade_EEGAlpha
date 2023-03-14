%% figure 3 co vs. ip
load('/Users/baiweil/AnalysisDesk/Spontaneous_alpha/data2plot/GA_saccSize_cvsi_TF-PO78.mat')
%%
GA_struct_cvsi= GA_struct;
mFig = true;
freq_i = [8:12];
% normlize the data 
power_co_norm = nan(size(GA_struct.power_co));
power_ip_norm = nan(size(GA_struct.power_co));

t_base= [-0.5 0.5];
t_base_pre = [-0.5 -0.15];

for subjInd = 1:23 
    for condInd = 1:3
        for freqInd = 1:50
            
            %% raw data
            t_sel = dsearchn(GA_struct.time',t_base')';
            trl_data =  squeeze(GA_struct.power_co(subjInd,condInd, freqInd,:)); 
            trl_base =  nanmean(trl_data([t_sel(1):t_sel(2)]));
            trl_norm = (trl_data - trl_base) ./(trl_data + trl_base) * 100;
            power_co_norm_epoch(subjInd,condInd, freqInd,:) = trl_norm;
            
            t_sel = dsearchn(GA_struct.time',t_base_pre')';
            trl_data =  squeeze(GA_struct.power_co(subjInd,condInd, freqInd,:)); 
            trl_base =  nanmean(trl_data([t_sel(1):t_sel(2)]));
            trl_norm = (trl_data - trl_base) ./(trl_data + trl_base) * 100;
            power_co_norm_pre(subjInd,condInd, freqInd,:) = trl_norm;
            
            t_sel = dsearchn(GA_struct.time',t_base')';
            trl_data =  squeeze(GA_struct.power_ip(subjInd,condInd, freqInd,:)); 
            trl_base =  nanmean(trl_data([t_sel(1):t_sel(2)]));
            trl_norm = (trl_data - trl_base) ./(trl_data + trl_base) * 100;
            power_ip_norm_epoch(subjInd,condInd, freqInd,:) = trl_norm;
            
            t_sel = dsearchn(GA_struct.time',t_base_pre')';
            trl_data =  squeeze(GA_struct.power_ip(subjInd,condInd, freqInd,:)); 
            trl_base =  nanmean(trl_data([t_sel(1):t_sel(2)]));
            trl_norm = (trl_data - trl_base) ./(trl_data + trl_base) * 100;
            power_ip_norm_pre(subjInd,condInd, freqInd,:) = trl_norm;
            
            %% log data
            
            t_sel = dsearchn(GA_struct.time',t_base')';
            trl_data =  squeeze(GA_struct.power_co_log(subjInd,condInd, freqInd,:)); 
            trl_base =  nanmean(trl_data([t_sel(1):t_sel(2)]));
            trl_norm = trl_data - trl_base;
            power_co_norm_epoch_log(subjInd,condInd, freqInd,:) = trl_norm;
            
            t_sel = dsearchn(GA_struct.time',t_base_pre')';
            trl_data =  squeeze(GA_struct.power_co_log(subjInd,condInd, freqInd,:)); 
            trl_base =  nanmean(trl_data([t_sel(1):t_sel(2)]));
            trl_norm = trl_data - trl_base;
            power_co_norm_pre_log(subjInd,condInd, freqInd,:) = trl_norm;
            
            t_sel = dsearchn(GA_struct.time',t_base')';
            trl_data =  squeeze(GA_struct.power_ip_log(subjInd,condInd, freqInd,:)); 
            trl_base =  nanmean(trl_data([t_sel(1):t_sel(2)]));
            trl_norm = trl_data - trl_base;
            power_ip_norm_epoch_log(subjInd,condInd, freqInd,:) = trl_norm;
            
            t_sel = dsearchn(GA_struct.time',t_base_pre')';
            trl_data =  squeeze(GA_struct.power_ip_log(subjInd,condInd, freqInd,:)); 
            trl_base =  nanmean(trl_data([t_sel(1):t_sel(2)]));
            trl_norm = trl_data - trl_base;
            power_ip_norm_pre_log(subjInd,condInd, freqInd,:) = trl_norm;
        end
    end 
end
    
GA_struct.power_co_norm_epoch = power_co_norm_epoch;
GA_struct.power_ip_norm_epoch = power_ip_norm_epoch;

GA_struct.power_co_norm_pre = power_co_norm_pre;
GA_struct.power_ip_norm_pre = power_ip_norm_pre;

GA_struct.power_co_norm_epoch_log = power_co_norm_epoch_log;
GA_struct.power_ip_norm_epoch_log = power_ip_norm_epoch_log;

GA_struct.power_co_norm_pre_log = power_co_norm_pre_log;
GA_struct.power_ip_norm_pre_log = power_ip_norm_pre_log;

%%

timeLim = [-0.5 0.5] ; 
figure('position', [100 100 1200 800], 'color', [1 1 1])

subplot(3,2,1)
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.power_co(:,3,3:40,:));
cfg.zli = [0 0.000001]* 1.05;
cfg.cluster_mask = false;
cfg.npermutations = 10
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
    set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
    set(gca,'TickDir','out');

if mFig; xticklabels({}); ylabel({}); yticklabels({}) 
else 
    colorbar('eastoutside')	
    ylabel('Freq')    
    xlabel('time from onset of saccade')
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
%title('con')
xlim([-0.5 0.5])
    
    
subplot(3,2,2) 
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.power_ip(:,3,3:40,:));
cfg.zli = [0 0.000001] * 1.05;
%cfg.cluster_mask = false;
cfg.npermutations =10;
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
    set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
    set(gca,'TickDir','out');

if mFig; xticklabels({}); yticklabels({}); ylabel({});
else 
    colorbar('eastoutside')	
    ylabel('Freq')    
    xlabel('time from onset of saccade')
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

%title('ipi')
xlim([-0.5 0.5])

%colormap('jet')

%
timeLim = [-0.5 0.5]; 
zli = [-0.3 0.3];

subplot(3,2,3)
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.power_co_norm_pre_log(:,3,3:40,:));
cfg.zli = zli;
cfg.cluster_mask = true;
cfg.npermutations =10000;
stat = eegFuture_tfPlot(cfg); hold on
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
    set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
    set(gca,'TickDir','out');

if mFig; xticklabels({}); ylabel({}); yticklabels({}) 
else 
    colorbar('eastoutside')	
    ylabel('Freq')    
    xlabel('time from onset of saccade')
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
xlim([-0.5 0.5])
    
    
subplot(3,2,4) 
cfg = [];
cfg.time = GA_struct.time;
cfg.time_i = timeLim;
cfg.freq = GA_struct.freq(3:40);
cfg.tf_data = squeeze(GA_struct.power_ip_norm_pre_log(:,3,3:40,:));
cfg.zli = zli;
cfg.cluster_mask = true;
cfg.npermutations =10000;
stat = eegFuture_tfPlot(cfg); hold on
stat_ip_tf = stat;
plot([0 0], [GA_struct.freq(3) GA_struct.freq(40)], '--k')
%plot(timeLim, [8 8], ':k');
%plot(timeLim, [12 12], ':k');
set(gca,'LineWidth',3)
    set(gca, 'FontSize', 20,'Fontname','Arial', 'LineWidth', 1.2);
    set(gca,'TickDir','out');

if mFig; xticklabels({}); yticklabels({}); ylabel({});
else 
    colorbar('eastoutside')	
    ylabel('Freq')    
    xlabel('time from onset of saccade')
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

xlim([-0.5 0.5])

%
yl = [-0.3 0.6];
data2plot_co= squeeze(nanmean(GA_struct.power_co_norm_pre_log(:,3,freq_i,:),3));
data2plot_ip= squeeze(nanmean(GA_struct.power_ip_norm_pre_log(:,3,freq_i,:),3));
subplot(3,2,5) 
m1= lineErr_plot(GA_struct.time,data2plot_co, 'k', 'ci'); hold on
xlim([-0.5 0.5])
ylim(yl)
plot([0 0], ylim, '--k')
plot(xlim, [0 0], '--k')

time2plot = GA_struct.time >= -0.5 & GA_struct.time <= 0.5;
compData = zeros(size(data2plot_co(:,time2plot)));
cfg = [];
cfg.xax = GA_struct.time(time2plot);
cfg.npermutations = 10;
cfg.clusterStatEvalaluationAlpha= 0.025;
cfg.nsub=23;
cfg.statMethod = 'montecarlo';
cfg.time_i = [-0.5,0.5];
state_t = cluster_perm_1D(cfg,data2plot_co(:,time2plot),compData);

mask_xxx = double(state_t.mask); mask_xxx(mask_xxx==0) = nan;
plot(GA_struct.time(time2plot), mask_xxx * -0.27, 'k', 'LineWidth', 3);
ax = gca; ax.XTick = [-0.5 : 0.25: 0.5]; ax.YTick = [-0.3 : 0.3: 0.6]; 

if mFig;  xlabel({}); yticklabels({});xticklabels({}) 
else 
    colorbar('eastoutside')	
    ylabel('ALpha lateralisation')    
    xlabel('Time from onset of saccade')
end
subplot(3,2,6) 
m2 = lineErr_plot(GA_struct.time,data2plot_ip, 'k', 'ci');
xlim([-0.5 0.5])
ylim(yl)
plot([0 0], ylim, '--k')
plot(xlim, [0 0], '--k')
cfg = [];
cfg.xax = GA_struct.time(time2plot);
cfg.npermutations = 10;
cfg.clusterStatEvalaluationAlpha= 0.025;
cfg.nsub=23;
cfg.statMethod = 'montecarlo';
cfg.time_i = [-0.5,0.5];
state_t = cluster_perm_1D(cfg,data2plot_ip(:,time2plot),compData);
stat_ip_time = state_t;
mask_xxx = double(state_t.mask); mask_xxx(mask_xxx==0) = nan;
plot(GA_struct.time(time2plot), mask_xxx * -0.27, 'k', 'LineWidth', 3);
ax = gca; ax.XTick = [-0.5 : 0.25: 0.5]; ax.YTick = [-0.3 : 0.3: 0.6]; 
if mFig;  yticklabels({}); xlabel({}); xticklabels({}) 
else 
    colorbar('eastoutside')	
    ylabel('ALpha lateralisation')    
    xlabel('Time from onset of saccade')
end
%% 
set(gcf, 'renderer', 'Painters'); 
saveas(gcf, ['/Users/baiweil/Library/Mobile Documents/com~apple~CloudDocs/Projects/Current Projects/ReAnalysis/EEG_eye_delay_NN/figs_ms/fig_sup_coVSip_return.eps'], 'epsc');