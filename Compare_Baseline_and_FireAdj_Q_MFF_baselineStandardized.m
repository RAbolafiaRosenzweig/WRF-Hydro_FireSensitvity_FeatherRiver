clc;clear all;close all;

cms_to_cfs_conversion_factor = 35.314666212661;
met_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/Feather_Forcing/';
WYs = 2000:2022;
nyears = length(WYs);


%% plot info during these years for: MIDDLE FORK
STA_IDs = {'DCS','F56','F57','GYB','ICR','JBR','LCB','MER','MFP','NFP','NYS','ORH','SPK','TM1','TM2','TM3','WFR','YPB','YRS'};
current_station = 'MER'; %station closest to outlet of middle fork feather
IDX_sta = strcmp(STA_IDs,current_station);
idx = find(IDX_sta==1);
store_sta_idx = idx;

%% get baseline sim data

%load in simulated Q - baseline:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Feather_Baseline/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
Baseline_Streamflow = Simulated_Q.Streamflow;
Baseline_Streamflow = Baseline_Streamflow(:,store_sta_idx);
sim_dates = Simulated_Q.dates;
sim_datevecs = datevec(sim_dates);

%% get realistic sim data
%load in simulated Q - realistic adjustments - CLEANER:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Realistic/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
Realistic_Streamflow = Simulated_Q.Streamflow;
Realistic_Streamflow = Realistic_Streamflow(:,store_sta_idx);

%% average sim data to daily:
[u,~,j] = unique(sim_datevecs(:,1:3),'rows','stable');
baseline_Q_daily = accumarray(j,Baseline_Streamflow,[],@nanmean);
realistic_Q_daily = accumarray(j,Realistic_Streamflow,[],@nanmean);
sim_dates_daily = datenum(u);

%% get obs data:
obs_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/Observations/CADWR/Streamflow_csv_files/';
store_obs_Q=[];
for WY=2000:2022
    infilename = sprintf('%s_WY%d.csv',current_station,WY);
    if exist([obs_dir,infilename],'file') > 0
        Data = readtable([obs_dir,infilename]);
        dates = Data.POSIXct;
        q = Data.obs;
        datenums = datenum(dates);
        datevecs = datevec(datenums);
        
        %aggregate to daily
        [u,~,j] = unique(datevecs(:,1:3),'rows','stable');
        Q_daily = accumarray(j,q,[],@nanmean);
        store_obs_Q = [store_obs_Q;datenum(u),Q_daily];
    end
end

%% temporally match obs and sim data:
[c,ia,ib] = intersect(store_obs_Q(:,1),sim_dates_daily,'rows');
matched_dates = c;
matched_datevecs = datevec(matched_dates);
Baseline_Q_daily_matched = baseline_Q_daily(ib);
Realistic_Q_daily_matched = realistic_Q_daily(ib);
obs_Q_daily_matched = store_obs_Q(ia,2)

%% convert to anomalies:
post_fire_start = datenum([2021 11 1]);
pre_fire_end = datenum([2020 8 1]);
idx_postfire = find(matched_dates >= post_fire_start);
idx_prefire = find(matched_dates < pre_fire_end);
idx_mid=find(matched_dates<post_fire_start & matched_dates>=pre_fire_end);
Baseline_Q_daily_matched(idx_mid) = NaN;
Realistic_Q_daily_matched(idx_mid) = NaN;
obs_Q_daily_matched(idx_mid) = NaN;

postfire_Realistic_Q_daily_matched = Realistic_Q_daily_matched(idx_postfire);
Realistic_Q_daily_matched = Baseline_Q_daily_matched;
Realistic_Q_daily_matched(idx_postfire) = postfire_Realistic_Q_daily_matched;

Baseline_Q_daily_anoms = (Baseline_Q_daily_matched - nanmean(Baseline_Q_daily_matched(idx_prefire)))./nanstd(Baseline_Q_daily_matched(idx_prefire));
Realistic_Q_daily_anoms = (Realistic_Q_daily_matched - nanmean(Realistic_Q_daily_matched(idx_prefire)))./nanstd(Realistic_Q_daily_matched(idx_prefire));
Obs_Q_daily_anoms = (obs_Q_daily_matched - nanmean(obs_Q_daily_matched(idx_prefire)))./nanstd(obs_Q_daily_matched(idx_prefire));

%% plot scatter plot comparrissons (baseline):
f=figure;
hold on
p1=plot(Baseline_Q_daily_anoms,Obs_Q_daily_anoms,'o','markerfacecolor',[0.4660, 0.6740, 0.1880],'markeredgecolor',[0.4660, 0.6740, 0.1880],'markersize',2);
p2=plot(Baseline_Q_daily_anoms(idx_postfire),Obs_Q_daily_anoms(idx_postfire),'bo','markerfacecolor','b','markersize',7);
p3=plot(Realistic_Q_daily_anoms(idx_postfire),Obs_Q_daily_anoms(idx_postfire),'ro','markerfacecolor','r','markersize',3);
plot(linspace(-100,100,10),linspace(-100,100,10),'--k')

grid on
box on
f.Position = [80   171   973   614];
xlabel('Simulated daily Q anomalies')
%ylabel('Observed daily Q anomalies')
set(gca,'fontsize',26)
xlim([-0.6 10])
ylim([-0.6 10])
[r_pre,p_pre] = corr(Baseline_Q_daily_anoms(idx_prefire) , Obs_Q_daily_anoms(idx_prefire),'type','pearson','rows','complete');
mean_bias_pre = nanmean(Baseline_Q_daily_anoms(idx_prefire) - Obs_Q_daily_anoms(idx_prefire));
prctile10_bias_pre = prctile(Baseline_Q_daily_anoms(idx_prefire) - Obs_Q_daily_anoms(idx_prefire),10);
prctile90_pre = prctile(Baseline_Q_daily_anoms(idx_prefire) - Obs_Q_daily_anoms(idx_prefire),90);

[r_post,p_post] = corr(Baseline_Q_daily_anoms(idx_postfire) , Obs_Q_daily_anoms(idx_postfire),'type','pearson','rows','complete');
mean_bias_post = nanmean(Baseline_Q_daily_anoms(idx_postfire) - Obs_Q_daily_anoms(idx_postfire));
prctile10_bias_post = prctile(Baseline_Q_daily_anoms(idx_postfire) - Obs_Q_daily_anoms(idx_postfire),10);
prctile90_post = prctile(Baseline_Q_daily_anoms(idx_postfire) - Obs_Q_daily_anoms(idx_postfire),90);

%baseline prefire stats
if p_pre>0.001
    text(0,9.5,sprintf('prefire r = %.2f (p=%.3f)',round(r_pre,2),round(p_pre,3)),'color',[0 0.5 0],'fontsize',26)
else
    text(0,9.5,sprintf('prefire r = %.2f (p<0.001)',round(r_pre,2)),'color',[0 0.5 0],'fontsize',26)
end
% % text(0,8.5,sprintf('prefire bias = %.3f [%.3f-%.3f]',round(mean_bias_pre,3),round(prctile10_bias_pre,3),round(prctile90_pre,3)),'color',[0 0.5 0],'fontsize',20)

%baseline postfire stats
if p_post>0.001
    text(5.5,2,sprintf('postfire r = %.2f (p=%.3f)',round(r_post,2),round(p_post,3)),'color','b','fontsize',26)
else
    text(5.5,2,sprintf('postfire r = %.2f (p<0.001)',round(r_post,2)),'color','b','fontsize',26)
end
% % text(5,2.8,sprintf('postfire bias = %.3f [%.3f-%.3f]',round(mean_bias_post,3),round(prctile10_bias_post,3),round(prctile90_post,3)),'color','b','fontsize',20)

%realistic postfire stats
[r_post,p_post] = corr(Realistic_Q_daily_anoms(idx_postfire) , Obs_Q_daily_anoms(idx_postfire),'type','pearson','rows','complete');
mean_bias_post = nanmean(Realistic_Q_daily_anoms(idx_postfire) - Obs_Q_daily_anoms(idx_postfire));
prctile10_bias_post = prctile(Realistic_Q_daily_anoms(idx_postfire) - Obs_Q_daily_anoms(idx_postfire),10);
prctile90_post = prctile(Realistic_Q_daily_anoms(idx_postfire) - Obs_Q_daily_anoms(idx_postfire),90);
if p_post>0.001
    text(5.5,1.2,sprintf('postfire r = %.2f (p=%.3f)',round(r_post,2),round(p_post,3)),'color','r','fontsize',26)
else
    text(5.5,1.2,sprintf('postfire r = %.2f (p<0.001)',round(r_post,2)),'color','r','fontsize',26)
end
% % text(5,0.8,sprintf('postfire bias = %.3f [%.3f-%.3f]',round(mean_bias_post,3),round(prctile10_bias_post,3),round(prctile90_post,3)),'color','r','fontsize',20)

legend([p1 p2 p3],{'prefire baseline','post-fire baseline','post-fire realistic'},'fontsize',40,'location','northeast')

saveas(f,'/Users/abolafia/ASO_Fire/Plots/Baseline_and_Realistic_MFF_Qanom_scatterprefireSTD.png')

%% create histogram plot of biases:
bias_pre_baseline = (Baseline_Q_daily_anoms(idx_prefire) - Obs_Q_daily_anoms(idx_prefire));
bias_post_baseline = (Baseline_Q_daily_anoms(idx_postfire) - Obs_Q_daily_anoms(idx_postfire));
bias_post_realistic = (Realistic_Q_daily_anoms(idx_postfire) - Obs_Q_daily_anoms(idx_postfire));

plot_data = nan(length(bias_pre_baseline),3);
plot_data(:,1) = bias_pre_baseline;
plot_data(1:length(bias_post_baseline),2) = bias_post_baseline;
plot_data(1:length(bias_post_baseline),3) = bias_post_realistic;

f=figure;
hold on
b=boxplot(plot_data,'Notch','off','whisker',1,'labels',{'prefire baseline','postfire baseline','postfire realistic'});
ylim([-0.25 0.25])
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   

outlier1 = a(1);
outlier2 = a(2);
outlier3 = a(3);
set(outlier1, 'Marker','o','markerfacecolor','r','markeredgecolor','r','markersize',1.5); 
set(outlier2, 'Marker','o','markerfacecolor','b','markeredgecolor','b','markersize',1.5); 
set(outlier3, 'Marker','o','markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0],'markersize',1.5); 

median1 = a(4);
median2 = a(5);
median3 = a(6);
set(median1, 'Color', 'r','linewidth',3); 
set(median2, 'Color', 'b','linewidth',3); 
set(median3, 'Color', [0 0.5 0],'linewidth',3); 

box1 = a(7);
box2 = a(8);
box3 = a(9);
set(box1, 'Color', 'r','linewidth',3); 
set(box2, 'Color', 'b','linewidth',3); 
set(box3, 'Color', [0 0.5 0],'linewidth',3); 

LL1 = a(10);
LL2 = a(11);
LL3 = a(12);
set(LL1, 'Color', 'r','linewidth',3); 
set(LL2, 'Color', 'b','linewidth',3); 
set(LL3, 'Color', [0 0.5 0],'linewidth',3); 

UL1 = a(13);
UL2 = a(14);
UL3 = a(15);
set(UL1, 'Color', 'r','linewidth',3); 
set(UL2, 'Color', 'b','linewidth',3); 
set(UL3, 'Color', [0 0.5 0],'linewidth',3); 

lw1 = a(16);
lw2 = a(17);
lw3 = a(18);
set(lw1, 'Color', 'r','linewidth',3); 
set(lw2, 'Color', 'b','linewidth',3); 
set(lw3, 'Color', [0 0.5 0],'linewidth',3); 

uw1 = a(19);
uw2 = a(20);
uw3 = a(21);
set(uw1, 'Color', 'r','linewidth',3); 
set(uw2, 'Color', 'b','linewidth',3); 
set(uw3, 'Color', [0 0.5 0],'linewidth',3); 

f.Position = [139         202        1074         595];
set(gca,'fontsize',22)
grid on 
box on
xlim([0.5 3.5])
plot(linspace(0,5,10),linspace(0,0,10),'--k','linewidth',1.5)

%show mean biases
plot(1,nanmean(bias_pre_baseline),'*','markersize',22,'color',[0 0.5 0]);
plot(2,nanmean(bias_post_baseline),'*','markersize',22,'color','b');
plot(3,nanmean(bias_post_realistic),'*','markersize',22,'color','r');
plot(1,nanmean(bias_pre_baseline),'.','markersize',22,'color',[0 0.5 0]);
plot(2,nanmean(bias_post_baseline),'.','markersize',22,'color','b');
plot(3,nanmean(bias_post_realistic),'.','markersize',22,'color','r');


annotation('textbox', [0.15 0.8 0.1 0.1], 'String', sprintf('prefire baseline mean bias = %.3f',round(nanmean(bias_pre_baseline),3)), 'EdgeColor', 'none','fontsize',18,'backgroundcolor','w','color',[0 0.5 0]);
annotation('textbox', [0.4 0.75 0.1 0.1], 'String', sprintf('postfire baseline mean bias = %.3f',round(nanmean(bias_post_baseline),3)), 'EdgeColor', 'none','fontsize',18,'backgroundcolor','w','color','b');
annotation('textbox', [0.6 0.7 0.1 0.1], 'String', sprintf('postfire realistic mean bias = %.3f',round(nanmean(bias_post_realistic),3)), 'EdgeColor', 'none','fontsize',18,'backgroundcolor','w','color','r');

%ylabel('daily anomaly Q biases')

saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_and_Realistic_MFF_Bias_BoxplotprefireSTD.png')

%% show time series for pre fire comparissons:
dates_prefire= matched_dates(idx_prefire);
prefire_obs = Obs_Q_daily_anoms(idx_prefire);
prefire_baseline = Baseline_Q_daily_anoms(idx_prefire);
% % prefire_obs = (obs_Q_daily_matched(idx_prefire) - mean(obs_Q_daily_matched(idx_prefire)));
% % prefire_baseline = (Baseline_Q_daily_matched(idx_prefire) - mean(Baseline_Q_daily_matched(idx_prefire)));

f=figure;
hold on
p1=plot(dates_prefire,prefire_obs,'-k','linewidth',3);
p2=plot(dates_prefire,prefire_baseline,'-b','linewidth',1.5);
grid on
box on
xlim([min(dates_prefire) max(dates_prefire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',40)
%ylabel('pre-fire Q anomalies')
% % legend([p1 p2],{'obs','baseline'},'location','best','fontsize',22)
f.Position = [-1803         325        1804         435];
%report R:
[R,p] = corr(prefire_obs,prefire_baseline,'type','pearson','rows','complete');
%report NSE:
NSE =nashsutcliffe([dates_prefire,prefire_baseline], [dates_prefire,prefire_obs]);
%report bias:
bias =nanmean(prefire_baseline - prefire_obs);
% if p>0.001
%     text(dates_prefire(1000),20,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias - %.2f',round(R^2,2),round(p,2),round(NSE,2),round(bias,2)),'fontsize',40,'color','b')
% else
%     text(dates_prefire(1000),20,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f; bias = %.2f',round(R^2,2),round(NSE,2),round(bias,2)),'fontsize',40,'color','b')
% end
xtickangle(35)
saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_MFF_Qanom_timeseries_prefireprefireSTD.png')

%% show time series for post fire comparissons:
dates_postfire= matched_dates(idx_postfire);
postfire_obs = Obs_Q_daily_anoms(idx_postfire);
postfire_baseline = Baseline_Q_daily_anoms(idx_postfire);
postfire_realistic = Realistic_Q_daily_anoms(idx_postfire);

% % postfire_obs = (obs_Q_daily_matched(idx_postfire) - mean(obs_Q_daily_matched(idx_postfire)));
% % postfire_baseline = (Baseline_Q_daily_matched(idx_postfire) - mean(Baseline_Q_daily_matched(idx_postfire)));
% % postfire_realistic = (Realistic_Q_daily_matched(idx_postfire) - mean(Realistic_Q_daily_matched(idx_postfire)));

f=figure;
hold on
p1=plot(dates_postfire,postfire_obs,'-k','linewidth',3);
p2=plot(dates_postfire,postfire_baseline,'-b','linewidth',1.5);
p3=plot(dates_postfire,postfire_realistic,'-r','linewidth',1.5);

grid on
box on
xlim([min(dates_postfire) max(dates_postfire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',40)
%ylabel('post-fire Q anomalies')
% % legend([p1 p2 p3],{'obs','baseline','Mod-params+GVF+Veg-class+Snow-alb'},'location','northeast','fontsize',22)
f.Position = [-1803         325         900         435];

%report R:
[R_baseline,p_baseline] = corr(postfire_obs,postfire_baseline,'type','pearson','rows','complete');
%report NSE:
NSE_baseline =nashsutcliffe([dates_postfire,postfire_baseline], [dates_postfire,postfire_obs]);
%report bias:
bias_baseline = nanmean(postfire_baseline - postfire_obs);
% if p>0.001
%     text(dates_postfire(150),1.5,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias = %.2f',round(R_baseline^2,2),round(p_baseline,2),round(NSE_baseline,2),round(bias_baseline,2)),'fontsize',20,'color','b')
% else
%     text(dates_postfire(150),1.5,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f;bias = %.2f',round(R_baseline^2,2),round(NSE_baseline,2),round(bias_baseline,2)),'fontsize',20,'color','b')
% end

%report R:
[R_realistic,p_realistic] = corr(postfire_obs,postfire_realistic,'type','pearson','rows','complete');
%report NSE:
NSE_realistic =nashsutcliffe([dates_postfire,postfire_realistic], [dates_postfire,postfire_obs]);
%report bias:
bias_realistic = nanmean(postfire_realistic - postfire_obs);
% if p>0.001
%     text(dates_postfire(150),1,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias = %.2f',round(R_realistic^2,2),round(p_realistic,2),round(NSE_realistic,2),round(bias_realistic,2)),'fontsize',20,'color','r')
% else
%     text(dates_postfire(150),1,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f; bias = %.2f',round(R_realistic^2,2),round(NSE_realistic,2),round(bias_realistic,2)),'fontsize',20,'color','r')
% end
xtickangle(35)
saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/MFF_Qanom_timeseries_postfireprefireSTD.png')

mean_postfire_anom = nanmean(postfire_obs);
max_postfire_anom = nanmax(postfire_obs);

%% use wilcoxon rank sum test to see difference in the biases:
realistic_bias = postfire_realistic - postfire_obs;
baseline_bias = postfire_baseline - postfire_obs;
% null: it is equally likely that a randomly selected value of one population will be lesser or greater than a randomly selected value from a second population
p_bias_2sided = ranksum(realistic_bias,baseline_bias);
sprintf('difference in bias is significant with p = %.3f for daily validation',p_bias_2sided)

%% use permutation test for 2 independent correlations:
delta_r = R_realistic - R_baseline;
n_perms = 10000;

L_baseline = length(postfire_baseline);
store_delta_r=[];
for j=1:n_perms
    permuted_baseline = postfire_baseline(randperm(L_baseline));
    permuted_enhanced = postfire_baseline(randperm(L_baseline));

    baseline_r_tmp = corr(permuted_baseline,postfire_obs,'type','pearson');
    realistic_r_tmp = corr(permuted_enhanced,postfire_obs,'type','pearson');
    
    store_delta_r = [store_delta_r; realistic_r_tmp - baseline_r_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_r) > abs(delta_r));
sprintf('difference in r is significant with p = %.3f for daily validation',length(idx)/n_perms)

%% use permutation test for 2 independent NSE's:
delta_nse = NSE_realistic - NSE_baseline;
n_perms = 10000;

store_delta_nse=[];
for j=1:n_perms
    permuted_baseline = postfire_baseline(randperm(L_baseline));
    permuted_enhanced = postfire_baseline(randperm(L_baseline));

    baseline_nse_tmp = nashsutcliffe([dates_postfire,permuted_baseline],[dates_postfire,postfire_obs]);
    realistic_nse_tmp = nashsutcliffe([dates_postfire,permuted_enhanced],[dates_postfire,postfire_obs]);
    
    store_delta_nse = [store_delta_nse; realistic_nse_tmp - baseline_nse_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_nse) > abs(delta_nse));
sprintf('difference in nse is significant with p = %.3f for daily validation',length(idx)/n_perms)

%% report baseline stats from years with a similar dryness as post-fire periods:
post_fire_start_day = datevec(dates_postfire(1));
npostfire = length(dates_postfire);

pre_fire_datevecs = datevec(dates_prefire);
idx1=find(pre_fire_datevecs(:,2) == post_fire_start_day(:,2) & pre_fire_datevecs(:,3) == post_fire_start_day(:,3));
store_stats=[];
for j=1:length(idx1)-1
    current_start_idx = idx1(j);
    current_end_idx = current_start_idx+(npostfire-1);
    current_baseline_data = prefire_baseline(current_start_idx:current_end_idx);
    current_obs_data = prefire_obs(current_start_idx:current_end_idx);
    current_dates =dates_prefire(current_start_idx:current_end_idx);
    
    Mean_obs = nanmean(current_obs_data);
    Peak_obs = max(current_obs_data);
    current_R2 = corr(current_baseline_data,current_obs_data,'rows','complete')^2;
    current_NSE = nashsutcliffe([current_dates,current_baseline_data], [current_dates,current_obs_data]);
    current_bias = nanmean(current_baseline_data -current_obs_data);
    store_stats = [store_stats;Mean_obs,current_R2,current_NSE,current_bias,Peak_obs];
end
[R_meananom_R2,p] = corr(store_stats(:,1) , store_stats(:,2));
sprintf('R(mean anom , R2) = %.4f,p=%.4f',R_meananom_R2,p)

R_meananom_NSE = corr(store_stats(:,1) , store_stats(:,3));
sprintf('R(mean anom , NSE) = %.4f',R_meananom_NSE)

search_rng = [mean_postfire_anom-0.2*mean_postfire_anom , mean_postfire_anom+0.2*mean_postfire_anom];
% idx=find(store_stats(:,1) >= search_rng(2) & store_stats(:,1) <= search_rng(1));
idx=find(store_stats(:,1) < 0 );
similar_stats = store_stats(idx,:);
mean_R2_similar = nanmean(similar_stats(:,2));
mean_NSE_similar = nanmean(similar_stats(:,3));
mean_bias_similar = nanmean(similar_stats(:,4));
sprintf('mean [R2, NSE, bias] for similar years = [%.4f,%.4f,%.4f]',mean_R2_similar,mean_NSE_similar,mean_bias_similar)

%create corresponding plot:
idx_low_flow = find(store_stats(:,1) < 0 & store_stats(:,5) > 1 & store_stats(:,5) < 3);
low_flow_stats = store_stats(idx_low_flow,:);

store=[];
for k=1:length(idx_low_flow)
    j = idx_low_flow(k);
    current_start_idx = idx1(j);
    current_end_idx = current_start_idx+(npostfire-1);
    current_baseline_data = prefire_baseline(current_start_idx:current_end_idx);
    current_obs_data = prefire_obs(current_start_idx:current_end_idx);
    current_dates =dates_prefire(current_start_idx:current_end_idx);
    disp('end date')
    datevec(current_dates(end))

    Mean_obs = nanmean(current_obs_data);
    max_obs = max(current_obs_data);
    assert(Mean_obs <0,'wrong period selected')
    store=[store;Mean_obs,max_obs];

    f=figure;
    hold on
    p1=plot(current_dates,current_obs_data,'-k','linewidth',3);
    p2=plot(current_dates,current_baseline_data,'-b','linewidth',1.5);
    grid on
    box on
    xlim([min(current_dates) max(current_dates)])
    ylim([-1 3.5])
    datetick('x','mm/yyyy','keepticks','keeplimits')
    set(gca,'fontsize',40)
    %ylabel('pre-fire Q anomalies')
    % % legend([p1 p2],{'obs','baseline'},'location','best','fontsize',22)
    f.Position = [-1803         325        900         435];
    %report R:
    [R,p] = corr(current_obs_data,current_baseline_data,'type','pearson');
    %report NSE:
    NSE =nashsutcliffe([current_dates,current_baseline_data], [current_dates,current_obs_data]);
    %report bias:
    bias =nanmean(current_baseline_data - current_obs_data);
    text(current_dates(10),3.2,sprintf('R^{2}=%.2f; NSE = %.2f; bias = %.2f',round(R^2,2),round(NSE,2),round(bias,2)),'fontsize',30,'color','b')
    % if p>0.001
    %     text(dates_prefire(1000),10,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias - %.2f',round(R^2,2),round(p,2),round(NSE,2),round(bias,2)),'fontsize',40,'color','b')
    % else
    %     text(dates_prefire(1000),10,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f; bias = %.2f',round(R^2,2),round(NSE,2),round(bias,2)),'fontsize',40,'color','b')
    % end
    xtickangle(35)
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_MFF_Qanom_timeseries_prefire_Standardized_similarDryperiod%d.png',k))
end
ppp
%% show time series for pre fire comparissons - monthly:
matched_datevecs = datevec(matched_dates);
[unique_dates , ~,j] = unique(matched_datevecs(:,1:2),'rows');
Baseline_Q_Mon = accumarray(j,Baseline_Q_daily_matched,[],@nanmean);
Realistic_Q_Mon = accumarray(j,Realistic_Q_daily_matched,[],@nanmean);
Obs_Q_Mon = accumarray(j,obs_Q_daily_matched,[],@nanmean);

Baseline_Q_monthly_anoms = (Baseline_Q_Mon - nanmean(Baseline_Q_Mon))./nanstd(Baseline_Q_Mon);
Realistic_Q_monthly_anoms = (Realistic_Q_Mon - nanmean(Realistic_Q_Mon))./nanstd(Realistic_Q_Mon);
Obs_Q_monthly_anoms = (Obs_Q_Mon - nanmean(Obs_Q_Mon))./nanstd(Obs_Q_Mon);

matched_dates_monthly = [unique_dates , ones(length(unique_dates),1)];
matched_dates_monthly = datenum(matched_dates_monthly);

idx_postfire = find(matched_dates_monthly >= post_fire_start);
idx_prefire = find(matched_dates_monthly < pre_fire_end);

dates_prefire= matched_dates_monthly(idx_prefire);
prefire_obs = Obs_Q_monthly_anoms(idx_prefire);
prefire_baseline = Baseline_Q_monthly_anoms(idx_prefire);

f=figure;
hold on
p1=plot(dates_prefire,prefire_obs,'-k','linewidth',3);
p2=plot(dates_prefire,prefire_baseline,'-b','linewidth',1.5);
grid on
box on
xlim([min(dates_prefire) max(dates_prefire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',40)
%ylabel('pre-fire Q anomalies')
% % legend([p1 p2],{'obs','baseline'},'location','best','fontsize',22)
f.Position = [-1803         325        1804         435];
%report R:
[R,p] = corr(prefire_obs,prefire_baseline,'type','pearson','rows','complete');
%report NSE:
NSE =nashsutcliffe([dates_prefire,prefire_baseline], [dates_prefire,prefire_obs]);
%report bias:
bias =nanmean(prefire_baseline - prefire_obs);
% if p>0.001
%     text(dates_prefire(50),6,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias - %.2f',round(R^2,2),round(p,2),round(NSE,2),round(bias,2)),'fontsize',40,'color','b')
% else
%     text(dates_prefire(50),6,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f; bias = %.2f',round(R^2,2),round(NSE,2),round(bias,2)),'fontsize',40,'color','b')
% end
saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_MFF_Qanom_timeseries_prefire_monthlyprefireSTD.png')

%% show time series for post fire comparissons - monthly:
dates_postfire= matched_dates_monthly(idx_postfire);
postfire_obs = Obs_Q_monthly_anoms(idx_postfire);
postfire_baseline = Baseline_Q_monthly_anoms(idx_postfire);
postfire_realistic = Realistic_Q_monthly_anoms(idx_postfire);

f=figure;
hold on
p1=plot(dates_postfire,postfire_obs,'-k','linewidth',3);
p2=plot(dates_postfire,postfire_baseline,'-b','linewidth',1.5);
p3=plot(dates_postfire,postfire_realistic,'-r','linewidth',1.5);

grid on
box on
xlim([min(dates_postfire) max(dates_postfire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',40)
%ylabel('post-fire Q anomalies')
% % leg=legend([p1 p2 p3],{'obs','baseline','Mod-params+GVF+Veg-class+Snow-alb'},'location','best','fontsize',22);
% % leg.Position = [0.6580    0.7219    0.2450    0.1908];
f.Position = [-1803         325         928         435];

%report R:
[R_baseline,p_baseline] = corr(postfire_obs,postfire_baseline,'type','pearson','rows','complete');
%report NSE:
NSE_baseline =nashsutcliffe([dates_postfire,postfire_baseline], [dates_postfire,postfire_obs]);
%report bias:
bias_baseline = nanmean(postfire_baseline - postfire_obs);
% if p>0.001
%     text(dates_postfire(6),0.5,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias = %.2f',round(R_baseline^2,2),round(p_baseline,2),round(NSE_baseline,2),round(bias_baseline,2)),'fontsize',20,'color','b')
% else
%     text(dates_postfire(6),0.5,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f; bias = %.2f',round(R_baseline^2,2),round(NSE_baseline,2),round(bias_baseline,2)),'fontsize',20,'color','b')
% end

%report R:
[R_realistic,p_realistic] = corr(postfire_obs,postfire_realistic,'type','pearson','rows','complete');
%report NSE:
NSE_realistic =nashsutcliffe([dates_postfire,postfire_realistic], [dates_postfire,postfire_obs]);
%report bias:
bias_realistic = nanmean(postfire_realistic - postfire_obs);
% if p>0.001
%     text(dates_postfire(6),0.3,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias = %.2f',round(R_realistic^2,2),round(p_realistic,2),round(NSE_realistic,2),round(bias_realistic,2)),'fontsize',20,'color','r')
% else
%     text(dates_postfire(6),0.3,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f; bias = %.2f',round(R_realistic^2,2),round(NSE_realistic,2),round(bias_realistic,2)),'fontsize',20,'color','r')
% end
saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/MFF_Qanom_timeseries_postfire_monthlyprefireSTD.png')