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
obs_Q_daily_matched = store_obs_Q(ia,2);

%% convert to anomalies:
post_fire_start = datenum([2021 11 1]);
pre_fire_end = datenum([2020 8 1]);
idx_prefire = find(matched_dates < pre_fire_end);
Baseline_Q_daily_matched = Baseline_Q_daily_matched(idx_prefire);
Realistic_Q_daily_matched = Realistic_Q_daily_matched(idx_prefire);
obs_Q_daily_matched = obs_Q_daily_matched(idx_prefire);

Baseline_Q_daily_anoms = (Baseline_Q_daily_matched - nanmean(Baseline_Q_daily_matched))./nanstd(Baseline_Q_daily_matched);
Realistic_Q_daily_anoms = (Realistic_Q_daily_matched - nanmean(Realistic_Q_daily_matched))./nanstd(Realistic_Q_daily_matched);
Obs_Q_daily_anoms = (obs_Q_daily_matched - nanmean(obs_Q_daily_matched))./nanstd(obs_Q_daily_matched);

%% show time series for pre fire comparissons - not anomalies:
dates_prefire= matched_dates(idx_prefire);

f=figure;
hold on
p1=plot(dates_prefire,obs_Q_daily_matched.*cms_to_cfs_conversion_factor,'-k','linewidth',3);
p2=plot(dates_prefire,Baseline_Q_daily_matched.*cms_to_cfs_conversion_factor,'-b','linewidth',1.5);
p3=plot(dates_prefire,Realistic_Q_daily_matched.*cms_to_cfs_conversion_factor,'-r','linewidth',1.5);

grid on
box on
xlim([min(dates_prefire) max(dates_prefire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',40)
%ylabel('pre-fire Q anomalies')
% % legend([p1 p2],{'obs','baseline'},'location','best','fontsize',22)
f.Position = [-1803         325        1804         435];
xtickangle(35)

%report R:
[R_baseline,p_baseline] = corr(obs_Q_daily_matched,Baseline_Q_daily_matched,'type','pearson');
[R_realistic,p_realistic] = corr(obs_Q_daily_matched,Realistic_Q_daily_matched,'type','pearson');
%report NSE:
NSE_baseline =nashsutcliffe([dates_prefire,Baseline_Q_daily_matched], [dates_prefire,obs_Q_daily_matched]);
NSE_realistic =nashsutcliffe([dates_prefire,Realistic_Q_daily_matched], [dates_prefire,obs_Q_daily_matched]);
saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_MFF_Q_timeseries_prefire_eval.png')

%% use permutation test for 2 independent correlations:
delta_r = R_realistic - R_baseline;
n_perms = 10000;

L_baseline = length(Baseline_Q_daily_matched);
store_delta_r=[];
for j=1:n_perms
    permuted_baseline = Baseline_Q_daily_matched(randperm(L_baseline));
    permuted_enhanced = Baseline_Q_daily_matched(randperm(L_baseline));

    baseline_r_tmp = corr(permuted_baseline,obs_Q_daily_matched,'type','pearson');
    realistic_r_tmp = corr(permuted_enhanced,obs_Q_daily_matched,'type','pearson');
    
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
    permuted_baseline = Baseline_Q_daily_matched(randperm(L_baseline));
    permuted_enhanced = Baseline_Q_daily_matched(randperm(L_baseline));

    baseline_nse_tmp = nashsutcliffe([dates_prefire,permuted_baseline],[dates_prefire,obs_Q_daily_matched]);
    realistic_nse_tmp = nashsutcliffe([dates_prefire,permuted_enhanced],[dates_prefire,obs_Q_daily_matched]);
    
    store_delta_nse = [store_delta_nse; realistic_nse_tmp - baseline_nse_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_nse) > abs(delta_nse));
sprintf('difference in nse is significant with p = %.3f for daily validation',length(idx)/n_perms)


%% show time series for pre fire comparissons - anomalies:
dates_prefire= matched_dates(idx_prefire);

f=figure;
hold on
p1=plot(dates_prefire,Obs_Q_daily_anoms,'-k','linewidth',3);
p2=plot(dates_prefire,Baseline_Q_daily_anoms,'-b','linewidth',1.5);
p3=plot(dates_prefire,Realistic_Q_daily_anoms,'-r','linewidth',1.5);

grid on
box on
xlim([min(dates_prefire) max(dates_prefire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',40)
%ylabel('pre-fire Q anomalies')
% % legend([p1 p2],{'obs','baseline'},'location','best','fontsize',22)
f.Position = [-1803         325        1804         435];
xtickangle(35)

%report R:
[R_baseline_anom,p_baseline_anom] = corr(Obs_Q_daily_anoms,Baseline_Q_daily_anoms,'type','pearson');
[R_realistic_anom,p_realistic_anom] = corr(Obs_Q_daily_anoms,Realistic_Q_daily_anoms,'type','pearson');
%report NSE:
NSE_baseline_anom =nashsutcliffe([dates_prefire,Baseline_Q_daily_anoms], [dates_prefire,Obs_Q_daily_anoms]);
NSE_realistic_anom =nashsutcliffe([dates_prefire,Realistic_Q_daily_anoms], [dates_prefire,Obs_Q_daily_anoms]);

saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_MFF_Qanom_timeseries_prefire_eval.png')

%% use permutation test for 2 independent correlations:
delta_r = R_realistic_anom - R_baseline_anom;
n_perms = 10000;

L_baseline = length(Baseline_Q_daily_anoms);
store_delta_r=[];
for j=1:n_perms
    permuted_baseline = Baseline_Q_daily_anoms(randperm(L_baseline));
    permuted_enhanced = Baseline_Q_daily_anoms(randperm(L_baseline));

    baseline_r_tmp = corr(permuted_baseline,Obs_Q_daily_anoms,'type','pearson');
    realistic_r_tmp = corr(permuted_enhanced,Obs_Q_daily_anoms,'type','pearson');
    
    store_delta_r = [store_delta_r; realistic_r_tmp - baseline_r_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_r) > abs(delta_r));
sprintf('difference in r anomalies is significant with p = %.3f for daily validation',length(idx)/n_perms)

%% use permutation test for 2 independent NSE's:
delta_nse = NSE_realistic_anom - NSE_baseline_anom;
n_perms = 10000;

store_delta_nse=[];
for j=1:n_perms
    permuted_baseline = Baseline_Q_daily_anoms(randperm(L_baseline));
    permuted_enhanced = Baseline_Q_daily_anoms(randperm(L_baseline));

    baseline_nse_tmp = nashsutcliffe([dates_prefire,permuted_baseline],[dates_prefire,Obs_Q_daily_anoms]);
    realistic_nse_tmp = nashsutcliffe([dates_prefire,permuted_enhanced],[dates_prefire,Obs_Q_daily_anoms]);
    
    store_delta_nse = [store_delta_nse; realistic_nse_tmp - baseline_nse_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_nse) > abs(delta_nse));
sprintf('difference in nse anomalies is significant with p = %.3f for daily validation',length(idx)/n_perms)

%% show time series for pre fire comparissons - monthly - not anomalies:
matched_datevecs = datevec(dates_prefire);
[unique_dates , ~,j] = unique(matched_datevecs(:,1:2),'rows');
Baseline_Q_Mon = accumarray(j,Baseline_Q_daily_matched,[],@nanmean);
Realistic_Q_Mon = accumarray(j,Realistic_Q_daily_matched,[],@nanmean);
Obs_Q_Mon = accumarray(j,obs_Q_daily_matched,[],@nanmean);

matched_dates_monthly = [unique_dates , ones(length(unique_dates),1)];
matched_dates_monthly = datenum(matched_dates_monthly);

Baseline_Q_monthly_anoms = (Baseline_Q_Mon - nanmean(Baseline_Q_Mon))./nanstd(Baseline_Q_Mon);
Realistic_Q_monthly_anoms = (Realistic_Q_Mon - nanmean(Realistic_Q_Mon))./nanstd(Realistic_Q_Mon);
Obs_Q_monthly_anoms = (Obs_Q_Mon - nanmean(Obs_Q_Mon))./nanstd(Obs_Q_Mon);

f=figure;
hold on
p1=plot(matched_dates_monthly,Obs_Q_Mon.*cms_to_cfs_conversion_factor,'-k','linewidth',3);
p2=plot(matched_dates_monthly,Baseline_Q_Mon.*cms_to_cfs_conversion_factor,'-b','linewidth',1.5);
p3=plot(matched_dates_monthly,Realistic_Q_Mon.*cms_to_cfs_conversion_factor,'-r','linewidth',1.5);
grid on
box on
xlim([min(matched_dates_monthly) max(matched_dates_monthly)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',40)
f.Position = [-1803         325        1804         435];
%report R:
[R_baseline_monthly,p_baseline_monthly] = corr(Obs_Q_Mon,Baseline_Q_Mon,'type','pearson');
[R_realistic_monthly,p_realistic_monthly] = corr(Obs_Q_Mon,Realistic_Q_Mon,'type','pearson');
%report NSE:
NSE_baseline_monthly =nashsutcliffe([matched_dates_monthly,Baseline_Q_Mon], [matched_dates_monthly,Obs_Q_Mon]);
NSE_realistic_monthly =nashsutcliffe([matched_dates_monthly,Realistic_Q_Mon], [matched_dates_monthly,Obs_Q_Mon]);

saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_MFF_Q_timeseries_prefire_monthly_Standardized_eval.png')


%% use permutation test for 2 independent correlations:
delta_r = R_realistic_monthly - R_baseline_monthly;
n_perms = 10000;

L_baseline = length(Baseline_Q_Mon);
store_delta_r=[];
for j=1:n_perms
    permuted_baseline = Baseline_Q_Mon(randperm(L_baseline));
    permuted_enhanced = Baseline_Q_Mon(randperm(L_baseline));

    baseline_r_tmp = corr(permuted_baseline,Obs_Q_Mon,'type','pearson');
    realistic_r_tmp = corr(permuted_enhanced,Obs_Q_Mon,'type','pearson');
    
    store_delta_r = [store_delta_r; realistic_r_tmp - baseline_r_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_r) > abs(delta_r));
sprintf('MONTHLY: difference in r is significant with p = %.3f for daily validation',length(idx)/n_perms)

%% use permutation test for 2 independent NSE's:
delta_nse = NSE_realistic_monthly - NSE_baseline_monthly;
n_perms = 10000;
L_baseline = length(Baseline_Q_Mon);
store_delta_nse=[];
for j=1:n_perms
    permuted_baseline = Baseline_Q_Mon(randperm(L_baseline));
    permuted_enhanced = Baseline_Q_Mon(randperm(L_baseline));

    baseline_nse_tmp = nashsutcliffe([matched_dates_monthly,permuted_baseline],[matched_dates_monthly,Obs_Q_Mon]);
    realistic_nse_tmp = nashsutcliffe([matched_dates_monthly,permuted_enhanced],[matched_dates_monthly,Obs_Q_Mon]);
    
    store_delta_nse = [store_delta_nse; realistic_nse_tmp - baseline_nse_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_nse) > abs(delta_nse));
sprintf('MONTHLY: difference in nse is significant with p = %.3f for daily validation',length(idx)/n_perms)

%% show time series for pre fire comparissons - monthly -  anomalies:
f=figure;
hold on
p1=plot(matched_dates_monthly,Obs_Q_monthly_anoms,'-k','linewidth',3);
p2=plot(matched_dates_monthly,Baseline_Q_monthly_anoms,'-b','linewidth',1.5);
p3=plot(matched_dates_monthly,Realistic_Q_monthly_anoms,'-r','linewidth',1.5);
grid on
box on
xlim([min(matched_dates_monthly) max(matched_dates_monthly)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',40)
f.Position = [-1803         325        1804         435];
%report R:
[R_baseline_monthly_anom,p_baseline_monthly_anom] = corr(Obs_Q_monthly_anoms,Baseline_Q_monthly_anoms,'type','pearson');
[R_realistic_monthly_anom,p_realistic_monthly_anom] = corr(Obs_Q_monthly_anoms,Realistic_Q_monthly_anoms,'type','pearson');
%report NSE:
NSE_baseline_monthly_anom =nashsutcliffe([matched_dates_monthly,Baseline_Q_monthly_anoms], [matched_dates_monthly,Obs_Q_monthly_anoms]);
NSE_realistic_monthly_anom =nashsutcliffe([matched_dates_monthly,Realistic_Q_monthly_anoms], [matched_dates_monthly,Obs_Q_monthly_anoms]);
saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_MFF_Qanom_timeseries_prefire_monthly_Standardized_eval.png')


%% use permutation test for 2 independent correlations:
delta_r_anom = R_realistic_monthly_anom - R_baseline_monthly_anom;
n_perms = 10000;

L_baseline = length(Baseline_Q_monthly_anoms);
store_delta_r=[];
for j=1:n_perms
    permuted_baseline = Baseline_Q_monthly_anoms(randperm(L_baseline));
    permuted_enhanced = Baseline_Q_monthly_anoms(randperm(L_baseline));

    baseline_r_tmp = corr(permuted_baseline,Obs_Q_monthly_anoms,'type','pearson');
    realistic_r_tmp = corr(permuted_enhanced,Obs_Q_monthly_anoms,'type','pearson');
    
    store_delta_r = [store_delta_r; realistic_r_tmp - baseline_r_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_r) > abs(delta_r_anom));
sprintf('MONTHLY anom: difference in r is significant with p = %.3f for daily validation',length(idx)/n_perms)

%% use permutation test for 2 independent NSE's:
delta_nse_anom = NSE_realistic_monthly_anom - NSE_baseline_monthly_anom;
n_perms = 10000;
L_baseline = length(Baseline_Q_monthly_anoms);
store_delta_nse=[];
for j=1:n_perms
    permuted_baseline = Baseline_Q_monthly_anoms(randperm(L_baseline));
    permuted_enhanced = Baseline_Q_monthly_anoms(randperm(L_baseline));

    baseline_nse_tmp = nashsutcliffe([matched_dates_monthly,permuted_baseline],[matched_dates_monthly,Obs_Q_monthly_anoms]);
    realistic_nse_tmp = nashsutcliffe([matched_dates_monthly,permuted_enhanced],[matched_dates_monthly,Obs_Q_monthly_anoms]);
    
    store_delta_nse = [store_delta_nse; realistic_nse_tmp - baseline_nse_tmp];
end
%what percent of the time is the difference bigger by random selection than in reality?
idx=find(abs(store_delta_nse) > abs(delta_nse_anom));
sprintf('MONTHLY anom: difference in nse is significant with p = %.3f for daily validation',length(idx)/n_perms)