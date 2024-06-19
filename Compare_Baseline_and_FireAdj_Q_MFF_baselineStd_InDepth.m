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

%% Get SWE data for comparisson:
baseline_filename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/Catchment_averaged_outputs/baseline_LSM_outputs_Catch_2.mat');
realistic_filename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/Catchment_averaged_outputs/realistic_LSM_outputs_Catch_2.mat');
%baseline data
baseline_data = load(baseline_filename);
baseline_data = baseline_data.Catchment_outputs;
dates = baseline_data.date;
SWE_dates = double(dates);
baseline_SWE = baseline_data.swe;
baseline_soilm1 = baseline_data.soilm1;
baseline_soilm2 = baseline_data.soilm2;
baseline_soilm3 = baseline_data.soilm3;
baseline_soilm4 = baseline_data.soilm4;
baseline_soilm = baseline_soilm1.*(0.1/2) + baseline_soilm2.*(0.3/2) + baseline_soilm4.*(0.6/2) + baseline_soilm4.*(1/2);

%realistic data
realistic_data = load(realistic_filename);
realistic_data = realistic_data.Catchment_outputs;
realistic_SWE = realistic_data.swe;
realistic_soilm1 = realistic_data.soilm1;
realistic_soilm2 = realistic_data.soilm2;
realistic_soilm3 = realistic_data.soilm3;
realistic_soilm4 = realistic_data.soilm4;
realistic_soilm = realistic_soilm1.*(0.1/2) + realistic_soilm2.*(0.3/2) + realistic_soilm4.*(0.6/2) + realistic_soilm4.*(1/2);

%% load in simulation grid:
geofilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/Fire_Adjusted/wrfinput_Feather_postifre_Grassland.nc';
LAT = ncread(geofilename,'XLAT');
LON = ncread(geofilename,'XLONG');
%% load in shapefile data:
%load in study catchment shapes:
Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
Catchments = Catchments(2);
ncatch = length(Catchments);
%% for each catchment create indice lists for sim grid points in resepctive catchments
%show shape boundaries:
n = length(Catchments);
for i=1:n
    current_shape =  Catchments(i);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    in = inpolygon(LON(:) , LAT(:) , catch_lon , catch_lat);
    idx_in = find(in==1);
    
    %this variable stores all sim pixels in the catchment
    store_r_c = [];
    for j=1:length(idx_in)
        [r,c] = ind2sub(size(LAT) , idx_in(j));
        store_r_c = [store_r_c;r,c];
    end
    store_catch_idx_in = store_r_c;
end
met_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/Feather_Forcing/';
store_Catch_PRCP=[];
for WY=2000:2022
    current_dates = datenum([WY-1 10 1]):datenum([WY 9 30]);

    infilename = sprintf('daily_2D_WY%04d.nc',WY);
    PRCP = ncread([met_dir,infilename],'RAINRATE');
    PRCP = PRCP.*3600.*24;
    S=size(PRCP);
    PRCP = reshape(PRCP,[S(1)*S(2),S(3)]);
    PRCP = PRCP(store_catch_idx_in,:);
    mean_catch_PRCP = nanmean(PRCP,1)';
    store_Catch_PRCP = [store_Catch_PRCP;current_dates',mean_catch_PRCP];
end

%% match pre-fire data:
dates_prefire= matched_dates(idx_prefire);
[c1,ia1,ib1] = intersect(SWE_dates,dates_prefire);
[c2,ia2,ib2] = intersect(store_Catch_PRCP(:,1),dates_prefire);

Baseline_SWE_prefire = baseline_SWE(ia1);
Realistic_SWE_prefire = realistic_SWE(ia1);
PRCP_prefire = store_Catch_PRCP(ia2,2);
Realistic_soilm_prefire = realistic_soilm(ia1);
Baseline_soilm_prefire = baseline_soilm(ia1);

%% match post-fire data
dates_postfire= matched_dates(idx_postfire);
[c1,ia1,ib1] = intersect(SWE_dates,dates_postfire);
[c2,ia2,ib2] = intersect(store_Catch_PRCP(:,1),dates_postfire);

Baseline_SWE_postfire = baseline_SWE(ia1);
Realistic_SWE_postfire = realistic_SWE(ia1);
PRCP_postfire = store_Catch_PRCP(ia2,2);
Realistic_soilm_postfire = realistic_soilm(ia1);
Baseline_soilm_postfire = baseline_soilm(ia1);

%% show time series for post fire comparissons:
postfire_obs = Obs_Q_daily_anoms(idx_postfire);
postfire_baseline = Baseline_Q_daily_anoms(idx_postfire);
postfire_realistic = Realistic_Q_daily_anoms(idx_postfire);

f=figure;
subplot(3,1,1)
hold on

% % yyaxis left
left_color='k';
axR=gca;
axR.YColor=left_color;
p1=plot(dates_postfire,postfire_obs,'-k','linewidth',3);
p2=plot(dates_postfire,postfire_baseline,'-b','linewidth',1.5);
p3=plot(dates_postfire,postfire_realistic,'-r','linewidth',1.5);
ylabel('Q-anomaly')
set(gca,'fontsize',36)
ylim([-1 2])

%include hyetograph:
yyaxis right
right_color='b';
axR=gca;
axR.YColor=right_color;
prcp_plot = stem(dates_postfire,PRCP_postfire,'b');
set(gca,'ydir','reverse');
ylabel({'Prcp.'; '(mm/day)'},'fontsize',28)
ylim([0 185])
grid on
box on
xlim([min(dates_postfire) max(dates_postfire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',36)
%ylabel('post-fire Q anomalies')
% % legend([p1 p2 p3],{'obs','baseline','Mod-params+GVF+Veg-class+Snow-alb'},'location','best','fontsize',22)
f.Position = [-1919        -299        1811        1096];

%report R:
[R_baseline,p_baseline] = corr(postfire_obs,postfire_baseline,'type','pearson');
%report NSE:
NSE_baseline =nashsutcliffe([dates_postfire,postfire_baseline], [dates_postfire,postfire_obs]);
%report bias:
bias_baseline = nanmean(postfire_baseline - postfire_obs);
% if p>0.001
%     text(dates_postfire(100),3.5,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias = %.2f',round(R_baseline^2,2),round(p_baseline,2),round(NSE_baseline,2),round(bias_baseline,2)),'fontsize',40,'color','b')
% else
%     text(dates_postfire(100),3.5,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f;bias = %.2f',round(R_baseline^2,2),round(NSE_baseline,2),round(bias_baseline,2)),'fontsize',40,'color','b')
% end

%report R:
[R_realistic,p_realistic] = corr(postfire_obs,postfire_realistic,'type','pearson');
%report NSE:
NSE_realistic =nashsutcliffe([dates_postfire,postfire_realistic], [dates_postfire,postfire_obs]);
%report bias:
bias_realistic = nanmean(postfire_realistic - postfire_obs);
% if p>0.001
%     text(dates_postfire(100),2.8,sprintf('R^{2}=%.2f (p=%.2f); NSE = %.2f; bias = %.2f',round(R_realistic^2,2),round(p_realistic,2),round(NSE_realistic,2),round(bias_realistic,2)),'fontsize',40,'color','r')
% else
%     text(dates_postfire(100),2.8,sprintf('R^{2}=%.2f (p<0.001); NSE = %.2f; bias = %.2f',round(R_realistic^2,2),round(NSE_realistic,2),round(bias_realistic,2)),'fontsize',40,'color','r')
% end
xtickangle(35)
% saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/NFF_Qanom_timeseries_postfire_baseline_Standardized.png')

mean_postfire_anom = nanmean(postfire_obs);
max_postfire_anom = max(postfire_obs);

subplot(3,1,2)
hold on
p4 = plot(dates_postfire,Realistic_SWE_postfire,'-r');
p5 = plot(dates_postfire,Baseline_SWE_postfire,'-b');
ylabel('SWE (mm)')
grid on
box on
xlim([min(dates_postfire) max(dates_postfire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',36)

subplot(3,1,3)
hold on
p4 = plot(dates_postfire,Realistic_soilm_postfire,'-r');
p5 = plot(dates_postfire,Baseline_soilm_postfire,'-b');
ylabel({'soil moisture'; '(mm/mm)'})
grid on
box on
xlim([min(dates_postfire) max(dates_postfire)])
datetick('x','mm/yyyy','keepticks','keeplimits')
set(gca,'fontsize',36)
saveas(f,'/Users/abolafia/ASO_Fire/Plots/PaperPlots/MFF_Qanom_timeseries_postfire_baseline_Standardized_InDepth.png')

