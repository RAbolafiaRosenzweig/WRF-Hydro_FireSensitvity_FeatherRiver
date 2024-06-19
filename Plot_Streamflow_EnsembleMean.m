clc;clear all;close all;

cms_to_cfs_conversion_factor = 35.314666212661;

ex_dates = datenum([2005 10 1]):datenum([2006 9 30]);
idx_season = find(ex_dates>=datenum([2005 12 1]) & ex_dates<=datenum([2006 3 31]));

met_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/Feather_Forcing/';

WYs = 2000:2022;
nyears = length(WYs);

%% plot hydrographs during these years for: MIDDLE FORK
STA_IDs = {'DCS','F56','F57','GYB','ICR','JBR','LCB','MER','MFP','NFP','NYS','ORH','SPK','TM1','TM2','TM3','WFR','YPB','YRS'};
% % stations = {'MER','MFP'}; %stations 
stations = {'MER'}; %station closest to outlet of middle fork feather

store_sta_idx=[];
for s=1:length(stations)
    IDX_sta = strcmp(STA_IDs,stations{1});
    idx = find(IDX_sta==1);
    store_sta_idx = [store_sta_idx;idx];
end

%load in simulated Q - baseline:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Feather_Baseline/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
Baseline_Streamflow = Simulated_Q.Streamflow;
Baseline_Streamflow = Baseline_Streamflow(:,store_sta_idx);
Baseline_Streamflow = nanmean(Baseline_Streamflow,2);
sim_dates = Simulated_Q.dates;
sim_datevecs = datevec(sim_dates);

%load in simulated Q - modified parameters
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParams_Streamflow = Simulated_Q.Streamflow;
ModParams_Streamflow = ModParams_Streamflow(:,store_sta_idx);
ModParams_Streamflow = nanmean(ModParams_Streamflow,2);
%load in simulated Q - modified parameters & veg class - BARE:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam_GVF/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParam_GVF = Simulated_Q.Streamflow;
ModParam_GVF = ModParam_GVF(:,store_sta_idx);
ModParam_GVF_streamflow = nanmean(ModParam_GVF,2);
%load in simulated Q - modified parameters & veg class - GRASS:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam_GVF_VegClass/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParam_GVF_VegClass = Simulated_Q.Streamflow;
ModParam_GVF_VegClass = ModParam_GVF_VegClass(:,store_sta_idx);
ModParam_GVF_VegClass_streamflow = nanmean(ModParam_GVF_VegClass,2);

%load in simulated Q - modified parameters & veg class & albedo - BARE:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Realistic/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
Realistic = Simulated_Q.Streamflow;
Realistic = Realistic(:,store_sta_idx);
Realistic_Streamflow = nanmean(Realistic,2);

%compute daily mean hydrographs 
nyears = length(WYs);
store_baseline_Q=[];
store_ModParam_Q=[];
store_ModParam_GVF_Q=[];
store_ModParam_GVF_VegClass_Q=[];
store_Realistic_Q=[];

for y = 1:nyears
    current_year =  WYs(y);
    current_date_idx = find(sim_dates>=datenum([current_year-1 10 1]) & sim_dates<=datenum([current_year 9 30]));
    current_datevecs = sim_datevecs(current_date_idx,:);
    
    baseline_Q =Baseline_Streamflow(current_date_idx);
    mod_param_Q =ModParams_Streamflow(current_date_idx);
    mod_param_GVF_Q =ModParam_GVF_streamflow(current_date_idx);
    mod_param_GVF_VegClass_Q =ModParam_GVF_VegClass_streamflow(current_date_idx);
    Realistic_Q = Realistic_Streamflow(current_date_idx);
    
    %compute daily mean:
    [u,~,j] = unique(current_datevecs(:,1:3),'rows','stable');
    
    baseline_Q = accumarray(j,baseline_Q,[],@nanmean);
    mod_param_Q = accumarray(j,mod_param_Q,[],@nanmean);
    mod_param_GVF_Q = accumarray(j,mod_param_GVF_Q,[],@nanmean);
    mod_param_GVF_VegClass_Q = accumarray(j,mod_param_GVF_VegClass_Q,[],@nanmean);
    Realistic_Q = accumarray(j,Realistic_Q,[],@nanmean);
    
    store_baseline_Q = [store_baseline_Q,baseline_Q(1:365)];
    store_ModParam_Q = [store_ModParam_Q,mod_param_Q(1:365)];
    store_ModParam_GVF_Q = [store_ModParam_GVF_Q,mod_param_GVF_Q(1:365)];
    store_ModParam_GVF_VegClass_Q = [store_ModParam_GVF_VegClass_Q,mod_param_GVF_VegClass_Q(1:365)];
    store_Realistic_Q = [store_Realistic_Q,Realistic_Q(1:365)];
end

store_baseline_Q = store_baseline_Q.*cms_to_cfs_conversion_factor;
store_ModParam_Q = store_ModParam_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_Q = store_ModParam_GVF_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_VegClass_Q = store_ModParam_GVF_VegClass_Q.*cms_to_cfs_conversion_factor;
store_Realistic_Q = store_Realistic_Q.*cms_to_cfs_conversion_factor;

baseline_ann_Q = nanmean(store_baseline_Q,1);
ModParam_ann_Q = nanmean(store_ModParam_Q,1);
ModParam_GVF_ann_Q = nanmean(store_ModParam_GVF_Q,1);
ModParam_GVF_vegclass_ann_Q = nanmean(store_ModParam_GVF_VegClass_Q,1);
realistic_ann_Q = nanmean(store_Realistic_Q,1);

delta_baseline = (baseline_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_realistic = (realistic_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam = (ModParam_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF = (ModParam_GVF_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF_vegclass = (ModParam_GVF_vegclass_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF_vegclass2 = (ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q)./ModParam_GVF_ann_Q;
delta_albedo = (realistic_ann_Q - ModParam_GVF_vegclass_ann_Q)./ModParam_GVF_vegclass_ann_Q;

sprintf('MFF: baseline result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_baseline)*100,min(delta_baseline)*100,max(delta_baseline)*100)
sprintf('MFF: ModParam  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam )*100,min(delta_ModParam )*100,max(delta_ModParam )*100)
sprintf('MFF: ModParam_GVF  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF )*100,min(delta_ModParam_GVF )*100,max(delta_ModParam_GVF )*100)
sprintf('MFF: ModParam_GVF_vegclass  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF_vegclass )*100,min(delta_ModParam_GVF_vegclass )*100,max(delta_ModParam_GVF_vegclass )*100)
sprintf('MFF: realistic result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_realistic)*100,min(delta_realistic)*100,max(delta_realistic)*100)
sprintf('MFF: veg class result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF_vegclass2)*100,min(delta_ModParam_GVF_vegclass2)*100,max(delta_ModParam_GVF_vegclass2)*100)
sprintf('MFF: albedo result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_albedo)*100,min(delta_albedo)*100,max(delta_albedo)*100)

%daily mod-param - baseline deltas
mean_baseline_Q = nanmean(store_baseline_Q(:));
delta_ModParam_daily_Q = (store_ModParam_Q - store_baseline_Q);
mean_rel_diff_ModParam_daily_Q = nanmean(delta_ModParam_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('MFF: DAILY mean diff from mod params is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_daily_Q*100,nanmean(delta_ModParam_daily_Q(:)), max(delta_ModParam_daily_Q(:)), length(idx))

%daily ModParam_GVF_ann_Q - baseline deltas
delta_ModParam_GVF_daily_Q = (store_ModParam_GVF_Q - store_baseline_Q);
mean_rel_diff_ModParam_GVF_daily_Q = nanmean(delta_ModParam_GVF_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_GVF_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('MFF: DAILY mean diff from modparams + GVF is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_GVF_daily_Q*100,nanmean(delta_ModParam_GVF_daily_Q(:)), max(delta_ModParam_GVF_daily_Q(:)), length(idx))

%daily ModParam_GVF_vegclass_ann_Q - baseline deltas
delta_ModParam_GVF_vegclass_daily_Q = (store_ModParam_GVF_VegClass_Q - store_baseline_Q);
mean_rel_diff_ModParam_GVF_vegclass_daily_Q = nanmean(delta_ModParam_GVF_vegclass_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_GVF_vegclass_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('MFF: DAILY mean diff from modparams + GVF +veg class is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_GVF_vegclass_daily_Q*100,nanmean(delta_ModParam_GVF_vegclass_daily_Q(:)), max(delta_ModParam_GVF_vegclass_daily_Q(:)), length(idx))

%compare winter ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q
idx_winter = 93:182;
store_ModParam_GVF_Q_winter = store_ModParam_GVF_Q(idx_winter,:);
store_ModParam_GVF_VegClass_Q_winter = store_ModParam_GVF_VegClass_Q(idx_winter,:);
pct_change_winter =  (store_ModParam_GVF_VegClass_Q_winter - store_ModParam_GVF_Q_winter)./store_ModParam_GVF_Q_winter;
mean_pct_change_winter = nanmean(pct_change_winter(:));
sprintf('MFF: DAILY mean diff from modparams + GVF +veg class vs GVF in winter is %.4f%%',mean_pct_change_winter*100)

%compare summer ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q
idx_summer = 213:273;
store_ModParam_GVF_Q_summer = store_ModParam_GVF_Q(idx_summer,:);
store_ModParam_GVF_VegClass_Q_summer = store_ModParam_GVF_VegClass_Q(idx_summer,:);
pct_change_summer =  (store_ModParam_GVF_VegClass_Q_summer - store_ModParam_GVF_Q_summer)./store_ModParam_GVF_Q_summer;
mean_pct_change_summer = nanmean(pct_change_summer(:));
sprintf('MFF: DAILY mean diff from modparams + GVF +veg class vs GVF in summer is %.4f%%',mean_pct_change_summer*100)

%compare winter realistic_ann_Q_winter - store_ModParam_GVF_VegClass_Q_winter
idx_winter = 93:182;
realistic_ann_Q_winter = store_Realistic_Q(idx_winter,:);
store_ModParam_GVF_VegClass_Q_winter = store_ModParam_GVF_VegClass_Q(idx_winter,:);
pct_change_winter =  (realistic_ann_Q_winter - store_ModParam_GVF_VegClass_Q_winter)./store_ModParam_GVF_VegClass_Q_winter;
mean_pct_change_winter = nanmean(pct_change_winter(:));
sprintf('MFF: DAILY mean diff from modparams + GVF +veg class +albedo vs GVF+class  in winter is %.4f%%',mean_pct_change_winter*100)

%compare summer realistic_ann_Q_summer - store_ModParam_GVF_VegClass_Q_summer
idx_summer = 183:273;
realistic_ann_Q_summer = store_Realistic_Q(idx_summer,:);
store_ModParam_GVF_VegClass_Q_summer = store_ModParam_GVF_VegClass_Q(idx_summer,:);
pct_change_summer =  (realistic_ann_Q_summer - store_ModParam_GVF_VegClass_Q_summer)./store_ModParam_GVF_VegClass_Q_summer;
mean_pct_change_summer = nanmean(pct_change_summer(:));
sprintf('MFF: DAILY mean diff from modparams + GVF +veg class +albedo vs GVF+class  in summer is %.4f%%',mean_pct_change_summer*100)

%compare pt1 realistic_ann_Q_pt1 - baseline_pt1
idx_pt1 = 1:212;
realistic_ann_Q_pt1 = store_Realistic_Q(idx_pt1,:);
store_baseline_Q_pt1 = store_baseline_Q(idx_pt1,:);
pct_change_pt1 =  (realistic_ann_Q_pt1 - store_baseline_Q_pt1)./store_baseline_Q_pt1;
mean_pct_change_pt1 = nanmean(pct_change_pt1(:));
sprintf('MFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in pt1 is %.4f%%',mean_pct_change_pt1*100)

%compare pt2 realistic_ann_Q_pt2 - baseline_pt2
idx_pt2 = 244:365;
realistic_ann_Q_pt2 = store_Realistic_Q(idx_pt2,:);
store_baseline_Q_pt2 = store_baseline_Q(idx_pt2,:);
pct_change_pt2 =  (realistic_ann_Q_pt2 - store_baseline_Q_pt2)./store_baseline_Q_pt2;
mean_pct_change_pt2 = nanmean(pct_change_pt2(:));
sprintf('MFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in pt2 is %.4f%%',mean_pct_change_pt2*100)

%compare may realistic_ann_Q_pt2 - baseline_pt2
idx_pt2 = 213:243;
realistic_ann_Q_pt2 = store_Realistic_Q(idx_pt2,:);
store_baseline_Q_pt2 = store_baseline_Q(idx_pt2,:);
pct_change_pt2 =  (realistic_ann_Q_pt2 - store_baseline_Q_pt2)./store_baseline_Q_pt2;
mean_pct_change_pt2 = nanmean(pct_change_pt2(:));
sprintf('MFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in May is %.4f%%',mean_pct_change_pt2*100)


baseline_hydrograph = nanmean(store_baseline_Q,2)';
mod_param_hydrograph = nanmean(store_ModParam_Q,2)';
mod_param_GVF_hydrograph = nanmean(store_ModParam_GVF_Q,2)';
mod_param_GVF_VegClass_hydrograph = nanmean(store_ModParam_GVF_VegClass_Q,2)';
Realistic_hydrograph = nanmean(store_Realistic_Q,2)';

%compute corresponding STDs:
baseline_hydrograph_std = std(store_baseline_Q','omitnan');
mod_param_hydrograph_std = std(store_ModParam_Q','omitnan');
store_ModParam_GVF_Q_std = std(store_ModParam_GVF_Q','omitnan');
store_ModParam_GVF_VegClass_Q_std = std(store_ModParam_GVF_VegClass_Q','omitnan');
Realistic_hydrograph_std = std(store_Realistic_Q','omitnan');

f=figure;
f.Position = [ -1535         -48        1139        1096];
subplot(2,1,1)
%title('Middle Fork (station MER)','fontsize',42)
hold on
%Plot basecase 
tpday = 365;
p01=plot(1:tpday,baseline_hydrograph(1:tpday),'-k','linewidth',3);
%Plot mod param case 
p02=plot(1:tpday,mod_param_hydrograph(1:tpday),'-b','linewidth',1.2);
%Plot mod param+gvf
p03=plot(1:tpday,mod_param_GVF_hydrograph(1:tpday),'-','linewidth',1.2,'color',[0 0.5 0]);
%Plot mod param+gvf+veg class
p04=plot(1:tpday,mod_param_GVF_VegClass_hydrograph(1:tpday),'--','linewidth',1.2,'color',[0 0.5 0]);
%realistic
p05=plot(1:tpday,Realistic_hydrograph(1:tpday),'-','linewidth',1.2,'color','r');
grid on
box on
ylabel({'multi-year mean daily'; 'Q (cfs) (2000-2022)'},'fontsize',42)
set(gca,'fontsize',42)
xlabel('day of water year','fontsize',42)
xlim([1 365])
ppp
%delta baseline plot:
subplot(2,1,2)
p01=plot(-1.*(1:tpday),baseline_hydrograph(1:tpday),'-k','linewidth',3); %for legend
hold on
%Plot mod param case 
p02=plot(1:tpday,mod_param_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-b','linewidth',2);
%Plot mod param+gvf
p03=plot(1:tpday,mod_param_GVF_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-','linewidth',1.2,'color',[0 0.5 0]);
%Plot mod param+gvf+veg class
p04=plot(1:tpday,mod_param_GVF_VegClass_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'--','linewidth',1.2,'color',[0 0.5 0]);
%realistic
p05=plot(1:tpday,Realistic_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-','linewidth',1.2,'color','r');
grid on
box on
ylabel({'\Delta multi-year daily Q'; '(cfs) (2000-2022)'},'fontsize',42)
set(gca,'fontsize',42)
xlabel('day of water year','fontsize',42)
xlim([1 365]);xticks([1:50:365]);xtickangle(90);
ylim([-200 800])
% leg = legend([p01 p02 p03 p04 p05],{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'},'fontsize',34,'location','best');
% leg.Position = [0.3738    0.3848    0.5755    0.2008];
saveas(f,'/Users/abolafia/ASO_Fire/Plots/Q_enmsemble_mean_2000_2022_MiddleFork.png')

%include bar plot comparing annual means:
mean_Q = [nanmean(baseline_hydrograph) , nanmean(mod_param_hydrograph) , nanmean(mod_param_GVF_hydrograph),nanmean(mod_param_GVF_VegClass_hydrograph) ,nanmean(Realistic_hydrograph)];
X=categorical({'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'});
X = reordercats(X,{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'});
f=figure;
b=bar(X,mean_Q,'FaceColor','k');
grid on
box on
set(gca,'fontsize',20)
f.Position = [-1536         232         643         515];
ylabel('mean water year Q(cfs)','fontsize',28)
saveas(f,'/Users/abolafia/ASO_Fire/Plots/Q_ensemble_ANN_mean_barplot_MiddleFork.png')

%% plot hydrographs during these years for: East Branch North FORK feather
STA_IDs = {'DCS','F56','F57','GYB','ICR','JBR','LCB','MER','MFP','NFP','NYS','ORH','SPK','TM1','TM2','TM3','WFR','YPB','YRS'};
stations = {'ICR'}; 

store_sta_idx=[];
for s=1:length(stations)
    IDX_sta = strcmp(STA_IDs,stations{1});
    idx = find(IDX_sta==1);
    store_sta_idx = [store_sta_idx;idx];
end

%load in simulated Q - baseline:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Feather_Baseline/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
Baseline_Streamflow = Simulated_Q.Streamflow;
Baseline_Streamflow = Baseline_Streamflow(:,store_sta_idx);
Baseline_Streamflow = nanmean(Baseline_Streamflow,2);
sim_dates = Simulated_Q.dates;
sim_datevecs = datevec(sim_dates);

%load in simulated Q - modified parameters
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParams_Streamflow = Simulated_Q.Streamflow;
ModParams_Streamflow = ModParams_Streamflow(:,store_sta_idx);
ModParams_Streamflow = nanmean(ModParams_Streamflow,2);
%load in simulated Q - modified parameters & veg class - BARE:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam_GVF/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParam_GVF = Simulated_Q.Streamflow;
ModParam_GVF = ModParam_GVF(:,store_sta_idx);
ModParam_GVF_streamflow = nanmean(ModParam_GVF,2);
%load in simulated Q - modified parameters & veg class - GRASS:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam_GVF_VegClass/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParam_GVF_VegClass = Simulated_Q.Streamflow;
ModParam_GVF_VegClass = ModParam_GVF_VegClass(:,store_sta_idx);
ModParam_GVF_VegClass_streamflow = nanmean(ModParam_GVF_VegClass,2);

%load in simulated Q - modified parameters & veg class & albedo - BARE:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Realistic/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
Realistic = Simulated_Q.Streamflow;
Realistic = Realistic(:,store_sta_idx);
Realistic_Streamflow = nanmean(Realistic,2);

%compute daily mean hydrographs 
nyears = length(WYs);
store_baseline_Q=[];
store_ModParam_Q=[];
store_ModParam_GVF_Q=[];
store_ModParam_GVF_VegClass_Q=[];
store_Realistic_Q=[];

for y = 1:nyears
    current_year =  WYs(y);
    current_date_idx = find(sim_dates>=datenum([current_year-1 10 1]) & sim_dates<=datenum([current_year 9 30]));
    current_datevecs = sim_datevecs(current_date_idx,:);
    
    baseline_Q =Baseline_Streamflow(current_date_idx);
    mod_param_Q =ModParams_Streamflow(current_date_idx);
    mod_param_GVF_Q =ModParam_GVF_streamflow(current_date_idx);
    mod_param_GVF_VegClass_Q =ModParam_GVF_VegClass_streamflow(current_date_idx);
    Realistic_Q = Realistic_Streamflow(current_date_idx);
    
    %compute daily mean:
    [u,~,j] = unique(current_datevecs(:,1:3),'rows','stable');
    
    baseline_Q = accumarray(j,baseline_Q,[],@nanmean);
    mod_param_Q = accumarray(j,mod_param_Q,[],@nanmean);
    mod_param_GVF_Q = accumarray(j,mod_param_GVF_Q,[],@nanmean);
    mod_param_GVF_VegClass_Q = accumarray(j,mod_param_GVF_VegClass_Q,[],@nanmean);
    Realistic_Q = accumarray(j,Realistic_Q,[],@nanmean);
    
    store_baseline_Q = [store_baseline_Q,baseline_Q(1:365)];
    store_ModParam_Q = [store_ModParam_Q,mod_param_Q(1:365)];
    store_ModParam_GVF_Q = [store_ModParam_GVF_Q,mod_param_GVF_Q(1:365)];
    store_ModParam_GVF_VegClass_Q = [store_ModParam_GVF_VegClass_Q,mod_param_GVF_VegClass_Q(1:365)];
    store_Realistic_Q = [store_Realistic_Q,Realistic_Q(1:365)];
end

store_baseline_Q = store_baseline_Q.*cms_to_cfs_conversion_factor;
store_ModParam_Q = store_ModParam_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_Q = store_ModParam_GVF_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_VegClass_Q = store_ModParam_GVF_VegClass_Q.*cms_to_cfs_conversion_factor;
store_Realistic_Q = store_Realistic_Q.*cms_to_cfs_conversion_factor;

baseline_ann_Q = nanmean(store_baseline_Q,1);
ModParam_ann_Q = nanmean(store_ModParam_Q,1);
ModParam_GVF_ann_Q = nanmean(store_ModParam_GVF_Q,1);
ModParam_GVF_vegclass_ann_Q = nanmean(store_ModParam_GVF_VegClass_Q,1);
realistic_ann_Q = nanmean(store_Realistic_Q,1);


delta_baseline = (baseline_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_realistic = (realistic_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam = (ModParam_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF = (ModParam_GVF_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF_vegclass = (ModParam_GVF_vegclass_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF_vegclass2 = (ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q)./ModParam_GVF_ann_Q;
delta_albedo = (realistic_ann_Q - ModParam_GVF_vegclass_ann_Q)./ModParam_GVF_vegclass_ann_Q;

sprintf('EBNFF: baseline result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_baseline)*100,min(delta_baseline)*100,max(delta_baseline)*100)
sprintf('EBNFF: ModParam  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam )*100,min(delta_ModParam )*100,max(delta_ModParam )*100)
sprintf('EBNFF: ModParam_GVF  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF )*100,min(delta_ModParam_GVF )*100,max(delta_ModParam_GVF )*100)
sprintf('EBNFF: ModParam_GVF_vegclass  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF_vegclass )*100,min(delta_ModParam_GVF_vegclass )*100,max(delta_ModParam_GVF_vegclass )*100)
sprintf('EBNFF: realistic result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_realistic)*100,min(delta_realistic)*100,max(delta_realistic)*100)
sprintf('EBNFF: veg class result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF_vegclass2)*100,min(delta_ModParam_GVF_vegclass2)*100,max(delta_ModParam_GVF_vegclass2)*100)
sprintf('EBNFF: albedo result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_albedo)*100,min(delta_albedo)*100,max(delta_albedo)*100)

%daily mod-param - baseline deltas
mean_baseline_Q = nanmean(store_baseline_Q(:));
delta_ModParam_daily_Q = (store_ModParam_Q - store_baseline_Q);
mean_rel_diff_ModParam_daily_Q = nanmean(delta_ModParam_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('EBNFF: DAILY mean diff from mod params is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_daily_Q*100,nanmean(delta_ModParam_daily_Q(:)), max(delta_ModParam_daily_Q(:)), length(idx))

%daily ModParam_GVF_ann_Q - baseline deltas
delta_ModParam_GVF_daily_Q = (store_ModParam_GVF_Q - store_baseline_Q);
mean_rel_diff_ModParam_GVF_daily_Q = nanmean(delta_ModParam_GVF_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_GVF_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('EBNFF: DAILY mean diff from modparams + GVF is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_GVF_daily_Q*100,nanmean(delta_ModParam_GVF_daily_Q(:)), max(delta_ModParam_GVF_daily_Q(:)), length(idx))

%daily ModParam_GVF_vegclass_ann_Q - baseline deltas
delta_ModParam_GVF_vegclass_daily_Q = (store_ModParam_GVF_VegClass_Q - store_baseline_Q);
mean_rel_diff_ModParam_GVF_vegclass_daily_Q = nanmean(delta_ModParam_GVF_vegclass_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_GVF_vegclass_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('EBNFF: DAILY mean diff from modparams + GVF +veg class is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_GVF_vegclass_daily_Q*100,nanmean(delta_ModParam_GVF_vegclass_daily_Q(:)), max(delta_ModParam_GVF_vegclass_daily_Q(:)), length(idx))

%compare winter ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q
idx_winter = 93:182;
store_ModParam_GVF_Q_winter = store_ModParam_GVF_Q(idx_winter,:);
store_ModParam_GVF_VegClass_Q_winter = store_ModParam_GVF_VegClass_Q(idx_winter,:);
pct_change_winter =  (store_ModParam_GVF_VegClass_Q_winter - store_ModParam_GVF_Q_winter)./store_ModParam_GVF_Q_winter;
mean_pct_change_winter = nanmean(pct_change_winter(:));
sprintf('EBNFF: DAILY mean diff from modparams + GVF +veg class vs GVF in winter is %.4f%%',mean_pct_change_winter*100)

%compare summer ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q
idx_summer = 213:273;
store_ModParam_GVF_Q_summer = store_ModParam_GVF_Q(idx_summer,:);
store_ModParam_GVF_VegClass_Q_summer = store_ModParam_GVF_VegClass_Q(idx_summer,:);
pct_change_summer =  (store_ModParam_GVF_VegClass_Q_summer - store_ModParam_GVF_Q_summer)./store_ModParam_GVF_Q_summer;
mean_pct_change_summer = nanmean(pct_change_summer(:));
sprintf('EBNFF: DAILY mean diff from modparams + GVF +veg class vs GVF in summer is %.4f%%',mean_pct_change_summer*100)

%compare winter realistic_ann_Q_winter - store_ModParam_GVF_VegClass_Q_winter
idx_winter = 93:182;
realistic_ann_Q_winter = store_Realistic_Q(idx_winter,:);
store_ModParam_GVF_VegClass_Q_winter = store_ModParam_GVF_VegClass_Q(idx_winter,:);
pct_change_winter =  (realistic_ann_Q_winter - store_ModParam_GVF_VegClass_Q_winter)./store_ModParam_GVF_VegClass_Q_winter;
mean_pct_change_winter = nanmean(pct_change_winter(:));
sprintf('EBNFF: DAILY mean diff from modparams + GVF +veg class +albedo vs GVF+class  in winter is %.4f%%',mean_pct_change_winter*100)

%compare summer realistic_ann_Q_summer - store_ModParam_GVF_VegClass_Q_summer
idx_summer = 183:273;
realistic_ann_Q_summer = store_Realistic_Q(idx_summer,:);
store_ModParam_GVF_VegClass_Q_summer = store_ModParam_GVF_VegClass_Q(idx_summer,:);
pct_change_summer =  (realistic_ann_Q_summer - store_ModParam_GVF_VegClass_Q_summer)./store_ModParam_GVF_VegClass_Q_summer;
mean_pct_change_summer = nanmean(pct_change_summer(:));
sprintf('EBNFF: DAILY mean diff from modparams + GVF +veg class +albedo vs GVF+class  in summer is %.4f%%',mean_pct_change_summer*100)

%compare pt1 realistic_ann_Q_pt1 - baseline_pt1
idx_pt1 = 1:212;
realistic_ann_Q_pt1 = store_Realistic_Q(idx_pt1,:);
store_baseline_Q_pt1 = store_baseline_Q(idx_pt1,:);
pct_change_pt1 =  (realistic_ann_Q_pt1 - store_baseline_Q_pt1)./store_baseline_Q_pt1;
mean_pct_change_pt1 = nanmean(pct_change_pt1(:));
sprintf('EBNFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in pt1 is %.4f%%',mean_pct_change_pt1*100)

%compare pt2 realistic_ann_Q_pt2 - baseline_pt2
idx_pt2 = 244:365;
realistic_ann_Q_pt2 = store_Realistic_Q(idx_pt2,:);
store_baseline_Q_pt2 = store_baseline_Q(idx_pt2,:);
pct_change_pt2 =  (realistic_ann_Q_pt2 - store_baseline_Q_pt2)./store_baseline_Q_pt2;
mean_pct_change_pt2 = nanmean(pct_change_pt2(:));
sprintf('EBNFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in pt2 is %.4f%%',mean_pct_change_pt2*100)

%compare may realistic_ann_Q_pt2 - baseline_pt2
idx_pt2 = 213:243;
realistic_ann_Q_pt2 = store_Realistic_Q(idx_pt2,:);
store_baseline_Q_pt2 = store_baseline_Q(idx_pt2,:);
pct_change_pt2 =  (realistic_ann_Q_pt2 - store_baseline_Q_pt2)./store_baseline_Q_pt2;
mean_pct_change_pt2 = nanmean(pct_change_pt2(:));
sprintf('EBNFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in May is %.4f%%',mean_pct_change_pt2*100)


baseline_hydrograph = nanmean(store_baseline_Q,2)';
mod_param_hydrograph = nanmean(store_ModParam_Q,2)';
mod_param_GVF_hydrograph = nanmean(store_ModParam_GVF_Q,2)';
mod_param_GVF_VegClass_hydrograph = nanmean(store_ModParam_GVF_VegClass_Q,2)';
Realistic_hydrograph = nanmean(store_Realistic_Q,2)';

%compute corresponding STDs:
baseline_hydrograph_std = std(store_baseline_Q','omitnan');
mod_param_hydrograph_std = std(store_ModParam_Q','omitnan');
store_ModParam_GVF_Q_std = std(store_ModParam_GVF_Q','omitnan');
store_ModParam_GVF_VegClass_Q_std = std(store_ModParam_GVF_VegClass_Q','omitnan');
Realistic_hydrograph_std = std(store_Realistic_Q','omitnan');

f=figure;
f.Position = [ -1535         -48        1139        1096];
subplot(2,1,1)
%title('East Branch North Fork Feather (station ICR)','fontsize',42)
hold on

%Plot basecase 
tpday = 365;
p01=plot(1:tpday,baseline_hydrograph(1:tpday),'-k','linewidth',3);
%Plot mod param case 
p02=plot(1:tpday,mod_param_hydrograph(1:tpday),'-b','linewidth',1.2);
%Plot mod param+gvf
p03=plot(1:tpday,mod_param_GVF_hydrograph(1:tpday),'-','linewidth',1.2,'color',[0 0.5 0]);
%Plot mod param+gvf+veg class
p04=plot(1:tpday,mod_param_GVF_VegClass_hydrograph(1:tpday),'--','linewidth',1.2,'color',[0 0.5 0]);
%realistic
p05=plot(1:tpday,Realistic_hydrograph(1:tpday),'-','linewidth',1.2,'color','r');
% % leg = legend([p01 p02 p03 p04 p05],{'Baseline','Modified params','Modified params+F_{veg}','Modified params+F_{veg}+Veg-class','Realistic (Modified params+F_{veg}+Veg-class+snow albedo)'},'fontsize',19,'location','best');
% % leg.Position = [0.6082    0.77    0.2911    0.15];
grid on
box on
ylabel({'multi-year mean daily'; 'Q (cfs) (2000-2022)'},'fontsize',42)
set(gca,'fontsize',42)
xlabel('day of water year','fontsize',42)
xlim([1 365])

%delta baseline plot:
subplot(2,1,2)
hold on
%Plot mod param case 
p02=plot(1:tpday,mod_param_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-b','linewidth',2);
%Plot mod param+gvf
p03=plot(1:tpday,mod_param_GVF_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-','linewidth',1.2,'color',[0 0.5 0]);
%Plot mod param+gvf+veg class
p04=plot(1:tpday,mod_param_GVF_VegClass_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'--','linewidth',1.2,'color',[0 0.5 0]);
%realistic
p05=plot(1:tpday,Realistic_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-','linewidth',1.2,'color','r');
grid on
box on
ylabel({'\Delta multi-year daily Q'; '(cfs) (2000-2022)'},'fontsize',42)
set(gca,'fontsize',42)
xlabel('day of water year','fontsize',42)
xlim([1 365])
saveas(f,'/Users/abolafia/ASO_Fire/Plots/Q_enmsemble_mean_2000_2022_EB_NF_Feather.png')

%include bar plot comparing annual means:
mean_Q = [nanmean(baseline_hydrograph) , nanmean(mod_param_hydrograph) , nanmean(mod_param_GVF_hydrograph),nanmean(mod_param_GVF_VegClass_hydrograph) ,nanmean(Realistic_hydrograph)];
X=categorical({'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'});
X = reordercats(X,{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'});
f=figure;
b=bar(X,mean_Q,'FaceColor','k');
grid on
box on
set(gca,'fontsize',20)
f.Position = [-1536         232         643         515];
ylabel('mean water year Q(cfs)','fontsize',28)
saveas(f,'/Users/abolafia/ASO_Fire/Plots/Q_ensemble_ANN_mean_barplot_EB_NF_Feather.png')

%% plot hydrographs during these years for the 3 burned catchments: North Fork Feather
STA_IDs = {'DCS','F56','F57','GYB','ICR','JBR','LCB','MER','MFP','NFP','NYS','ORH','SPK','TM1','TM2','TM3','WFR','YPB','YRS'};
stations = {'NFP'}; 

store_sta_idx=[];
for s=1:length(stations)
    IDX_sta = strcmp(STA_IDs,stations{1});
    idx = find(IDX_sta==1);
    store_sta_idx = [store_sta_idx;idx];
end

%load in simulated Q - baseline:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Feather_Baseline/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
Baseline_Streamflow = Simulated_Q.Streamflow;
Baseline_Streamflow = Baseline_Streamflow(:,store_sta_idx);
Baseline_Streamflow = nanmean(Baseline_Streamflow,2);
sim_dates = Simulated_Q.dates;
sim_datevecs = datevec(sim_dates);

%load in simulated Q - modified parameters
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParams_Streamflow = Simulated_Q.Streamflow;
ModParams_Streamflow = ModParams_Streamflow(:,store_sta_idx);
ModParams_Streamflow = nanmean(ModParams_Streamflow,2);
%load in simulated Q - modified parameters & veg class - BARE:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam_GVF/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParam_GVF = Simulated_Q.Streamflow;
ModParam_GVF = ModParam_GVF(:,store_sta_idx);
ModParam_GVF_streamflow = nanmean(ModParam_GVF,2);
%load in simulated Q - modified parameters & veg class - GRASS:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam_GVF_VegClass/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
ModParam_GVF_VegClass = Simulated_Q.Streamflow;
ModParam_GVF_VegClass = ModParam_GVF_VegClass(:,store_sta_idx);
ModParam_GVF_VegClass_streamflow = nanmean(ModParam_GVF_VegClass,2);

%load in simulated Q - modified parameters & veg class & albedo - BARE:
Simulated_Q = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Realistic/CADWR_Q/Streamflow_for_CADWR_Stations.mat');
Simulated_Q = Simulated_Q.Simulated_Streamflow;
Realistic = Simulated_Q.Streamflow;
Realistic = Realistic(:,store_sta_idx);
Realistic_Streamflow = nanmean(Realistic,2);

%compute daily mean hydrographs 
nyears = length(WYs);
store_baseline_Q=[];
store_ModParam_Q=[];
store_ModParam_GVF_Q=[];
store_ModParam_GVF_VegClass_Q=[];
store_Realistic_Q=[];

for y = 1:nyears
    current_year =  WYs(y);
    current_date_idx = find(sim_dates>=datenum([current_year-1 10 1]) & sim_dates<=datenum([current_year 9 30]));
    current_datevecs = sim_datevecs(current_date_idx,:);
    
    baseline_Q =Baseline_Streamflow(current_date_idx);
    mod_param_Q =ModParams_Streamflow(current_date_idx);
    mod_param_GVF_Q =ModParam_GVF_streamflow(current_date_idx);
    mod_param_GVF_VegClass_Q =ModParam_GVF_VegClass_streamflow(current_date_idx);
    Realistic_Q = Realistic_Streamflow(current_date_idx);
    
    %compute daily mean:
    [u,~,j] = unique(current_datevecs(:,1:3),'rows','stable');
    
    baseline_Q = accumarray(j,baseline_Q,[],@nanmean);
    mod_param_Q = accumarray(j,mod_param_Q,[],@nanmean);
    mod_param_GVF_Q = accumarray(j,mod_param_GVF_Q,[],@nanmean);
    mod_param_GVF_VegClass_Q = accumarray(j,mod_param_GVF_VegClass_Q,[],@nanmean);
    Realistic_Q = accumarray(j,Realistic_Q,[],@nanmean);
    
    store_baseline_Q = [store_baseline_Q,baseline_Q(1:365)];
    store_ModParam_Q = [store_ModParam_Q,mod_param_Q(1:365)];
    store_ModParam_GVF_Q = [store_ModParam_GVF_Q,mod_param_GVF_Q(1:365)];
    store_ModParam_GVF_VegClass_Q = [store_ModParam_GVF_VegClass_Q,mod_param_GVF_VegClass_Q(1:365)];
    store_Realistic_Q = [store_Realistic_Q,Realistic_Q(1:365)];
end

store_baseline_Q = store_baseline_Q.*cms_to_cfs_conversion_factor;
store_ModParam_Q = store_ModParam_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_Q = store_ModParam_GVF_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_VegClass_Q = store_ModParam_GVF_VegClass_Q.*cms_to_cfs_conversion_factor;
store_Realistic_Q = store_Realistic_Q.*cms_to_cfs_conversion_factor;

baseline_ann_Q = nanmean(store_baseline_Q,1);
ModParam_ann_Q = nanmean(store_ModParam_Q,1);
ModParam_GVF_ann_Q = nanmean(store_ModParam_GVF_Q,1);
ModParam_GVF_vegclass_ann_Q = nanmean(store_ModParam_GVF_VegClass_Q,1);
realistic_ann_Q = nanmean(store_Realistic_Q,1);


delta_baseline = (baseline_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_realistic = (realistic_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam = (ModParam_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF = (ModParam_GVF_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF_vegclass = (ModParam_GVF_vegclass_ann_Q - baseline_ann_Q)./baseline_ann_Q;
delta_ModParam_GVF_vegclass2 = (ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q)./ModParam_GVF_ann_Q;
delta_albedo = (realistic_ann_Q - ModParam_GVF_vegclass_ann_Q)./ModParam_GVF_vegclass_ann_Q;

sprintf('NFF: baseline result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_baseline)*100,min(delta_baseline)*100,max(delta_baseline)*100)
sprintf('NFF: ModParam  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam )*100,min(delta_ModParam )*100,max(delta_ModParam )*100)
sprintf('NFF: ModParam_GVF  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF )*100,min(delta_ModParam_GVF )*100,max(delta_ModParam_GVF )*100)
sprintf('NFF: ModParam_GVF_vegclass  result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF_vegclass )*100,min(delta_ModParam_GVF_vegclass )*100,max(delta_ModParam_GVF_vegclass )*100)
sprintf('NFF: realistic result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_realistic)*100,min(delta_realistic)*100,max(delta_realistic)*100)
sprintf('NFF: vegclass result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_ModParam_GVF_vegclass2)*100,min(delta_ModParam_GVF_vegclass2)*100,max(delta_ModParam_GVF_vegclass2)*100)
sprintf('NFF: albedo result in a mean annual delta of %.4f%% ranging from %.4f%%-%.4f%%', mean(delta_albedo)*100,min(delta_albedo)*100,max(delta_albedo)*100)

%daily mod-param - baseline deltas
mean_baseline_Q = nanmean(store_baseline_Q(:));
delta_ModParam_daily_Q = (store_ModParam_Q - store_baseline_Q);
mean_rel_diff_ModParam_daily_Q = nanmean(delta_ModParam_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('NFF: DAILY mean diff from mod params is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_daily_Q*100,nanmean(delta_ModParam_daily_Q(:)), max(delta_ModParam_daily_Q(:)), length(idx))

%daily ModParam_GVF_ann_Q - baseline deltas
delta_ModParam_GVF_daily_Q = (store_ModParam_GVF_Q - store_baseline_Q);
mean_rel_diff_ModParam_GVF_daily_Q = nanmean(delta_ModParam_GVF_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_GVF_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('NFF: DAILY mean diff from modparams + GVF is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_GVF_daily_Q*100,nanmean(delta_ModParam_GVF_daily_Q(:)), max(delta_ModParam_GVF_daily_Q(:)), length(idx))

%daily ModParam_GVF_vegclass_ann_Q - baseline deltas
delta_ModParam_GVF_vegclass_daily_Q = (store_ModParam_GVF_VegClass_Q - store_baseline_Q);
mean_rel_diff_ModParam_GVF_vegclass_daily_Q = nanmean(delta_ModParam_GVF_vegclass_daily_Q(:))./mean_baseline_Q;
idx=find(delta_ModParam_GVF_vegclass_daily_Q(:) > (0.1*mean_baseline_Q));
sprintf('NFF: DAILY mean diff from modparams + GVF +veg class is %.4f%% (%.4f cfs) of baseline Q and max diff is %.4f cfs, and exceed 10%% of mean daily Q %d times',mean_rel_diff_ModParam_GVF_vegclass_daily_Q*100,nanmean(delta_ModParam_GVF_vegclass_daily_Q(:)), max(delta_ModParam_GVF_vegclass_daily_Q(:)), length(idx))

%compare winter ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q
idx_winter = 93:182;
store_ModParam_GVF_Q_winter = store_ModParam_GVF_Q(idx_winter,:);
store_ModParam_GVF_VegClass_Q_winter = store_ModParam_GVF_VegClass_Q(idx_winter,:);
pct_change_winter =  (store_ModParam_GVF_VegClass_Q_winter - store_ModParam_GVF_Q_winter)./store_ModParam_GVF_Q_winter;
mean_pct_change_winter = nanmean(pct_change_winter(:));
sprintf('NFF: DAILY mean diff from modparams + GVF +veg class vs GVF in winter is %.4f%%',mean_pct_change_winter*100)

%compare summer ModParam_GVF_vegclass_ann_Q - ModParam_GVF_ann_Q
idx_summer = 213:273;
store_ModParam_GVF_Q_summer = store_ModParam_GVF_Q(idx_summer,:);
store_ModParam_GVF_VegClass_Q_summer = store_ModParam_GVF_VegClass_Q(idx_summer,:);
pct_change_summer =  (store_ModParam_GVF_VegClass_Q_summer - store_ModParam_GVF_Q_summer)./store_ModParam_GVF_Q_summer;
mean_pct_change_summer = nanmean(pct_change_summer(:));
sprintf('NFF: DAILY mean diff from modparams + GVF +veg class vs GVF in summer is %.4f%%',mean_pct_change_summer*100)

%compare winter realistic_ann_Q_winter - store_ModParam_GVF_VegClass_Q_winter
idx_winter = 93:182;
realistic_ann_Q_winter = store_Realistic_Q(idx_winter,:);
store_ModParam_GVF_VegClass_Q_winter = store_ModParam_GVF_VegClass_Q(idx_winter,:);
pct_change_winter =  (realistic_ann_Q_winter - store_ModParam_GVF_VegClass_Q_winter)./store_ModParam_GVF_VegClass_Q_winter;
mean_pct_change_winter = nanmean(pct_change_winter(:));
sprintf('NFF: DAILY mean diff from modparams + GVF +veg class +albedo vs GVF+class  in winter is %.4f%%',mean_pct_change_winter*100)

%compare summer realistic_ann_Q_summer - store_ModParam_GVF_VegClass_Q_summer
idx_summer = 183:273;
realistic_ann_Q_summer = store_Realistic_Q(idx_summer,:);
store_ModParam_GVF_VegClass_Q_summer = store_ModParam_GVF_VegClass_Q(idx_summer,:);
pct_change_summer =  (realistic_ann_Q_summer - store_ModParam_GVF_VegClass_Q_summer)./store_ModParam_GVF_VegClass_Q_summer;
mean_pct_change_summer = nanmean(pct_change_summer(:));
sprintf('NFF: DAILY mean diff from modparams + GVF +veg class +albedo vs GVF+class  in summer is %.4f%%',mean_pct_change_summer*100)

%compare pt1 realistic_ann_Q_pt1 - baseline_pt1
idx_pt1 = 1:212;
realistic_ann_Q_pt1 = store_Realistic_Q(idx_pt1,:);
store_baseline_Q_pt1 = store_baseline_Q(idx_pt1,:);
pct_change_pt1 =  (realistic_ann_Q_pt1 - store_baseline_Q_pt1)./store_baseline_Q_pt1;
mean_pct_change_pt1 = nanmean(pct_change_pt1(:));
sprintf('NFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in pt1 is %.4f%%',mean_pct_change_pt1*100)

%compare pt2 realistic_ann_Q_pt2 - baseline_pt2
idx_pt2 = 244:365;
realistic_ann_Q_pt2 = store_Realistic_Q(idx_pt2,:);
store_baseline_Q_pt2 = store_baseline_Q(idx_pt2,:);
pct_change_pt2 =  (realistic_ann_Q_pt2 - store_baseline_Q_pt2)./store_baseline_Q_pt2;
mean_pct_change_pt2 = nanmean(pct_change_pt2(:));
sprintf('NFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in pt2 is %.4f%%',mean_pct_change_pt2*100)

%compare may realistic_ann_Q_pt2 - baseline_pt2
idx_pt2 = 213:243;
realistic_ann_Q_pt2 = store_Realistic_Q(idx_pt2,:);
store_baseline_Q_pt2 = store_baseline_Q(idx_pt2,:);
pct_change_pt2 =  (realistic_ann_Q_pt2 - store_baseline_Q_pt2)./store_baseline_Q_pt2;
mean_pct_change_pt2 = nanmean(pct_change_pt2(:));
sprintf('NFF: DAILY mean diff from modparams + GVF +veg class +albedo vs baseline in May is %.4f%%',mean_pct_change_pt2*100)

baseline_hydrograph = nanmean(store_baseline_Q,2)';
mod_param_hydrograph = nanmean(store_ModParam_Q,2)';
mod_param_GVF_hydrograph = nanmean(store_ModParam_GVF_Q,2)';
mod_param_GVF_VegClass_hydrograph = nanmean(store_ModParam_GVF_VegClass_Q,2)';
Realistic_hydrograph = nanmean(store_Realistic_Q,2)';

%compute corresponding STDs:
baseline_hydrograph_std = std(store_baseline_Q','omitnan');
mod_param_hydrograph_std = std(store_ModParam_Q','omitnan');
store_ModParam_GVF_Q_std = std(store_ModParam_GVF_Q','omitnan');
store_ModParam_GVF_VegClass_Q_std = std(store_ModParam_GVF_VegClass_Q','omitnan');
Realistic_hydrograph_std = std(store_Realistic_Q','omitnan');

f=figure;
f.Position = [ -1535         -48        1139        1096];
subplot(2,1,1)
%title('North Fork Feather (station NFP)','fontsize',42)
hold on

%Plot basecase 
tpday = 365;
p01=plot(1:tpday,baseline_hydrograph(1:tpday),'-k','linewidth',3);
%Plot mod param case 
p02=plot(1:tpday,mod_param_hydrograph(1:tpday),'-b','linewidth',1.2);
%Plot mod param+gvf
p03=plot(1:tpday,mod_param_GVF_hydrograph(1:tpday),'-','linewidth',1.2,'color',[0 0.5 0]);
%Plot mod param+gvf+veg class
p04=plot(1:tpday,mod_param_GVF_VegClass_hydrograph(1:tpday),'--','linewidth',1.2,'color',[0 0.5 0]);
%realistic
p05=plot(1:tpday,Realistic_hydrograph(1:tpday),'-','linewidth',1.2,'color','r');
% % leg = legend([p01 p02 p03 p04 p05],{'Baseline','Modified params','Modified params+F_{veg}','Modified params+F_{veg}+Veg-class','Realistic (Modified params+F_{veg}+Veg-class+snow albedo)'},'fontsize',19,'location','best');
% % leg.Position = [0.6082    0.77    0.2911    0.15];
grid on
box on
ylabel({'multi-year mean daily'; 'Q (cfs) (2000-2022)'},'fontsize',42)
set(gca,'fontsize',42)
xlabel('day of water year','fontsize',42)
xlim([1 365])

%delta baseline plot:
subplot(2,1,2)
hold on
%Plot mod param case 
p02=plot(1:tpday,mod_param_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-b','linewidth',2);
%Plot mod param+gvf
p03=plot(1:tpday,mod_param_GVF_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-','linewidth',1.2,'color',[0 0.5 0]);
%Plot mod param+gvf+veg class
p04=plot(1:tpday,mod_param_GVF_VegClass_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'--','linewidth',1.2,'color',[0 0.5 0]);
%realistic
p05=plot(1:tpday,Realistic_hydrograph(1:tpday)-baseline_hydrograph(1:tpday),'-','linewidth',1.2,'color','r');
grid on
box on
ylabel({'\Delta multi-year daily Q'; '(cfs) (2000-2022)'},'fontsize',42)
set(gca,'fontsize',42)
xlabel('day of water year','fontsize',42)
xlim([1 365])
saveas(f,'/Users/abolafia/ASO_Fire/Plots/Q_enmsemble_mean_2000_2022_NFF.png')

%include bar plot comparing annual means:
mean_Q = [nanmean(baseline_hydrograph) , nanmean(mod_param_hydrograph) , nanmean(mod_param_GVF_hydrograph),nanmean(mod_param_GVF_VegClass_hydrograph) ,nanmean(Realistic_hydrograph)];
X=categorical({'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'});
X = reordercats(X,{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'});
f=figure;
b=bar(X,mean_Q,'FaceColor','k');
grid on
box on
set(gca,'fontsize',20)
f.Position = [-1536         232         643         515];
ylabel('mean water year Q(cfs)','fontsize',28)
saveas(f,'/Users/abolafia/ASO_Fire/Plots/Q_ensemble_ANN_mean_barplot_NFF.png')
