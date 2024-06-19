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
store_dates=[];
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

    store_baseline_Q = [store_baseline_Q;baseline_Q];
    store_ModParam_Q = [store_ModParam_Q;mod_param_Q];
    store_ModParam_GVF_Q = [store_ModParam_GVF_Q;mod_param_GVF_Q];
    store_ModParam_GVF_VegClass_Q = [store_ModParam_GVF_VegClass_Q;mod_param_GVF_VegClass_Q];
    store_Realistic_Q = [store_Realistic_Q;Realistic_Q];

    store_dates = [store_dates;u];
end

store_baseline_Q = store_baseline_Q.*cms_to_cfs_conversion_factor;
store_ModParam_Q = store_ModParam_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_Q = store_ModParam_GVF_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_VegClass_Q = store_ModParam_GVF_VegClass_Q.*cms_to_cfs_conversion_factor;
store_Realistic_Q = store_Realistic_Q.*cms_to_cfs_conversion_factor;


geofilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/Fire_Adjusted/wrfinput_Feather_postifre_Grassland.nc';
LAT = ncread(geofilename,'XLAT');
LON = ncread(geofilename,'XLONG');
%load in study catchment shapes:
Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
Catchments = Catchments(2);
ncatch = length(Catchments);
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
    store_Catch_PRCP = [store_Catch_PRCP;mean_catch_PRCP];
end

store_datenums = datenum(store_dates);
for WY=2000:2022
    current_idx=find( store_datenums>=datenum([WY-1 10 1]) & store_datenums<=datenum([WY 9 30]));
    f=figure;
    f.Position = [-1919          71         720         987];
    subplot(2,1,1)
    title('Middle Fork (station MER)','fontsize',40)
    hold on
    plot_dates = datenum(store_datenums(current_idx));
    %Plot basecase
    tpday = 365;
    p01=plot(plot_dates,store_baseline_Q(current_idx),'-k','linewidth',2);
    %Plot mod param case
    p02=plot(plot_dates,store_ModParam_Q(current_idx),'-b','linewidth',1.2);
    %Plot mod param+gvf
    p03=plot(plot_dates,store_ModParam_GVF_Q(current_idx),'-','linewidth',1.2,'color',[0 0.5 0]);
    %Plot mod param+gvf+veg class
    p04=plot(plot_dates,store_ModParam_GVF_VegClass_Q(current_idx),'--','linewidth',1.2,'color',[0 0.5 0]);
    %realistic
    p05=plot(plot_dates,store_Realistic_Q(current_idx),'-','linewidth',1.2,'color','r');
%     if WY==2000
%         leg = legend([p01 p02 p03 p04 p05],{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'},'fontsize',28);
%         leg.Position = [ 0.1993    0.7495    0.7597    0.1733];
%     end
    grid on
    box on
    ylabel({'daily Q'; '(cfs) (2000-2022)'},'fontsize',30)
    set(gca,'fontsize',30)
    xlim([min(plot_dates) max(plot_dates)])
    datetick('x','mm/yyyy','keepticks','keeplimits')
    xtickangle(30)
    title(sprintf('WY: %04d',WY),'fontsize',30)
    ylim([0 max(store_Realistic_Q)])

    %include hyetograph:
    yyaxis right
    right_color='b';
    axR=gca;
    axR.YColor=right_color;
    prcp_plot = stem(plot_dates,store_Catch_PRCP(current_idx),'b');
    set(gca,'ydir','reverse');
    ylabel({'Prcp.'; '(mm/day)'},'fontsize',28)
    ylim([0 max(store_Catch_PRCP).*1.5])
    xlim([min(plot_dates) max(plot_dates)])


    %delta baseline plot:
    subplot(2,1,2)
    hold on
    %Plot mod param case
    p02=plot(plot_dates,store_ModParam_Q(current_idx)-store_baseline_Q(current_idx),'-b','linewidth',2);
    %Plot mod param+gvf
    p03=plot(plot_dates,store_ModParam_GVF_Q(current_idx)-store_baseline_Q(current_idx),'-','linewidth',1.2,'color',[0 0.5 0]);
    %Plot mod param+gvf+veg class
    p04=plot(plot_dates,store_ModParam_GVF_VegClass_Q(current_idx)-store_baseline_Q(current_idx),'--','linewidth',1.2,'color',[0 0.5 0]);
    %realistic
    p05=plot(plot_dates,store_Realistic_Q(current_idx)-store_baseline_Q(current_idx),'-','linewidth',1.2,'color','r');
    grid on
    box on
    ylabel('\Delta daily Q (cfs)','fontsize',30)
    set(gca,'fontsize',30)

    xlim([min(plot_dates) max(plot_dates)])
    datetick('x','mm/yyyy','keepticks','keeplimits')
    xtickangle(30)
    ylim([-1000 max(store_Realistic_Q-store_baseline_Q)])

    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/Q_total_timeseries_2000_2022_MiddleFork_WY%04d.png',WY))
end

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

    store_baseline_Q = [store_baseline_Q;baseline_Q];
    store_ModParam_Q = [store_ModParam_Q;mod_param_Q];
    store_ModParam_GVF_Q = [store_ModParam_GVF_Q;mod_param_GVF_Q];
    store_ModParam_GVF_VegClass_Q = [store_ModParam_GVF_VegClass_Q;mod_param_GVF_VegClass_Q];
    store_Realistic_Q = [store_Realistic_Q;Realistic_Q];

end

store_baseline_Q = store_baseline_Q.*cms_to_cfs_conversion_factor;
store_ModParam_Q = store_ModParam_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_Q = store_ModParam_GVF_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_VegClass_Q = store_ModParam_GVF_VegClass_Q.*cms_to_cfs_conversion_factor;
store_Realistic_Q = store_Realistic_Q.*cms_to_cfs_conversion_factor;

geofilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/Fire_Adjusted/wrfinput_Feather_postifre_Grassland.nc';
LAT = ncread(geofilename,'XLAT');
LON = ncread(geofilename,'XLONG');
%load in study catchment shapes:
Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
Catchments = Catchments(1);
ncatch = length(Catchments);
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
    store_Catch_PRCP = [store_Catch_PRCP;mean_catch_PRCP];
end

for WY=2000:2022
    current_idx=find( store_datenums>=datenum([WY-1 10 1]) & store_datenums<=datenum([WY 9 30]));
    f=figure;
    f.Position = [-1919          71         720         987];
    subplot(2,1,1)
    hold on
    plot_dates = datenum(store_datenums(current_idx));
    %Plot basecase
    tpday = 365;
    p01=plot(plot_dates,store_baseline_Q(current_idx),'-k','linewidth',2);
    %Plot mod param case
    p02=plot(plot_dates,store_ModParam_Q(current_idx),'-b','linewidth',1.2);
    %Plot mod param+gvf
    p03=plot(plot_dates,store_ModParam_GVF_Q(current_idx),'-','linewidth',1.2,'color',[0 0.5 0]);
    %Plot mod param+gvf+veg class
    p04=plot(plot_dates,store_ModParam_GVF_VegClass_Q(current_idx),'--','linewidth',1.2,'color',[0 0.5 0]);
    %realistic
    p05=plot(plot_dates,store_Realistic_Q(current_idx),'-','linewidth',1.2,'color','r');
%     if WY==2005
%         leg = legend([p01 p02 p03 p04 p05],{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'},'fontsize',28);
%         leg.Position = [ 0.1993    0.7495    0.7597    0.1733];
%     end
    grid on
    box on
    ylabel({'daily Q'; '(cfs) (2000-2022)'},'fontsize',30)
    set(gca,'fontsize',30)
    xlim([min(plot_dates) max(plot_dates)])
    datetick('x','mm/yyyy','keepticks','keeplimits')
    xtickangle(30)
    title(sprintf('WY: %04d',WY),'fontsize',30)
    ylim([0 max(store_Realistic_Q)])

    %include hyetograph:
    yyaxis right
    right_color='b';
    axR=gca;
    axR.YColor=right_color;
    prcp_plot = stem(plot_dates,store_Catch_PRCP(current_idx),'b');
    set(gca,'ydir','reverse');
    ylabel({'Prcp.'; '(mm/day)'},'fontsize',28)
    ylim([0 max(store_Catch_PRCP).*1.5])
    xlim([min(plot_dates) max(plot_dates)])


    %delta baseline plot:
    subplot(2,1,2)
    hold on
    %Plot mod param case
    p02=plot(plot_dates,store_ModParam_Q(current_idx)-store_baseline_Q(current_idx),'-b','linewidth',2);
    %Plot mod param+gvf
    p03=plot(plot_dates,store_ModParam_GVF_Q(current_idx)-store_baseline_Q(current_idx),'-','linewidth',1.2,'color',[0 0.5 0]);
    %Plot mod param+gvf+veg class
    p04=plot(plot_dates,store_ModParam_GVF_VegClass_Q(current_idx)-store_baseline_Q(current_idx),'--','linewidth',1.2,'color',[0 0.5 0]);
    %realistic
    p05=plot(plot_dates,store_Realistic_Q(current_idx)-store_baseline_Q(current_idx),'-','linewidth',1.2,'color','r');
    grid on
    box on
    ylabel('\Delta daily Q (cfs)','fontsize',30)
    set(gca,'fontsize',30)

    xlim([min(plot_dates) max(plot_dates)])
    datetick('x','mm/yyyy','keepticks','keeplimits')
    xtickangle(30)
    ylim([-1000 max(store_Realistic_Q-store_baseline_Q)])
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/Q_total_timeseries_2000_2022_EBNFF_WY%04d.png',WY))
end

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

    store_baseline_Q = [store_baseline_Q;baseline_Q];
    store_ModParam_Q = [store_ModParam_Q;mod_param_Q];
    store_ModParam_GVF_Q = [store_ModParam_GVF_Q;mod_param_GVF_Q];
    store_ModParam_GVF_VegClass_Q = [store_ModParam_GVF_VegClass_Q;mod_param_GVF_VegClass_Q];
    store_Realistic_Q = [store_Realistic_Q;Realistic_Q];

    store_dates = [store_dates;u];
end

store_baseline_Q = store_baseline_Q.*cms_to_cfs_conversion_factor;
store_ModParam_Q = store_ModParam_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_Q = store_ModParam_GVF_Q.*cms_to_cfs_conversion_factor;
store_ModParam_GVF_VegClass_Q = store_ModParam_GVF_VegClass_Q.*cms_to_cfs_conversion_factor;
store_Realistic_Q = store_Realistic_Q.*cms_to_cfs_conversion_factor;

geofilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/Fire_Adjusted/wrfinput_Feather_postifre_Grassland.nc';
LAT = ncread(geofilename,'XLAT');
LON = ncread(geofilename,'XLONG');
%load in study catchment shapes:
Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
Catchments = Catchments(3);
ncatch = length(Catchments);
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
    store_Catch_PRCP = [store_Catch_PRCP;mean_catch_PRCP];
end

for WY=2000:2022
    current_idx=find( store_datenums>=datenum([WY-1 10 1]) & store_datenums<=datenum([WY 9 30]));
    f=figure;
    f.Position = [-1919          71         720         987];
    subplot(2,1,1)
    hold on
    plot_dates = datenum(store_datenums(current_idx));
    %Plot basecase
    tpday = 365;
    p01=plot(plot_dates,store_baseline_Q(current_idx),'-k','linewidth',2);
    %Plot mod param case
    p02=plot(plot_dates,store_ModParam_Q(current_idx),'-b','linewidth',1.2);
    %Plot mod param+gvf
    p03=plot(plot_dates,store_ModParam_GVF_Q(current_idx),'-','linewidth',1.2,'color',[0 0.5 0]);
    %Plot mod param+gvf+veg class
    p04=plot(plot_dates,store_ModParam_GVF_VegClass_Q(current_idx),'--','linewidth',1.2,'color',[0 0.5 0]);
    %realistic
    p05=plot(plot_dates,store_Realistic_Q(current_idx),'-','linewidth',1.2,'color','r');
%     if WY==2005
%         leg = legend([p01 p02 p03 p04 p05],{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'},'fontsize',28);
%         leg.Position = [ 0.1993    0.7495    0.7597    0.1733];
%     end
    grid on
    box on
    ylabel({'daily Q'; '(cfs) (2000-2022)'},'fontsize',30)
    set(gca,'fontsize',30)
    xlim([min(plot_dates) max(plot_dates)])
    datetick('x','mm/yyyy','keepticks','keeplimits')
    xtickangle(30)
    title(sprintf('WY: %04d',WY),'fontsize',30)
    ylim([0 max(store_Realistic_Q)])

    %include hyetograph:
    yyaxis right
    right_color='b';
    axR=gca;
    axR.YColor=right_color;
    prcp_plot = stem(plot_dates,store_Catch_PRCP(current_idx),'b');
    set(gca,'ydir','reverse');
    ylabel({'Prcp.'; '(mm/day)'},'fontsize',28)
    ylim([0 max(store_Catch_PRCP).*1.5])
    xlim([min(plot_dates) max(plot_dates)])


    %delta baseline plot:
    subplot(2,1,2)
    hold on
    %Plot mod param case
    p02=plot(plot_dates,store_ModParam_Q(current_idx)-store_baseline_Q(current_idx),'-b','linewidth',2);
    %Plot mod param+gvf
    p03=plot(plot_dates,store_ModParam_GVF_Q(current_idx)-store_baseline_Q(current_idx),'-','linewidth',1.2,'color',[0 0.5 0]);
    %Plot mod param+gvf+veg class
    p04=plot(plot_dates,store_ModParam_GVF_VegClass_Q(current_idx)-store_baseline_Q(current_idx),'--','linewidth',1.2,'color',[0 0.5 0]);
    %realistic
    p05=plot(plot_dates,store_Realistic_Q(current_idx)-store_baseline_Q(current_idx),'-','linewidth',1.2,'color','r');
    grid on
    box on
    ylabel('\Delta daily Q (cfs)','fontsize',30)
    set(gca,'fontsize',30)

    xlim([min(plot_dates) max(plot_dates)])
    datetick('x','mm/yyyy','keepticks',['kee3v' ...
        '?Å¸ŒÍ≈
wåQ4
\3]w2qa%$^%H3gw 
    xtickangle(30)
    ylim([-1000 max(store_Realistic_Q-store_baseline_Q)])
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/Q_total_timeseries_2000_2022_NFF_WY%04d.png',WY))
end
