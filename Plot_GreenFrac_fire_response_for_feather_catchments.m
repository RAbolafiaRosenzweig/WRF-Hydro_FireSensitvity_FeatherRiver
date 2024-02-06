clc;clear all;close all;

%% load in LAI grid:
BA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200101.nc','XLAT_M');
BA_latvec = BA_lat(:);

BA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200101.nc','XLONG_M');
BA_lonvec = BA_lon(:);

%% load in fire perimeters:
Study_fires = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Fires_of_interest/firep_2018_2021.shp');
shape_info_fires = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Fires_of_interest/firep_2018_2021.shp');
p1_fire = shape_info_fires.CoordinateReferenceSystem;
nfire = length(Study_fires);

%% trim LAI grid to the fire perimeter bounding box:
lats=[];
lons=[];
for fi= 1:length(Study_fires)
    current_fire = Study_fires(fi);
    Shape_X=current_fire.X;
    Shape_Y = current_fire.Y;
    [fire_lat,fire_lon] = projinv(p1_fire,Shape_X,Shape_Y);
    idx_nan=find(isnan(fire_lat) | isnan(fire_lon));
    nshapes = length(idx_nan);
    lats = [lats;fire_lat'];
    lons = [lons;fire_lon'];
end

latlim = [min(lats) max(lats)];
lonlim = [min(lons) max(lons)];

IDX_BB = find(BA_latvec>=latlim(1) & BA_latvec<=latlim(2) & BA_lonvec>=lonlim(1) & BA_lonvec<=lonlim(2));

BA_latvec_BB = BA_latvec(IDX_BB);
BA_lonvec_BB = BA_lonvec(IDX_BB);

%% determine which points on lai grid are burned:
% % burn_years=[];
% % for fi= 1:length(Study_fires)
% %     current_fire = Study_fires(fi);
% %     Shape_X=current_fire.X;
% %     Shape_Y = current_fire.Y;
% %     [fire_lat,fire_lon] = projinv(p1_fire,Shape_X,Shape_Y);
% %     idx_nan=find(isnan(fire_lat) | isnan(fire_lon));
% %     nshapes = length(idx_nan);
% %
% %     start_iter = 1;
% %     store_burned_idx=[];
% %     for S = 1:nshapes
% %         current_fire_lat = fire_lat(start_iter:idx_nan(S)-1);
% %         current_fire_lon = fire_lon(start_iter:idx_nan(S)-1);
% %         [in] = inpolygon(BA_latvec_BB,BA_lonvec_BB,current_fire_lat,current_fire_lon);
% %         IDX_in = find(in==1);
% %         store_burned_idx = [store_burned_idx;IDX_in];
% %         start_iter = idx_nan(S)+1;
% %     end
% %     store_burned_idx = unique(store_burned_idx);
% %     store_burned_idx = sort(store_burned_idx);
% %     Burn_IDX_by_fire{fi} = store_burned_idx;
% %     burn_years = [burn_years;str2num(current_fire.YEAR_)];
% % end
% % Burned_LAI_IDX_byFire.Burn_IDX_by_fire = Burn_IDX_by_fire;
% % Burned_LAI_IDX_byFire.burn_years = burn_years;
% % save('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/LAI_info/Burned_LAI_IDX_byFire.mat','Burned_LAI_IDX_byFire', '-v7.3');
Burned_LAI_IDX_byFire = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/LAI_info/Burned_LAI_IDX_byFire.mat');
Burned_LAI_IDX_byFire=Burned_LAI_IDX_byFire.Burned_LAI_IDX_byFire;
Burn_IDX_by_fire = Burned_LAI_IDX_byFire.Burn_IDX_by_fire;
burn_years = Burned_LAI_IDX_byFire.burn_years;


%% collect LAI time series for pixels within each fire event:
start_day = datenum([2000 1 1]);
end_day = datenum([2022 9 30]);
datelist = start_day:end_day;
ndates = length(datelist);
datevecs = datevec(datelist);
MOD_GreenFrac_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/MOD15A2_LAI/Aggregated_8Day_Grids_GREENFRAC/';

% % store_dates=[];
% % store_GreenFrac_timeseries{1}= [];store_GreenFrac_timeseries{2}= [];store_GreenFrac_timeseries{3}= [];store_GreenFrac_timeseries{4}= [];store_GreenFrac_timeseries{5}= [];
% % for d=1:ndates
% %     year = datevecs(d,1);
% %     month = datevecs(d,2);
% %     day = datevecs(d,3);
% %     [year month day]
% %     infilename = sprintf('Gridded_MOD15A2_%04d%02d%02d.mat',year,month,day);
% %     if exist([MOD_GreenFrac_dir,infilename],'file') > 0
% %         store_dates = [store_dates;year,month,day];
% %
% %         GreenFrac_data = load([MOD_GreenFrac_dir,infilename]);
% %         GreenFrac_data=GreenFrac_data.MOD15A2_Fpar;
% %         GreenFrac_data = GreenFrac_data(:);
% %         GreenFrac_data_BB = GreenFrac_data(IDX_BB);
% %         for fi=1:nfire
% %            current_fi_idx =  Burn_IDX_by_fire{fi};
% %            current_GreenFrac = GreenFrac_data_BB(current_fi_idx)./100;
% %            current_mean_GreenFrac = nanmean(current_GreenFrac);
% %            store_GreenFrac_timeseries{fi} = [store_GreenFrac_timeseries{fi};current_mean_GreenFrac];
% %         end
% %     end
% % end
% % GreenFrac_timeseries_by_Fire.Fpar = store_GreenFrac_timeseries;
% % GreenFrac_timeseries_by_Fire.dates = store_dates;
% % save('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/LAI_info/GreenFrac_timeseries_by_Fire.mat','GreenFrac_timeseries_by_Fire', '-v7.3');
GreenFrac_timeseries_by_Fire = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/LAI_info/GreenFrac_timeseries_by_Fire.mat');
GreenFrac_timeseries_by_Fire = GreenFrac_timeseries_by_Fire.GreenFrac_timeseries_by_Fire;
GreenFrac_dates = GreenFrac_timeseries_by_Fire.dates;
GreenFrac_timeseries = GreenFrac_timeseries_by_Fire.Fpar;

%% Plot annual GreenFrac time series grouped by fire:
f=figure;
f.Position=[-1919         389        1920         408];
store_GreenFrac_reductions=[];
for fi=1:nfire
    fire_name = Study_fires(fi).FIRE_NAME;
    current_GreenFrac = GreenFrac_timeseries{fi};
    fire_year = burn_years(fi);
    
    %aggregate to annual time series:
    [u,~,j] = unique(GreenFrac_dates(:,1),'rows','stable');
    GreenFrac_Ann = accumarray(j,current_GreenFrac,[],@max);
    
    idx_prefire = find(u<fire_year);
    min_GreenFrac_prefire = min(GreenFrac_Ann(idx_prefire));
    mean_GreenFrac_prefire = nanmean(GreenFrac_Ann(idx_prefire));
    
    %plot data:
    subplot(1,5,fi)
    hold on
    title(fire_name,'fontsize',24)
    b1=bar(u,GreenFrac_Ann,'facecolor',[0.5 0.5 0.5],'edgecolor',[0.0 0.5 0.0],'barwidth',0.25);
    ylabel('Max Annual Veg. Fraction','fontsize',28)
    set(gca,'fontsize',28)
    xlim([min(u)-0.5 max(u)+0.5])
    xticks([2001:3:2022])
    xtickangle(90)
    ylim([0 , max(GreenFrac_Ann)])
    set(gca,'fontsize',22)
    grid on; box on
    %plot ref lines for min and mean pre-fire GreenFrac:
    plot(u,linspace(min_GreenFrac_prefire,min_GreenFrac_prefire,length(u)),'--r','linewidth',2.5)
    plot(u,linspace(mean_GreenFrac_prefire,mean_GreenFrac_prefire,length(u)),'--k','linewidth',2.5)
    %plot refline for fire year:
    %shade the burn year:
    X_left = ones(10,1)*(fire_year-0.5);
    X_right = ones(10,1)*(fire_year+0.5);
    X = [X_left;X_right];
    Vertical=linspace(0,20,10);
    inBetween=[Vertical,fliplr(Vertical)];
    fill(X,inBetween,[1, 0, 0],'FaceAlpha',0.2,'edgecolor','r')
    
    %report mean LAI reduction:
    pre_fire_Greenfrac = mean_GreenFrac_prefire;
    post_fire_Greenfrac = nanmean(GreenFrac_Ann(idx_prefire(end)+2:end));
    GreenFrac_reduction = (post_fire_Greenfrac - pre_fire_Greenfrac)./pre_fire_Greenfrac;
    store_GreenFrac_reductions = [store_GreenFrac_reductions;GreenFrac_reduction];
end

outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/GreenFrac_annual_timeseries_by_fire.png');
saveas(f,outfilename)

%% report a grid of mean burn FPAR reductions:
% % start_day = datenum([2000 1 1]);
% % end_day = datenum([2022 9 30]);
% % datelist = start_day:end_day;
% % ndates = length(datelist);
% % datevecs = datevec(datelist);
% % MOD_GreenFrac_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/MOD15A2_LAI/Aggregated_8Day_Grids_GREENFRAC/';
% % store_GreenFrac_timeseries{1}= [];store_GreenFrac_timeseries{2}= [];store_GreenFrac_timeseries{3}= [];store_GreenFrac_timeseries{4}= [];store_GreenFrac_timeseries{5}= [];
% % store_dates=[];
% % iter=0;
% % for d=1:ndates
% %     year = datevecs(d,1);
% %     month = datevecs(d,2);
% %     day = datevecs(d,3);
% %     [year month day]
% %     infilename = sprintf('Gridded_MOD15A2_%04d%02d%02d.mat',year,month,day);
% %     if exist([MOD_GreenFrac_dir,infilename],'file') > 0
% %         store_dates = [store_dates;year,month,day];
% %         iter = iter+1;
% %         GreenFrac_data = load([MOD_GreenFrac_dir,infilename]);
% %         GreenFrac_data=GreenFrac_data.MOD15A2_Fpar;
% %         GreenFrac_data = GreenFrac_data(:);
% %         GreenFrac_data_BB = GreenFrac_data(IDX_BB);
% %         for fi=1:nfire
% %            current_fi_idx =  Burn_IDX_by_fire{fi};
% %            %store coordinates
% %            if iter==1
% %                 BA_latvec_fire{fi} = BA_latvec_BB(current_fi_idx);
% %                 BA_lonvec_fire{fi} = BA_lonvec_BB(current_fi_idx);
% %            end
% %            %store fpar
% %            current_GreenFrac = GreenFrac_data_BB(current_fi_idx)./100;
% %            store_GreenFrac_timeseries{fi} = [store_GreenFrac_timeseries{fi},current_GreenFrac];
% %         end
% %     end
% % end
% % GreenFrac_timeseries_by_Fire.Fpar = store_GreenFrac_timeseries;
% % GreenFrac_timeseries_by_Fire.dates = store_dates;
% % GreenFrac_timeseries_by_Fire.lats = BA_latvec_fire;
% % GreenFrac_timeseries_by_Fire.lons = BA_lonvec_fire;
% % save('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/LAI_info/Gridded_GreenFrac_timeseries_by_Fire.mat','GreenFrac_timeseries_by_Fire', '-v7.3');

%loop through fires and match to the simulation grid:
FPAR_timeseries_by_Fire = load('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/LAI_info/Gridded_GreenFrac_timeseries_by_Fire.mat');
fire_lats = FPAR_timeseries_by_Fire.GreenFrac_timeseries_by_Fire.lats;
fire_lons = FPAR_timeseries_by_Fire.GreenFrac_timeseries_by_Fire.lons;
FPAR_dates = FPAR_timeseries_by_Fire.GreenFrac_timeseries_by_Fire.dates;
FPAR = FPAR_timeseries_by_Fire.GreenFrac_timeseries_by_Fire.Fpar;

LAT = double(ncread('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/Fire_Adjusted/wrfinput.nc','XLAT'));
LONG = double(ncread('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/Fire_Adjusted/wrfinput.nc','XLONG'));
GreenFrac_reduction_grid = nan(size(LAT));
for fi= 1:length(Study_fires)
    fi
    %get data for this fire
    current_FPAR = FPAR{fi};
    current_fire_lats = fire_lats{fi};
    current_fire_lons = fire_lons{fi};
    current_dates = FPAR_dates;
    fire_year = burn_years(fi);
    
    %calculate FPAR pre- to post-fire change:
    idx_prefire = find(current_dates(:,1) < fire_year);
    idx_postfire = find(current_dates(:,1) > fire_year);
    FPAR_prefire = current_FPAR(:,idx_prefire);
    FPAR_postfire = current_FPAR(:,idx_postfire);
    
    pre_fire_dates = current_dates(idx_prefire,:);
    post_fire_dates = current_dates(idx_postfire,:);
    
    %find median(max ann Fveg):
    unique_prefire_years = unique(pre_fire_dates(:,1));
    store_max_prefire_FPAR=[];
    for yi=1:length(unique_prefire_years)
        current_year = unique_prefire_years(yi);
        idx_year = find(pre_fire_dates(:,1) == current_year);
        current_FPAR = FPAR_prefire(:,idx_year);
        max_FPAR = max(current_FPAR,[],2);
        store_max_prefire_FPAR = [store_max_prefire_FPAR,max_FPAR];
    end
    Pre_fire_FPAR_max = nanmedian(store_max_prefire_FPAR,2);
    
    unique_postfire_years = unique(post_fire_dates(:,1));
    store_max_postfire_FPAR=[];
    for yi=1:length(unique_postfire_years)
        current_year = unique_postfire_years(yi);
        idx_year = find(post_fire_dates(:,1) == current_year);
        current_FPAR = FPAR_postfire(:,idx_year);
        max_FPAR = max(current_FPAR,[],2);
        store_max_postfire_FPAR = [store_max_postfire_FPAR,max_FPAR];
    end
    post_fire_FPAR_max = max(store_max_postfire_FPAR,[],2);
    
    Percent_Change = (post_fire_FPAR_max - Pre_fire_FPAR_max)./Pre_fire_FPAR_max;
    %constrain to 0:
    idx = find(Percent_Change > 0);
    Percent_Change(idx) = 0;
    %spatiall match fire perimeter to AORC grid:
    fire_coords = [current_fire_lats';current_fire_lons'];
    geo_coords = [LAT(:)' ; LONG(:)'];
    IDX_NN=nearestneighbour(fire_coords,geo_coords);
    
    GreenFrac_reduction_grid(IDX_NN) = Percent_Change;
end


%write outputs:
outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/LAI_info/GreenFrac_reduction_grid.csv');
dlmwrite(outfilename,GreenFrac_reduction_grid,'delimiter',',','precision',8);

% Create a spatial plot of this grid:
Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
ncatch = length(Catchments);

%plot the shapes
f=figure;
sz=1;
geoscatter(LAT(1),LONG(1),sz,[0.9 0.9 0.9],'.');
hold on
%format
set(gca,'fontsize',25)
%show state lines:
states=shaperead('usastatehi', 'UseGeoCoords', true);
for i=1:49
    lat = states(i).Lat;
    lon = states(i).Lon;
    geoplot(lat,lon,'LineWidth',1.5,'color','k','linewidth',2);
end

latlim=[39.2, 40+40/60];
lonlim = [-123+45/60 -120];
geolimits(latlim,lonlim)

nlevels=14;
cmap1 = flipud(summer(nlevels/2));
cmap3 = copper(nlevels/2);
cmap=[cmap3;cmap1];
colormap(cmap)

idx0=find(GreenFrac_reduction_grid(:)>=0 | isnan(GreenFrac_reduction_grid(:)));
C=GreenFrac_reduction_grid;
C(idx0)=NaN;
sz=40;
geoscatter(LAT(:),LONG(:),sz,C(:),'.');
c=colorbar;
title('GVF reduction factor')

%show shape boundaries:
n = length(Catchments);
for i=1:3
    current_shape =  Catchments(i);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [lat,lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(lat==90 | isnan(lat) | isnan(lon));
    lat(idx)=[];
    lon(idx)=[];
    geoplot(lat,lon,'-k','linewidth',2)
end

outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/Spatial_GreenFrac_reduction.png');
saveas(f,outfilename)