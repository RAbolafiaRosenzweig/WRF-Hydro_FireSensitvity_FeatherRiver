clc;clear all;close all;

%load in study shapes:
Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
ncatch = length(Catchments);

Study_fires = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Fires_of_interest/firep_2018_2021.shp');
shape_info_fires = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Fires_of_interest/firep_2018_2021.shp');
p1_fire = shape_info_fires.CoordinateReferenceSystem;
nfire = length(Study_fires);

%plot the shapes
f=figure;
f.Position = [305    96   821   589];

%overlay statelines
states=shaperead('usastatehi', 'UseGeoCoords', true);
for s=1:length(states)
    current_state = states(s);
    lats = current_state.Lat;
    lons = current_state.Lon;
    idx=find(isnan(lats) | isnan(lons));
    if length(idx) >0
        lats(idx)=[];
        lons(idx)=[];
    end
    geoplot(lats,lons,'-k','linewidth',4)
    if s==1
        hold on
    end

end
% p2=geoshow(states,'FaceColor','none','Linewidth',5);
%show shape boundaries:

geobasemap topographic
geolimits([39+25/60 , 40+50/60],[-122 -119.85])
% % latlim=[38+46/60 , 40+40/60];
% % lonlim = [-123+45/60 -119.75];
% % xlim(lonlim)
% % ylim(latlim)

COLORS = [0.6350, 0.0780, 0.1840;...
    0.8500, 0.3250, 0.0980;...
    1, 0, 0;...
    0.75, 0.75, 0;...
    0.25, 0.25, 0.25];

COLORS = hot(30);
% % for fi= 1:length(Study_fires)
% %     current_fire = Study_fires(fi);
% %     Shape_X=current_fire.X;
% %     Shape_Y = current_fire.Y;
% %     [fire_lat,fire_lon] = projinv(p1_fire,Shape_X,Shape_Y);
% %     plot(fire_lon,fire_lat,'-','color',COLORS(1,:),'linewidth',2)
% % end
iter=0;
for fi= 1:length(Study_fires)
    current_fire = Study_fires(fi);
    Shape_X=current_fire.X;
    Shape_Y = current_fire.Y;
    [fire_lat,fire_lon] = projinv(p1_fire,Shape_X,Shape_Y);
    idx_nan=find(isnan(fire_lat) | isnan(fire_lon));
    nshapes = length(idx_nan);

    start_iter = 1;
    for S = 1:nshapes
        current_fire_lat = fire_lat(start_iter:idx_nan(S)-1);
        current_fire_lon = fire_lon(start_iter:idx_nan(S)-1);
%         c=geoplot(current_fire_lat, current_fire_lon, 'r');
        fire_poly = geopolyshape(current_fire_lat,current_fire_lon);
        geoplot(fire_poly,'facecolor',COLORS(7,:),'facealpha',0.5)
 
% % %         fire_poly = polyshape(current_fire_lon,current_fire_lat);
% %         plot(fire_poly,'facecolor',COLORS(7,:),'edgecolor',COLORS(7,:),'FaceAlpha',0.1,'linewidth',3)
        start_iter = idx_nan(S)+1;
    end
    iter = iter;
end

set(gca,'fontsize',24')

n = length(Catchments);
for i=1:n-1
    current_shape =  Catchments(i);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [lat,lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(lat==90);
    lat(idx)=[];
    lon(idx)=[];
    geoplot(lat,lon,'-k','linewidth',6)
    %     plot(lon,lat,'-k','linewidth',4)
end

% % %plot gages:
% % %get simulated Q data:
% % Daily_sim_Q = readtable('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/Feather_Baseline/Baseline_Streamflow.Rdata_daily.csv');
% % Dates = Daily_sim_Q.UTC_date;
% % sim_IDs = Daily_sim_Q.site_no;
% % q_cms_sim = Daily_sim_Q.q_cms;
% %
% % Site_IDs = Daily_sim_Q.site_no;
% % unique_IDS = unique(Site_IDs);
% % n_IDs = length(unique_IDS);
% %
% % %only consider data in the catchment polygons:
% % Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
% % shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
% % p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
% % ncatch = length(Catchments);
% %
% % % get meta data for USGS gages:
% % USGS_metadata = readtable('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Observations/USGS/metadata/streamflow_usgs_metadata.csv');
% % IDs = USGS_metadata.site_no;
% % USGS_lats =USGS_metadata.dec_lat_va;
% % USGS_lons =USGS_metadata.dec_long_va;
% %
% %
% % store_data=[];
% % for s=1:n_IDs
% %     s
% %     current_ID = unique_IDS(s);
% %     idx = find(IDs==current_ID);
% %     lat = USGS_lats(idx);
% %     lon = USGS_lons(idx);
% %     %check if this coord is in catchment polygons:
% %     for c=1:ncatch
% %         current_shape =  Catchments(c);
% %         Shape_X=current_shape.X;
% %         Shape_Y = current_shape.Y;
% %         [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
% %         [in,on] = inpolygon(lat,lon,catch_lat,catch_lon);
% %
% %         if max([in;on]) > 0
% %             store_data = [store_data;lat,lon,current_ID];
% %         end
% %     end
% % end
% %
% % plot(store_data(:,2),store_data(:,1),'bo','markersize',6,'markerfacecolor','b')

% % %plot SNOTEL sites:
% % SNOTEL_metadata = readtable('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Observations/SNOTEL/obsSnoData_bcqc_MetaData.csv');
% % snotel_lat = SNOTEL_metadata.lat_bcqc;
% % snotel_lon = SNOTEL_metadata.lon_bcqc;
% %
% % plot(snotel_lon,snotel_lat,'bx','markersize',8)
% %
% %
% % GAGESII_metadata = shaperead('/Volumes/Pruina_External_Elements/Fire_ModelErrors/GAGESII/gagesII_9322_point_shapefile/gagesII_9322_sept30_2011.shp');
% % shape_info = shapeinfo('/Volumes/Pruina_External_Elements/Fire_ModelErrors/GAGESII/gagesII_9322_point_shapefile/gagesII_9322_sept30_2011.shp');
% % p1 = shape_info.CoordinateReferenceSystem;
% % ngages = length(GAGESII_metadata);
% % iter=0;
% % store_coord=[];
% % for s=1:ngages
% %     %define which ones are in the catchment
% %     [lat,lon] = projinv(p1,GAGESII_metadata(s).X,GAGESII_metadata(s).Y);
% %     for c=1:ncatch
% %         current_shape =  Catchments(c);
% %         Shape_X=current_shape.X;
% %         Shape_Y = current_shape.Y;
% %         [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
% %         [in,on] = inpolygon(lat,lon,catch_lat,catch_lon);
% %         if max([in;on]) > 0
% %             iter=iter+1;
% %             plot(lon,lat,'rx','markersize',8)
% %             store_ID{iter} = GAGESII_metadata(s).STAID;
% %             store_coord = [store_coord;lat,lon];
% %         end
% %     end
% % end

%plot CA station data as well:
STA_IDs = {'DCS','F56','F57','GYB','ICR','JBR','LCB','MER','MFP','NFP','NYS','ORH','SPK','TM1','TM2','TM3','WFR','YPB','YRS'};
CA_station_data = readtable('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Observations/USGS/CA_data/CDEC_Hourly-RealTime_stationInfo.csv');
CA_lats = CA_station_data.Latitude;
CA_lons = CA_station_data.Longitude;
ngages = length(STA_IDs);
store_lat_lon=[];
for s=1:ngages
    current_station = STA_IDs{s};
    if s==5 || s== 8 || s== 10 %only plot gages of interest
        IDX = strcmp(CA_station_data.ID,current_station);
        idx = find(IDX==1);
        lat = CA_lats(idx);
        lon = CA_lons(idx);
        geoscatter(lat,lon,100,'markerfacecolor','k','markeredgecolor','k')
% %         plot(lon,lat,'ko','markersize',20,'markerfacecolor','k')

        if s<10
            text(lon+0.08,lat,current_station,'fontsize',20)
        else
            text(lon-0.17,lat,current_station,'fontsize',20)
        end

        store_lat_lon = [store_lat_lon;idx,lat,lon];
    end
end

%add feather river line:
rivers = shaperead('/Users/abolafia/Downloads/rv16my07/rv16my07.shp');
for r=1:length(rivers)
    current_name = rivers(r).PNAME;
    current_x = rivers(r).X;
    current_y = rivers(r).Y;
    median_x =nanmedian(current_x);
    median_y = nanmedian(current_y);
    if median_x>=-121.6 && median_x<=-120.2 && median_y>=39.4 && median_y<=40.55
        geoplot(current_y,current_x,'linewidth',1,'color',[0.3010, 0.7450, 0.9330])
        if length(current_name) >=7
            IDX = current_name(1:7) == 'FEATHER';
            if min(IDX)==1
                geoplot(current_y,current_x,'linewidth',3,'color','b')
            end
        end
    end
end

%add lake Oroville:
CA_lakes = readgeotable('/Volumes/Pruina_External_Elements/Shapefiles/California_Lakes/California_Lakes.shp');
lake_id = 18498;
Shapes = CA_lakes.Shape;
Oroville = Shapes(lake_id);
a=geoplot(Oroville,'facecolor',[0, 0.4470, 0.7410],'facealpha',1);

%save
outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/StudyDomain.png');
saveas(f,outfilename)

bsi_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/';
fire_names = {'NorthComplex_2020_MTBS', 'Camp_2018_MTBS','Walker_2019_MTBS','Sugar_2021_BAER','Dixie_2021_BAER'};
suffix = '_1000m.nc';
titles = {'North Complex (Aug 2020)' , 'Camp (Nov 2018)','Walker (Sept 2019)','Sugar (Jul 2021)','Dixie (Jul 2021)'};
Study_fires = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Fires_of_interest/firep_2018_2021.shp');
shape_info_fires = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Fires_of_interest/firep_2018_2021.shp');
p1_fire = shape_info_fires.CoordinateReferenceSystem;
nfire = length(Study_fires);

COLORS = hot(30);
cmap = [0 , 0.5 , 0 ; ...
    1.0000    1.0000    0.5000;...
    1.0000    0.4286         0;...
    0.4286         0         0];

feather_geofilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/Fire_Adjusted/wrfinput.nc';
LAT = ncread(feather_geofilename,'XLAT');
LON = ncread(feather_geofilename,'XLONG');
for fi = 1:length(fire_names)
    current_filename = [bsi_dir,fire_names{fi},suffix];
    bsi = ncread(current_filename,'Band1');
    % %    lat = ncread(current_filename,'y');
    % %    lon = ncread(current_filename,'x');
    % %    [LAT,LON] = meshgrid(lat,lon);
    Z=nan(size(bsi));
    idx0 = find(bsi==0);
    idx1 = find(bsi==1);
    idx2 = find(bsi==2);
    idx3 = find(bsi==3);
    idx4 = find(bsi==4);
    %    Z(idx0) = 1;
    Z(idx1) = 1;
    Z(idx2) = 2;
    Z(idx3) = 3;
    Z(idx4) = 4;



    f=figure;
    hold on

    current_fire = Study_fires(fi);
    Shape_X=current_fire.X;
    Shape_Y = current_fire.Y;
    [fire_lat,fire_lon] = projinv(p1_fire,Shape_X,Shape_Y);
    plot(fire_lon,fire_lat,'-','color',COLORS(7,:),'linewidth',2)

    colormap(cmap)
    geoshow(LAT,LON,Z,'DisplayType','surface')
    c=colorbar;
    c.Ticks=[];
    set(gca,'fontsize',24)
    grid on
    box on
    ylabel('latitude')
    xlabel('longitude')
    title(titles{fi},'fontsize',24)

    outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/BSI_map_%s.png',fire_names{fi});
    saveas(f,outfilename)
end
