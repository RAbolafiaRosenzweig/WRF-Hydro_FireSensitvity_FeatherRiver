clc;clear all;close all;

%% first part of code creates spatially averaged time series of data (commented out b/c it only needs to run 1 time)

%% load in simulation grid:
geofilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/Fire_Adjusted/wrfinput_Feather_postifre_Grassland.nc';
LAT = ncread(geofilename,'XLAT');
LON = ncread(geofilename,'XLONG');

%% load in burn data:
bsi = ncread('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/geogrid_bsi.nc','bsi');
% idx_burned = find(bsi > 0);

%% load in shapefile data:
%load in study catchment shapes:
Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
ncatch = length(Catchments);

%load in fire shapes:
Study_fires = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Fires_of_interest/firep_2018_2021.shp');
shape_info_fires = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Fires_of_interest/firep_2018_2021.shp');
p1_fire = shape_info_fires.CoordinateReferenceSystem;
nfire = length(Study_fires);

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
    store_catch_idx_in{i} = store_r_c;

    %create a new index for only the burned pixels:
    catch_bsi = bsi(idx_in);
    idx_burned = find(catch_bsi > 1);
    idx_in_burned = idx_in(idx_burned);
    store_r_c = [];
    for j=1:length(idx_in_burned)
        [r,c] = ind2sub(size(LAT) , idx_in_burned(j));
        store_r_c = [store_r_c;r,c];
    end
    store_burned_catch_idx_in{i} = store_r_c;
end

%% get lsm outputs
baseline_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/Feather_Baseline/';
modified_param_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam/';
modified_param_GVF_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam_GVF/';
modified_param_GVF_VegClass_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/ModParam_GVF_VegClass/';
realistic_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Realistic/';


baseline_Albedo_all = [];
baseline_ECAN_all = [];
baseline_EDIR_all = [];
baseline_ETRAN_all = [];
baseline_ET_all = [];
baseline_SWE_all = [];
baseline_LAI_all = [];
baseline_UGDRNOFF_all = [];
baseline_SOIL_M1_all = [];
baseline_SOIL_M2_all = [];
baseline_SOIL_M3_all = [];
baseline_SOIL_M4_all = [];

modified_param_Albedo_all = [];
modified_param_ECAN_all = [];
modified_param_EDIR_all = [];
modified_param_ETRAN_all = [];
modified_param_ET_all = [];
modified_param_SWE_all = [];
modified_param_LAI_all = [];
modified_param_UGDRNOFF_all = [];
modified_param_SOIL_M1_all = [];
modified_param_SOIL_M2_all = [];
modified_param_SOIL_M3_all = [];
modified_param_SOIL_M4_all = [];

modified_param_GVF_Albedo_all = [];
modified_param_GVF_ECAN_all = [];
modified_param_GVF_EDIR_all = [];
modified_param_GVF_ETRAN_all = [];
modified_param_GVF_ET_all = [];
modified_param_GVF_SWE_all = [];
modified_param_GVF_LAI_all = [];
modified_param_GVF_UGDRNOFF_all = [];
modified_param_GVF_SOIL_M1_all = [];
modified_param_GVF_SOIL_M2_all = [];
modified_param_GVF_SOIL_M3_all = [];
modified_param_GVF_SOIL_M4_all = [];

modified_param_GVF_VegClass_Albedo_all = [];
modified_param_GVF_VegClass_ECAN_all = [];
modified_param_GVF_VegClass_EDIR_all = [];
modified_param_GVF_VegClass_ETRAN_all = [];
modified_param_GVF_VegClass_ET_all = [];
modified_param_GVF_VegClass_SWE_all = [];
modified_param_GVF_VegClass_LAI_all = [];
modified_param_GVF_VegClass_UGDRNOFF_all = [];
modified_param_GVF_VegClass_SOIL_M1_all = [];
modified_param_GVF_VegClass_SOIL_M2_all = [];
modified_param_GVF_VegClass_SOIL_M3_all = [];
modified_param_GVF_VegClass_SOIL_M4_all = [];

realistic_Albedo_all = [];
realistic_ECAN_all = [];
realistic_EDIR_all = [];
realistic_ETRAN_all = [];
realistic_ET_all = [];
realistic_SWE_all = [];
realistic_LAI_all = [];
realistic_UGDRNOFF_all = [];
realistic_SOIL_M1_all = [];
realistic_SOIL_M2_all = [];
realistic_SOIL_M3_all = [];
realistic_SOIL_M4_all = [];

layers = [0.1,0.3,0.6,1];
for WY=2000:2022
    infilename = sprintf('daily_WY%04d.nc',WY);

    %get baseline data
    baseline_Albedo = ncread([baseline_dir,infilename],'ALBEDO');
    baseline_ECAN = ncread([baseline_dir,infilename],'ECAN');
    baseline_EDIR = ncread([baseline_dir,infilename],'EDIR');
    baseline_ETRAN = ncread([baseline_dir,infilename],'ETRAN');
    baseline_ET = ncread([baseline_dir,infilename],'ET');
    baseline_SWE = ncread([baseline_dir,infilename],'SWE');
    baseline_LAI = ncread([baseline_dir,infilename],'LAI');
    baseline_UGDRNOFF = ncread([baseline_dir,infilename],'UGDRNOFF');
    baseline_SOIL_M = ncread([baseline_dir,infilename],'SOIL_M');
    
    %get modified param data
    modified_param_Albedo = ncread([modified_param_dir,infilename],'ALBEDO');
    modified_param_ECAN = ncread([modified_param_dir,infilename],'ECAN');
    modified_param_EDIR = ncread([modified_param_dir,infilename],'EDIR');
    modified_param_ETRAN = ncread([modified_param_dir,infilename],'ETRAN');
    modified_param_ET = ncread([modified_param_dir,infilename],'ET');
    modified_param_SWE = ncread([modified_param_dir,infilename],'SWE');
    modified_param_LAI = ncread([modified_param_dir,infilename],'LAI');
    modified_param_UGDRNOFF = ncread([modified_param_dir,infilename],'UGDRNOFF');
    modified_param_SOIL_M = ncread([modified_param_dir,infilename],'SOIL_M');
    
    %get modified param & veg class data - GRASS
    modified_param_GVF_Albedo = ncread([modified_param_GVF_dir,infilename],'ALBEDO');
    modified_param_GVF_ECAN = ncread([modified_param_GVF_dir,infilename],'ECAN');
    modified_param_GVF_EDIR = ncread([modified_param_GVF_dir,infilename],'EDIR');
    modified_param_GVF_ETRAN = ncread([modified_param_GVF_dir,infilename],'ETRAN');
    modified_param_GVF_ET = ncread([modified_param_GVF_dir,infilename],'ET');
    modified_param_GVF_SWE = ncread([modified_param_GVF_dir,infilename],'SWE');
    modified_param_GVF_LAI = ncread([modified_param_GVF_dir,infilename],'LAI');
    modified_param_GVF_UGDRNOFF = ncread([modified_param_GVF_dir,infilename],'UGDRNOFF');
    modified_param_GVF_SOIL_M = ncread([modified_param_GVF_dir,infilename],'SOIL_M');
    
    %get modified param & veg class data - BARE
    modified_param_GVF_VegClass_Albedo = ncread([modified_param_GVF_VegClass_dir,infilename],'ALBEDO');
    modified_param_GVF_VegClass_ECAN = ncread([modified_param_GVF_VegClass_dir,infilename],'ECAN');
    modified_param_GVF_VegClass_EDIR = ncread([modified_param_GVF_VegClass_dir,infilename],'EDIR');
    modified_param_GVF_VegClass_ETRAN = ncread([modified_param_GVF_VegClass_dir,infilename],'ETRAN');
    modified_param_GVF_VegClass_ET = ncread([modified_param_GVF_VegClass_dir,infilename],'ET');
    modified_param_GVF_VegClass_SWE = ncread([modified_param_GVF_VegClass_dir,infilename],'SWE');
    modified_param_GVF_VegClass_LAI = ncread([modified_param_GVF_VegClass_dir,infilename],'LAI');
    modified_param_GVF_VegClass_UGDRNOFF = ncread([modified_param_GVF_VegClass_dir,infilename],'UGDRNOFF');
    modified_param_GVF_VegClass_SOIL_M = ncread([modified_param_GVF_VegClass_dir,infilename],'SOIL_M');
    
    %get modified realistic
    realistic_Albedo = ncread([realistic_dir,infilename],'ALBEDO');
    realistic_ECAN = ncread([realistic_dir,infilename],'ECAN');
    realistic_EDIR = ncread([realistic_dir,infilename],'EDIR');
    realistic_ETRAN = ncread([realistic_dir,infilename],'ETRAN');
    realistic_ET = ncread([realistic_dir,infilename],'ET');
    realistic_SWE = ncread([realistic_dir,infilename],'SWE');
    realistic_LAI = ncread([realistic_dir,infilename],'LAI');
    realistic_UGDRNOFF = ncread([realistic_dir,infilename],'UGDRNOFF');
    realistic_SOIL_M = ncread([realistic_dir,infilename],'SOIL_M');
    
    
% %     store data
    baseline_Albedo_all = cat(3,baseline_Albedo_all,baseline_Albedo);
    baseline_ECAN_all = cat(3,baseline_ECAN_all,baseline_ECAN);
    baseline_EDIR_all = cat(3,baseline_EDIR_all,baseline_EDIR);
    baseline_ETRAN_all = cat(3,baseline_ETRAN_all,baseline_ETRAN);
    baseline_ET_all = cat(3,baseline_ET_all,baseline_ET);
    baseline_SWE_all = cat(3,baseline_SWE_all,baseline_SWE);
    baseline_LAI_all = cat(3,baseline_LAI_all,baseline_LAI);
    baseline_UGDRNOFF_all = cat(3,baseline_UGDRNOFF_all,baseline_UGDRNOFF);
    baseline_SOIL_M1_all = cat(3,baseline_SOIL_M1_all,squeeze(baseline_SOIL_M(:,1,:,:)));
    baseline_SOIL_M2_all = cat(3,baseline_SOIL_M2_all,squeeze(baseline_SOIL_M(:,2,:,:)));
    baseline_SOIL_M3_all = cat(3,baseline_SOIL_M3_all,squeeze(baseline_SOIL_M(:,3,:,:)));
    baseline_SOIL_M4_all = cat(3,baseline_SOIL_M4_all,squeeze(baseline_SOIL_M(:,4,:,:)));
    
    modified_param_Albedo_all = cat(3,modified_param_Albedo_all,modified_param_Albedo);
    modified_param_ECAN_all = cat(3,modified_param_ECAN_all,modified_param_ECAN);
    modified_param_EDIR_all = cat(3,modified_param_EDIR_all,modified_param_EDIR);
    modified_param_ETRAN_all = cat(3,modified_param_ETRAN_all,modified_param_ETRAN);
    modified_param_ET_all = cat(3,modified_param_ET_all,modified_param_ET);
    modified_param_SWE_all = cat(3,modified_param_SWE_all,modified_param_SWE);
    modified_param_LAI_all = cat(3,modified_param_LAI_all,modified_param_LAI);
    modified_param_UGDRNOFF_all = cat(3,modified_param_UGDRNOFF_all,modified_param_UGDRNOFF);
    modified_param_SOIL_M1_all = cat(3,modified_param_SOIL_M1_all,squeeze(modified_param_SOIL_M(:,1,:,:)));
    modified_param_SOIL_M2_all = cat(3,modified_param_SOIL_M2_all,squeeze(modified_param_SOIL_M(:,2,:,:)));
    modified_param_SOIL_M3_all = cat(3,modified_param_SOIL_M3_all,squeeze(modified_param_SOIL_M(:,3,:,:)));
    modified_param_SOIL_M4_all = cat(3,modified_param_SOIL_M4_all,squeeze(modified_param_SOIL_M(:,4,:,:)));

    modified_param_GVF_Albedo_all = cat(3,modified_param_GVF_Albedo_all,modified_param_GVF_Albedo);
    modified_param_GVF_ECAN_all = cat(3,modified_param_GVF_ECAN_all,modified_param_GVF_ECAN);
    modified_param_GVF_EDIR_all = cat(3,modified_param_GVF_EDIR_all,modified_param_GVF_EDIR);
    modified_param_GVF_ETRAN_all = cat(3,modified_param_GVF_ETRAN_all,modified_param_GVF_ETRAN);
    modified_param_GVF_ET_all = cat(3,modified_param_GVF_ET_all,modified_param_GVF_ET);
    modified_param_GVF_SWE_all = cat(3,modified_param_GVF_SWE_all,modified_param_GVF_SWE);
    modified_param_GVF_LAI_all = cat(3,modified_param_GVF_LAI_all,modified_param_GVF_LAI);
    modified_param_GVF_UGDRNOFF_all = cat(3,modified_param_GVF_UGDRNOFF_all,modified_param_GVF_UGDRNOFF);
    modified_param_GVF_SOIL_M1_all = cat(3,modified_param_GVF_SOIL_M1_all,squeeze(modified_param_GVF_SOIL_M(:,1,:,:)));
    modified_param_GVF_SOIL_M2_all = cat(3,modified_param_GVF_SOIL_M2_all,squeeze(modified_param_GVF_SOIL_M(:,2,:,:)));
    modified_param_GVF_SOIL_M3_all = cat(3,modified_param_GVF_SOIL_M3_all,squeeze(modified_param_GVF_SOIL_M(:,3,:,:)));
    modified_param_GVF_SOIL_M4_all = cat(3,modified_param_GVF_SOIL_M4_all,squeeze(modified_param_GVF_SOIL_M(:,4,:,:)));
    
    modified_param_GVF_VegClass_Albedo_all = cat(3,modified_param_GVF_VegClass_Albedo_all,modified_param_GVF_VegClass_Albedo);
    modified_param_GVF_VegClass_ECAN_all = cat(3,modified_param_GVF_VegClass_ECAN_all,modified_param_GVF_VegClass_ECAN);
    modified_param_GVF_VegClass_EDIR_all = cat(3,modified_param_GVF_VegClass_EDIR_all,modified_param_GVF_VegClass_EDIR);
    modified_param_GVF_VegClass_ETRAN_all = cat(3,modified_param_GVF_VegClass_ETRAN_all,modified_param_GVF_VegClass_ETRAN);
    modified_param_GVF_VegClass_ET_all = cat(3,modified_param_GVF_VegClass_ET_all,modified_param_GVF_VegClass_ET);
    modified_param_GVF_VegClass_SWE_all = cat(3,modified_param_GVF_VegClass_SWE_all,modified_param_GVF_VegClass_SWE);
    modified_param_GVF_VegClass_LAI_all = cat(3,modified_param_GVF_VegClass_LAI_all,modified_param_GVF_VegClass_LAI);
    modified_param_GVF_VegClass_UGDRNOFF_all = cat(3,modified_param_GVF_VegClass_UGDRNOFF_all,modified_param_GVF_VegClass_UGDRNOFF);
    modified_param_GVF_VegClass_SOIL_M1_all = cat(3,modified_param_GVF_VegClass_SOIL_M1_all,squeeze(modified_param_GVF_VegClass_SOIL_M(:,1,:,:)));
    modified_param_GVF_VegClass_SOIL_M2_all = cat(3,modified_param_GVF_VegClass_SOIL_M2_all,squeeze(modified_param_GVF_VegClass_SOIL_M(:,2,:,:)));
    modified_param_GVF_VegClass_SOIL_M3_all = cat(3,modified_param_GVF_VegClass_SOIL_M3_all,squeeze(modified_param_GVF_VegClass_SOIL_M(:,3,:,:)));
    modified_param_GVF_VegClass_SOIL_M4_all = cat(3,modified_param_GVF_VegClass_SOIL_M4_all,squeeze(modified_param_GVF_VegClass_SOIL_M(:,4,:,:)));
    
    realistic_Albedo_all = cat(3,realistic_Albedo_all,realistic_Albedo);
    realistic_ECAN_all = cat(3,realistic_ECAN_all,realistic_ECAN);
    realistic_EDIR_all = cat(3,realistic_EDIR_all,realistic_EDIR);
    realistic_ETRAN_all = cat(3,realistic_ETRAN_all,realistic_ETRAN);
    realistic_ET_all = cat(3,realistic_ET_all,realistic_ET);
    realistic_SWE_all = cat(3,realistic_SWE_all,realistic_SWE);
    realistic_LAI_all = cat(3,realistic_LAI_all,realistic_LAI);
    realistic_UGDRNOFF_all = cat(3,realistic_UGDRNOFF_all,realistic_UGDRNOFF);
    realistic_SOIL_M1_all = cat(3,realistic_SOIL_M1_all,squeeze(realistic_SOIL_M(:,1,:,:)));
    realistic_SOIL_M2_all = cat(3,realistic_SOIL_M2_all,squeeze(realistic_SOIL_M(:,2,:,:)));
    realistic_SOIL_M3_all = cat(3,realistic_SOIL_M3_all,squeeze(realistic_SOIL_M(:,3,:,:)));
    realistic_SOIL_M4_all = cat(3,realistic_SOIL_M4_all,squeeze(realistic_SOIL_M(:,4,:,:)));
end


%% store lsm outputs by each catchment
dates = [datenum([1999 10 1]):datenum([2022 9 30])]';
n = length(Catchments);
for i=1:n
    catch_idx = store_catch_idx_in{i};
    %% baseline data
    baseline_store_Albedo_catch=[];
    baseline_store_ECAN_catch = [];
    baseline_store_EDIR_catch = [];
    baseline_store_ETRAN_catch = [];
    baseline_store_ET_catch = [];
    baseline_store_SWE_catch = [];
    baseline_store_LAI_catch = [];
    baseline_store_UGDRNOFF_catch = [];
    baseline_store_SOIL_M1_catch = [];
    baseline_store_SOIL_M2_catch = [];
    baseline_store_SOIL_M3_catch = [];
    baseline_store_SOIL_M4_catch = [];
    for j=1:length(catch_idx)
        current_loc = catch_idx(j,:);
        baseline_Albedo_catch = squeeze(baseline_Albedo_all(current_loc(1),current_loc(2),:));
        baseline_ECAN_catch = squeeze(baseline_ECAN_all(current_loc(1),current_loc(2),:));
        baseline_EDIR_catch = squeeze(baseline_EDIR_all(current_loc(1),current_loc(2),:));
        baseline_ETRAN_catch = squeeze(baseline_ETRAN_all(current_loc(1),current_loc(2),:));
        baseline_ET_catch = squeeze(baseline_ET_all(current_loc(1),current_loc(2),:));
        baseline_SWE_catch = squeeze(baseline_SWE_all(current_loc(1),current_loc(2),:));
        baseline_LAI_catch = squeeze(baseline_LAI_all(current_loc(1),current_loc(2),:));
        baseline_UGDRNOFF_catch = squeeze(baseline_UGDRNOFF_all(current_loc(1),current_loc(2),:));
        baseline_SOIL_M1_catch = squeeze(baseline_SOIL_M1_all(current_loc(1),current_loc(2),:));
        baseline_SOIL_M2_catch = squeeze(baseline_SOIL_M2_all(current_loc(1),current_loc(2),:));
        baseline_SOIL_M3_catch = squeeze(baseline_SOIL_M3_all(current_loc(1),current_loc(2),:));
        baseline_SOIL_M4_catch = squeeze(baseline_SOIL_M4_all(current_loc(1),current_loc(2),:));

        baseline_store_Albedo_catch = [baseline_store_Albedo_catch,baseline_Albedo_catch];
        baseline_store_ECAN_catch = [baseline_store_ECAN_catch,baseline_ECAN_catch];
        baseline_store_EDIR_catch = [baseline_store_EDIR_catch,baseline_EDIR_catch];
        baseline_store_ETRAN_catch = [baseline_store_ETRAN_catch,baseline_ETRAN_catch];
        baseline_store_ET_catch = [baseline_store_ET_catch,baseline_ET_catch];
        baseline_store_SWE_catch = [baseline_store_SWE_catch,baseline_SWE_catch];
        baseline_store_LAI_catch = [baseline_store_LAI_catch,baseline_LAI_catch];
        baseline_store_UGDRNOFF_catch = [baseline_store_UGDRNOFF_catch,baseline_UGDRNOFF_catch];
        baseline_store_SOIL_M1_catch = [baseline_store_SOIL_M1_catch,baseline_SOIL_M1_catch];
        baseline_store_SOIL_M2_catch = [baseline_store_SOIL_M2_catch,baseline_SOIL_M2_catch];
        baseline_store_SOIL_M3_catch = [baseline_store_SOIL_M3_catch,baseline_SOIL_M3_catch];
        baseline_store_SOIL_M4_catch = [baseline_store_SOIL_M4_catch,baseline_SOIL_M4_catch];
    end
    mean_catch_Albedo = nanmean(baseline_store_Albedo_catch')';
    mean_catch_ECAN = nanmean(baseline_store_ECAN_catch')';
    mean_catch_EDIR = nanmean(baseline_store_EDIR_catch')';
    mean_catch_ETRAN = nanmean(baseline_store_ETRAN_catch')';
    mean_catch_ET= nanmean(baseline_store_ET_catch')';
    mean_catch_SWE = nanmean(baseline_store_SWE_catch')';
    mean_catch_LAI = nanmean(baseline_store_LAI_catch')';
    mean_catch_UGDRNOFF = nanmean(baseline_store_UGDRNOFF_catch')';
    mean_catch_SOIL_M1 = nanmean(baseline_store_SOIL_M1_catch')';
    mean_catch_SOIL_M2 = nanmean(baseline_store_SOIL_M2_catch')';
    mean_catch_SOIL_M3 = nanmean(baseline_store_SOIL_M3_catch')';
    mean_catch_SOIL_M4 = nanmean(baseline_store_SOIL_M4_catch')';
    %write outputs to table:
    out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
    colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
    outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/baseline_LSM_outputs_Catch_%d.mat',i);
    Catchment_outputs = array2table(out_data,'VariableNames',colnames);
    Catchment_outputs.Properties.VariableNames=colnames;
    save(outfilename,'Catchment_outputs','-v7.3');

    %% mod param data
    modified_store_Albedo_catch=[];
    modified_store_ECAN_catch = [];
    modified_store_EDIR_catch = [];
    modified_store_ETRAN_catch = [];
    modified_store_ET_catch = [];
    modified_store_SWE_catch = [];
    modified_store_LAI_catch = [];
    modified_store_UGDRNOFF_catch = [];
    modified_store_SOIL_M1_catch = [];
    modified_store_SOIL_M2_catch = [];
    modified_store_SOIL_M3_catch = [];
    modified_store_SOIL_M4_catch = [];
    for j=1:length(catch_idx)
        current_loc = catch_idx(j,:);
        modified_Albedo_catch = squeeze(modified_param_Albedo_all(current_loc(1),current_loc(2),:));
        modified_ECAN_catch = squeeze(modified_param_ECAN_all(current_loc(1),current_loc(2),:));
        modified_EDIR_catch = squeeze(modified_param_EDIR_all(current_loc(1),current_loc(2),:));
        modified_ETRAN_catch = squeeze(modified_param_ETRAN_all(current_loc(1),current_loc(2),:));
        modified_ET_catch = squeeze(modified_param_ET_all(current_loc(1),current_loc(2),:));
        modified_SWE_catch = squeeze(modified_param_SWE_all(current_loc(1),current_loc(2),:));
        modified_LAI_catch = squeeze(modified_param_LAI_all(current_loc(1),current_loc(2),:));
        modified_UGDRNOFF_catch = squeeze(modified_param_UGDRNOFF_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M1_catch = squeeze(modified_param_SOIL_M1_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M2_catch = squeeze(modified_param_SOIL_M2_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M3_catch = squeeze(modified_param_SOIL_M3_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M4_catch = squeeze(modified_param_SOIL_M4_all(current_loc(1),current_loc(2),:));

        modified_store_Albedo_catch = [modified_store_Albedo_catch,modified_Albedo_catch];
        modified_store_ECAN_catch = [modified_store_ECAN_catch,modified_ECAN_catch];
        modified_store_EDIR_catch = [modified_store_EDIR_catch,modified_EDIR_catch];
        modified_store_ETRAN_catch = [modified_store_ETRAN_catch,modified_ETRAN_catch];
        modified_store_ET_catch = [modified_store_ET_catch,modified_ET_catch];
        modified_store_SWE_catch = [modified_store_SWE_catch,modified_SWE_catch];
        modified_store_LAI_catch = [modified_store_LAI_catch,modified_LAI_catch];
        modified_store_UGDRNOFF_catch = [modified_store_UGDRNOFF_catch,modified_UGDRNOFF_catch];
        modified_store_SOIL_M1_catch = [modified_store_SOIL_M1_catch,modified_SOIL_M1_catch];
        modified_store_SOIL_M2_catch = [modified_store_SOIL_M2_catch,modified_SOIL_M2_catch];
        modified_store_SOIL_M3_catch = [modified_store_SOIL_M3_catch,modified_SOIL_M3_catch];
        modified_store_SOIL_M4_catch = [modified_store_SOIL_M4_catch,modified_SOIL_M4_catch];
    end
    mean_catch_Albedo = nanmean(modified_store_Albedo_catch')';
    mean_catch_ECAN = nanmean(modified_store_ECAN_catch')';
    mean_catch_EDIR = nanmean(modified_store_EDIR_catch')';
    mean_catch_ETRAN = nanmean(modified_store_ETRAN_catch')';
    mean_catch_ET= nanmean(modified_store_ET_catch')';
    mean_catch_SWE = nanmean(modified_store_SWE_catch')';
    mean_catch_LAI = nanmean(modified_store_LAI_catch')';
    mean_catch_UGDRNOFF = nanmean(modified_store_UGDRNOFF_catch')';
    mean_catch_SOIL_M1 = nanmean(modified_store_SOIL_M1_catch')';
    mean_catch_SOIL_M2 = nanmean(modified_store_SOIL_M2_catch')';
    mean_catch_SOIL_M3 = nanmean(modified_store_SOIL_M3_catch')';
    mean_catch_SOIL_M4 = nanmean(modified_store_SOIL_M4_catch')';
    %write outputs to table:
    out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
    colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
    outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/modified_param_LSM_outputs_Catch_%d.mat',i);
    Catchment_outputs = array2table(out_data,'VariableNames',colnames);
    Catchment_outputs.Properties.VariableNames=colnames;
    save(outfilename,'Catchment_outputs','-v7.3');
    
    %% mod param & GVF
    modified_store_Albedo_catch=[];
    modified_store_ECAN_catch = [];
    modified_store_EDIR_catch = [];
    modified_store_ETRAN_catch = [];
    modified_store_ET_catch = [];
    modified_store_SWE_catch = [];
    modified_store_LAI_catch = [];
    modified_store_UGDRNOFF_catch = [];
    modified_store_SOIL_M1_catch = [];
    modified_store_SOIL_M2_catch = [];
    modified_store_SOIL_M3_catch = [];
    modified_store_SOIL_M4_catch = [];
    for j=1:length(catch_idx)
        current_loc = catch_idx(j,:);
        modified_Albedo_catch = squeeze(modified_param_GVF_Albedo_all(current_loc(1),current_loc(2),:));
        modified_ECAN_catch = squeeze(modified_param_GVF_ECAN_all(current_loc(1),current_loc(2),:));
        modified_EDIR_catch = squeeze(modified_param_GVF_EDIR_all(current_loc(1),current_loc(2),:));
        modified_ETRAN_catch = squeeze(modified_param_GVF_ETRAN_all(current_loc(1),current_loc(2),:));
        modified_ET_catch = squeeze(modified_param_GVF_ET_all(current_loc(1),current_loc(2),:));
        modified_SWE_catch = squeeze(modified_param_GVF_SWE_all(current_loc(1),current_loc(2),:));
        modified_LAI_catch = squeeze(modified_param_GVF_LAI_all(current_loc(1),current_loc(2),:));
        modified_UGDRNOFF_catch = squeeze(modified_param_GVF_UGDRNOFF_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M1_catch = squeeze(modified_param_GVF_SOIL_M1_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M2_catch = squeeze(modified_param_GVF_SOIL_M2_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M3_catch = squeeze(modified_param_GVF_SOIL_M3_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M4_catch = squeeze(modified_param_GVF_SOIL_M4_all(current_loc(1),current_loc(2),:));

        modified_store_Albedo_catch = [modified_store_Albedo_catch,modified_Albedo_catch];
        modified_store_ECAN_catch = [modified_store_ECAN_catch,modified_ECAN_catch];
        modified_store_EDIR_catch = [modified_store_EDIR_catch,modified_EDIR_catch];
        modified_store_ETRAN_catch = [modified_store_ETRAN_catch,modified_ETRAN_catch];
        modified_store_ET_catch = [modified_store_ET_catch,modified_ET_catch];
        modified_store_SWE_catch = [modified_store_SWE_catch,modified_SWE_catch];
        modified_store_LAI_catch = [modified_store_LAI_catch,modified_LAI_catch];
        modified_store_UGDRNOFF_catch = [modified_store_UGDRNOFF_catch,modified_UGDRNOFF_catch];
        modified_store_SOIL_M1_catch = [modified_store_SOIL_M1_catch,modified_SOIL_M1_catch];
        modified_store_SOIL_M2_catch = [modified_store_SOIL_M2_catch,modified_SOIL_M2_catch];
        modified_store_SOIL_M3_catch = [modified_store_SOIL_M3_catch,modified_SOIL_M3_catch];
        modified_store_SOIL_M4_catch = [modified_store_SOIL_M4_catch,modified_SOIL_M4_catch];
    end
    mean_catch_Albedo = nanmean(modified_store_Albedo_catch')';
    mean_catch_ECAN = nanmean(modified_store_ECAN_catch')';
    mean_catch_EDIR = nanmean(modified_store_EDIR_catch')';
    mean_catch_ETRAN = nanmean(modified_store_ETRAN_catch')';
    mean_catch_ET= nanmean(modified_store_ET_catch')';
    mean_catch_SWE = nanmean(modified_store_SWE_catch')';
    mean_catch_LAI = nanmean(modified_store_LAI_catch')';
    mean_catch_UGDRNOFF = nanmean(modified_store_UGDRNOFF_catch')';
    mean_catch_SOIL_M1 = nanmean(modified_store_SOIL_M1_catch')';
    mean_catch_SOIL_M2 = nanmean(modified_store_SOIL_M2_catch')';
    mean_catch_SOIL_M3 = nanmean(modified_store_SOIL_M3_catch')';
    mean_catch_SOIL_M4 = nanmean(modified_store_SOIL_M4_catch')';
    %write outputs to table:
    out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
    colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
    outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/modified_param_GVF_LSM_outputs_Catch_%d.mat',i);
    Catchment_outputs = array2table(out_data,'VariableNames',colnames);
    Catchment_outputs.Properties.VariableNames=colnames;
    save(outfilename,'Catchment_outputs','-v7.3');
    
     %% mod param +GVF +veg class
    modified_store_Albedo_catch=[];
    modified_store_ECAN_catch = [];
    modified_store_EDIR_catch = [];
    modified_store_ETRAN_catch = [];
    modified_store_ET_catch = [];
    modified_store_SWE_catch = [];
    modified_store_LAI_catch = [];
    modified_store_UGDRNOFF_catch = [];
    modified_store_SOIL_M1_catch = [];
    modified_store_SOIL_M2_catch = [];
    modified_store_SOIL_M3_catch = [];
    modified_store_SOIL_M4_catch = [];
    for j=1:length(catch_idx)
        current_loc = catch_idx(j,:);
        modified_Albedo_catch = squeeze(modified_param_GVF_VegClass_Albedo_all(current_loc(1),current_loc(2),:));
        modified_ECAN_catch = squeeze(modified_param_GVF_VegClass_ECAN_all(current_loc(1),current_loc(2),:));
        modified_EDIR_catch = squeeze(modified_param_GVF_VegClass_EDIR_all(current_loc(1),current_loc(2),:));
        modified_ETRAN_catch = squeeze(modified_param_GVF_VegClass_ETRAN_all(current_loc(1),current_loc(2),:));
        modified_ET_catch = squeeze(modified_param_GVF_VegClass_ET_all(current_loc(1),current_loc(2),:));
        modified_SWE_catch = squeeze(modified_param_GVF_VegClass_SWE_all(current_loc(1),current_loc(2),:));
        modified_LAI_catch = squeeze(modified_param_GVF_VegClass_LAI_all(current_loc(1),current_loc(2),:));
        modified_UGDRNOFF_catch = squeeze(modified_param_GVF_VegClass_UGDRNOFF_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M1_catch = squeeze(modified_param_GVF_VegClass_SOIL_M1_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M2_catch = squeeze(modified_param_GVF_VegClass_SOIL_M2_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M3_catch = squeeze(modified_param_GVF_VegClass_SOIL_M3_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M4_catch = squeeze(modified_param_GVF_VegClass_SOIL_M4_all(current_loc(1),current_loc(2),:));

        modified_store_Albedo_catch = [modified_store_Albedo_catch,modified_Albedo_catch];
        modified_store_ECAN_catch = [modified_store_ECAN_catch,modified_ECAN_catch];
        modified_store_EDIR_catch = [modified_store_EDIR_catch,modified_EDIR_catch];
        modified_store_ETRAN_catch = [modified_store_ETRAN_catch,modified_ETRAN_catch];
        modified_store_ET_catch = [modified_store_ET_catch,modified_ET_catch];
        modified_store_SWE_catch = [modified_store_SWE_catch,modified_SWE_catch];
        modified_store_LAI_catch = [modified_store_LAI_catch,modified_LAI_catch];
        modified_store_UGDRNOFF_catch = [modified_store_UGDRNOFF_catch,modified_UGDRNOFF_catch];
        modified_store_SOIL_M1_catch = [modified_store_SOIL_M1_catch,modified_SOIL_M1_catch];
        modified_store_SOIL_M2_catch = [modified_store_SOIL_M2_catch,modified_SOIL_M2_catch];
        modified_store_SOIL_M3_catch = [modified_store_SOIL_M3_catch,modified_SOIL_M3_catch];
        modified_store_SOIL_M4_catch = [modified_store_SOIL_M4_catch,modified_SOIL_M4_catch];
    end
    mean_catch_Albedo = nanmean(modified_store_Albedo_catch')';
    mean_catch_ECAN = nanmean(modified_store_ECAN_catch')';
    mean_catch_EDIR = nanmean(modified_store_EDIR_catch')';
    mean_catch_ETRAN = nanmean(modified_store_ETRAN_catch')';
    mean_catch_ET= nanmean(modified_store_ET_catch')';
    mean_catch_SWE = nanmean(modified_store_SWE_catch')';
    mean_catch_LAI = nanmean(modified_store_LAI_catch')';
    mean_catch_UGDRNOFF = nanmean(modified_store_UGDRNOFF_catch')';
    mean_catch_SOIL_M1 = nanmean(modified_store_SOIL_M1_catch')';
    mean_catch_SOIL_M2 = nanmean(modified_store_SOIL_M2_catch')';
    mean_catch_SOIL_M3 = nanmean(modified_store_SOIL_M3_catch')';
    mean_catch_SOIL_M4 = nanmean(modified_store_SOIL_M4_catch')';
    %write outputs to table:
    out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
    colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
    outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/modified_param_GVF_VegClass_LSM_outputs_Catch_%d.mat',i);
    Catchment_outputs = array2table(out_data,'VariableNames',colnames);
    Catchment_outputs.Properties.VariableNames=colnames;
    save(outfilename,'Catchment_outputs','-v7.3');
    
    
    
    %% realistic adj
    modified_store_Albedo_catch=[];
    modified_store_ECAN_catch = [];
    modified_store_EDIR_catch = [];
    modified_store_ETRAN_catch = [];
    modified_store_ET_catch = [];
    modified_store_SWE_catch = [];
    modified_store_LAI_catch = [];
    modified_store_UGDRNOFF_catch = [];
    modified_store_SOIL_M1_catch = [];
    modified_store_SOIL_M2_catch = [];
    modified_store_SOIL_M3_catch = [];
    modified_store_SOIL_M4_catch = [];
    for j=1:length(catch_idx)
        current_loc = catch_idx(j,:);
        modified_Albedo_catch = squeeze(realistic_Albedo_all(current_loc(1),current_loc(2),:));
        modified_ECAN_catch = squeeze(realistic_ECAN_all(current_loc(1),current_loc(2),:));
        modified_EDIR_catch = squeeze(realistic_EDIR_all(current_loc(1),current_loc(2),:));
        modified_ETRAN_catch = squeeze(realistic_ETRAN_all(current_loc(1),current_loc(2),:));
        modified_ET_catch = squeeze(realistic_ET_all(current_loc(1),current_loc(2),:));
        modified_SWE_catch = squeeze(realistic_SWE_all(current_loc(1),current_loc(2),:));
        modified_LAI_catch = squeeze(realistic_LAI_all(current_loc(1),current_loc(2),:));
        modified_UGDRNOFF_catch = squeeze(realistic_UGDRNOFF_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M1_catch = squeeze(realistic_SOIL_M1_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M2_catch = squeeze(realistic_SOIL_M2_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M3_catch = squeeze(realistic_SOIL_M3_all(current_loc(1),current_loc(2),:));
        modified_SOIL_M4_catch = squeeze(realistic_SOIL_M4_all(current_loc(1),current_loc(2),:));

        modified_store_Albedo_catch = [modified_store_Albedo_catch,modified_Albedo_catch];
        modified_store_ECAN_catch = [modified_store_ECAN_catch,modified_ECAN_catch];
        modified_store_EDIR_catch = [modified_store_EDIR_catch,modified_EDIR_catch];
        modified_store_ETRAN_catch = [modified_store_ETRAN_catch,modified_ETRAN_catch];
        modified_store_ET_catch = [modified_store_ET_catch,modified_ET_catch];
        modified_store_SWE_catch = [modified_store_SWE_catch,modified_SWE_catch];
        modified_store_LAI_catch = [modified_store_LAI_catch,modified_LAI_catch];
        modified_store_UGDRNOFF_catch = [modified_store_UGDRNOFF_catch,modified_UGDRNOFF_catch];
        modified_store_SOIL_M1_catch = [modified_store_SOIL_M1_catch,modified_SOIL_M1_catch];
        modified_store_SOIL_M2_catch = [modified_store_SOIL_M2_catch,modified_SOIL_M2_catch];
        modified_store_SOIL_M3_catch = [modified_store_SOIL_M3_catch,modified_SOIL_M3_catch];
        modified_store_SOIL_M4_catch = [modified_store_SOIL_M4_catch,modified_SOIL_M4_catch];
    end
    mean_catch_Albedo = nanmean(modified_store_Albedo_catch')';
    mean_catch_ECAN = nanmean(modified_store_ECAN_catch')';
    mean_catch_EDIR = nanmean(modified_store_EDIR_catch')';
    mean_catch_ETRAN = nanmean(modified_store_ETRAN_catch')';
    mean_catch_ET= nanmean(modified_store_ET_catch')';
    mean_catch_SWE = nanmean(modified_store_SWE_catch')';
    mean_catch_LAI = nanmean(modified_store_LAI_catch')';
    mean_catch_UGDRNOFF = nanmean(modified_store_UGDRNOFF_catch')';
    mean_catch_SOIL_M1 = nanmean(modified_store_SOIL_M1_catch')';
    mean_catch_SOIL_M2 = nanmean(modified_store_SOIL_M2_catch')';
    mean_catch_SOIL_M3 = nanmean(modified_store_SOIL_M3_catch')';
    mean_catch_SOIL_M4 = nanmean(modified_store_SOIL_M4_catch')';
    %write outputs to table:
    out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
    colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
    outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/realistic_LSM_outputs_Catch_%d.mat',i);
    Catchment_outputs = array2table(out_data,'VariableNames',colnames);
    Catchment_outputs.Properties.VariableNames=colnames;
    save(outfilename,'Catchment_outputs','-v7.3');
end

% % % % % %% store lsm outputs by each catchment only averaged over burned pixels
% % % % % dates = [datenum([1999 10 1]):datenum([2022 9 30])]';
% % % % % n = length(Catchments);
% % % % % for i=1:n-1
% % % % %     burned_catch_idx = store_burned_catch_idx_in{i};
% % % % %     
% % % % %     %% baseline data
% % % % %     baseline_store_Albedo_catch=[];
% % % % %     baseline_store_ECAN_catch = [];
% % % % %     baseline_store_EDIR_catch = [];
% % % % %     baseline_store_ETRAN_catch = [];
% % % % %     baseline_store_ET_catch = [];
% % % % %     baseline_store_SWE_catch = [];
% % % % %     baseline_store_LAI_catch = [];
% % % % %     baseline_store_UGDRNOFF_catch = [];
% % % % %     baseline_store_SOIL_M1_catch = [];
% % % % %     baseline_store_SOIL_M2_catch = [];
% % % % %     baseline_store_SOIL_M3_catch = [];
% % % % %     baseline_store_SOIL_M4_catch = [];
% % % % %     for j=1:length(burned_catch_idx)
% % % % %         current_loc = burned_catch_idx(j,:);
% % % % %         baseline_Albedo_catch = squeeze(baseline_Albedo_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_ECAN_catch = squeeze(baseline_ECAN_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_EDIR_catch = squeeze(baseline_EDIR_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_ETRAN_catch = squeeze(baseline_ETRAN_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_ET_catch = squeeze(baseline_ET_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_SWE_catch = squeeze(baseline_SWE_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_LAI_catch = squeeze(baseline_LAI_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_UGDRNOFF_catch = squeeze(baseline_UGDRNOFF_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_SOIL_M1_catch = squeeze(baseline_SOIL_M1_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_SOIL_M2_catch = squeeze(baseline_SOIL_M2_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_SOIL_M3_catch = squeeze(baseline_SOIL_M3_all(current_loc(1),current_loc(2),:));
% % % % %         baseline_SOIL_M4_catch = squeeze(baseline_SOIL_M4_all(current_loc(1),current_loc(2),:));
% % % % % 
% % % % %         baseline_store_Albedo_catch = [baseline_store_Albedo_catch,baseline_Albedo_catch];
% % % % %         baseline_store_ECAN_catch = [baseline_store_ECAN_catch,baseline_ECAN_catch];
% % % % %         baseline_store_EDIR_catch = [baseline_store_EDIR_catch,baseline_EDIR_catch];
% % % % %         baseline_store_ETRAN_catch = [baseline_store_ETRAN_catch,baseline_ETRAN_catch];
% % % % %         baseline_store_ET_catch = [baseline_store_ET_catch,baseline_ET_catch];
% % % % %         baseline_store_SWE_catch = [baseline_store_SWE_catch,baseline_SWE_catch];
% % % % %         baseline_store_LAI_catch = [baseline_store_LAI_catch,baseline_LAI_catch];
% % % % %         baseline_store_UGDRNOFF_catch = [baseline_store_UGDRNOFF_catch,baseline_UGDRNOFF_catch];
% % % % %         baseline_store_SOIL_M1_catch = [baseline_store_SOIL_M1_catch,baseline_SOIL_M1_catch];
% % % % %         baseline_store_SOIL_M2_catch = [baseline_store_SOIL_M2_catch,baseline_SOIL_M2_catch];
% % % % %         baseline_store_SOIL_M3_catch = [baseline_store_SOIL_M3_catch,baseline_SOIL_M3_catch];
% % % % %         baseline_store_SOIL_M4_catch = [baseline_store_SOIL_M4_catch,baseline_SOIL_M4_catch];
% % % % %     end
% % % % %     mean_catch_Albedo = nanmean(baseline_store_Albedo_catch')';
% % % % %     mean_catch_ECAN = nanmean(baseline_store_ECAN_catch')';
% % % % %     mean_catch_EDIR = nanmean(baseline_store_EDIR_catch')';
% % % % %     mean_catch_ETRAN = nanmean(baseline_store_ETRAN_catch')';
% % % % %     mean_catch_ET= nanmean(baseline_store_ET_catch')';
% % % % %     mean_catch_SWE = nanmean(baseline_store_SWE_catch')';
% % % % %     mean_catch_LAI = nanmean(baseline_store_LAI_catch')';
% % % % %     mean_catch_UGDRNOFF = nanmean(baseline_store_UGDRNOFF_catch')';
% % % % %     mean_catch_SOIL_M1 = nanmean(baseline_store_SOIL_M1_catch')';
% % % % %     mean_catch_SOIL_M2 = nanmean(baseline_store_SOIL_M2_catch')';
% % % % %     mean_catch_SOIL_M3 = nanmean(baseline_store_SOIL_M3_catch')';
% % % % %     mean_catch_SOIL_M4 = nanmean(baseline_store_SOIL_M4_catch')';
% % % % %     
% % % % %     %write outputs to table:
% % % % %     out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
% % % % %     colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
% % % % %     outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/baseline_LSM_outputs_BurnedPixels_Catch_%d.mat',i);
% % % % %     Catchment_outputs = array2table(out_data,'VariableNames',colnames);
% % % % %     Catchment_outputs.Properties.VariableNames=colnames;
% % % % %     save(outfilename,'Catchment_outputs','-v7.3');
% % % % %     
% % % % %     %% mod param
% % % % %     modified_store_Albedo_catch=[];
% % % % %     modified_store_ECAN_catch = [];
% % % % %     modified_store_EDIR_catch = [];
% % % % %     modified_store_ETRAN_catch = [];
% % % % %     modified_store_ET_catch = [];
% % % % %     modified_store_SWE_catch = [];
% % % % %     modified_store_LAI_catch = [];
% % % % %     modified_store_UGDRNOFF_catch = [];
% % % % %     modified_store_SOIL_M1_catch = [];
% % % % %     modified_store_SOIL_M2_catch = [];
% % % % %     modified_store_SOIL_M3_catch = [];
% % % % %     modified_store_SOIL_M4_catch = [];
% % % % %     for j=1:length(burned_catch_idx)
% % % % %         current_loc = burned_catch_idx(j,:);
% % % % %         modified_Albedo_catch = squeeze(modified_param_Albedo_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ECAN_catch = squeeze(modified_param_ECAN_all(current_loc(1),current_loc(2),:));
% % % % %         modified_EDIR_catch = squeeze(modified_param_EDIR_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ETRAN_catch = squeeze(modified_param_ETRAN_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ET_catch = squeeze(modified_param_ET_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SWE_catch = squeeze(modified_param_SWE_all(current_loc(1),current_loc(2),:));
% % % % %         modified_LAI_catch = squeeze(modified_param_LAI_all(current_loc(1),current_loc(2),:));
% % % % %         modified_UGDRNOFF_catch = squeeze(modified_param_UGDRNOFF_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M1_catch = squeeze(modified_param_SOIL_M1_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M2_catch = squeeze(modified_param_SOIL_M2_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M3_catch = squeeze(modified_param_SOIL_M3_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M4_catch = squeeze(modified_param_SOIL_M4_all(current_loc(1),current_loc(2),:));
% % % % % 
% % % % %         modified_store_Albedo_catch = [modified_store_Albedo_catch,modified_Albedo_catch];
% % % % %         modified_store_ECAN_catch = [modified_store_ECAN_catch,modified_ECAN_catch];
% % % % %         modified_store_EDIR_catch = [modified_store_EDIR_catch,modified_EDIR_catch];
% % % % %         modified_store_ETRAN_catch = [modified_store_ETRAN_catch,modified_ETRAN_catch];
% % % % %         modified_store_ET_catch = [modified_store_ET_catch,modified_ET_catch];
% % % % %         modified_store_SWE_catch = [modified_store_SWE_catch,modified_SWE_catch];
% % % % %         modified_store_LAI_catch = [modified_store_LAI_catch,modified_LAI_catch];
% % % % %         modified_store_UGDRNOFF_catch = [modified_store_UGDRNOFF_catch,modified_UGDRNOFF_catch];
% % % % %         modified_store_SOIL_M1_catch = [modified_store_SOIL_M1_catch,modified_SOIL_M1_catch];
% % % % %         modified_store_SOIL_M2_catch = [modified_store_SOIL_M2_catch,modified_SOIL_M2_catch];
% % % % %         modified_store_SOIL_M3_catch = [modified_store_SOIL_M3_catch,modified_SOIL_M3_catch];
% % % % %         modified_store_SOIL_M4_catch = [modified_store_SOIL_M4_catch,modified_SOIL_M4_catch];
% % % % %     end
% % % % %     mean_catch_Albedo = nanmean(modified_store_Albedo_catch')';
% % % % %     mean_catch_ECAN = nanmean(modified_store_ECAN_catch')';
% % % % %     mean_catch_EDIR = nanmean(modified_store_EDIR_catch')';
% % % % %     mean_catch_ETRAN = nanmean(modified_store_ETRAN_catch')';
% % % % %     mean_catch_ET= nanmean(modified_store_ET_catch')';
% % % % %     mean_catch_SWE = nanmean(modified_store_SWE_catch')';
% % % % %     mean_catch_LAI = nanmean(modified_store_LAI_catch')';
% % % % %     mean_catch_UGDRNOFF = nanmean(modified_store_UGDRNOFF_catch')';
% % % % %     mean_catch_SOIL_M1 = nanmean(modified_store_SOIL_M1_catch')';
% % % % %     mean_catch_SOIL_M2 = nanmean(modified_store_SOIL_M2_catch')';
% % % % %     mean_catch_SOIL_M3 = nanmean(modified_store_SOIL_M3_catch')';
% % % % %     mean_catch_SOIL_M4 = nanmean(modified_store_SOIL_M4_catch')';
% % % % %     %write outputs to table:
% % % % %     out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
% % % % %     colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
% % % % %     outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/modified_param_LSM_outputs_BurnedPixels_Catch_%d.mat',i);
% % % % %     Catchment_outputs = array2table(out_data,'VariableNames',colnames);
% % % % %     Catchment_outputs.Properties.VariableNames=colnames;
% % % % %     save(outfilename,'Catchment_outputs','-v7.3');
% % % % %     
% % % % %     %% mod param +GVF
% % % % %     modified_store_Albedo_catch=[];
% % % % %     modified_store_ECAN_catch = [];
% % % % %     modified_store_EDIR_catch = [];
% % % % %     modified_store_ETRAN_catch = [];
% % % % %     modified_store_ET_catch = [];
% % % % %     modified_store_SWE_catch = [];
% % % % %     modified_store_LAI_catch = [];
% % % % %     modified_store_UGDRNOFF_catch = [];
% % % % %     modified_store_SOIL_M1_catch = [];
% % % % %     modified_store_SOIL_M2_catch = [];
% % % % %     modified_store_SOIL_M3_catch = [];
% % % % %     modified_store_SOIL_M4_catch = [];
% % % % %     for j=1:length(burned_catch_idx)
% % % % %         current_loc = burned_catch_idx(j,:);
% % % % %         modified_Albedo_catch = squeeze(modified_param_GVF_Albedo_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ECAN_catch = squeeze(modified_param_GVF_ECAN_all(current_loc(1),current_loc(2),:));
% % % % %         modified_EDIR_catch = squeeze(modified_param_GVF_EDIR_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ETRAN_catch = squeeze(modified_param_GVF_ETRAN_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ET_catch = squeeze(modified_param_GVF_ET_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SWE_catch = squeeze(modified_param_GVF_SWE_all(current_loc(1),current_loc(2),:));
% % % % %         modified_LAI_catch = squeeze(modified_param_GVF_LAI_all(current_loc(1),current_loc(2),:));
% % % % %         modified_UGDRNOFF_catch = squeeze(modified_param_GVF_UGDRNOFF_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M1_catch = squeeze(modified_param_GVF_SOIL_M1_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M2_catch = squeeze(modified_param_GVF_SOIL_M2_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M3_catch = squeeze(modified_param_GVF_SOIL_M3_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M4_catch = squeeze(modified_param_GVF_SOIL_M4_all(current_loc(1),current_loc(2),:));
% % % % % 
% % % % %         modified_store_Albedo_catch = [modified_store_Albedo_catch,modified_Albedo_catch];
% % % % %         modified_store_ECAN_catch = [modified_store_ECAN_catch,modified_ECAN_catch];
% % % % %         modified_store_EDIR_catch = [modified_store_EDIR_catch,modified_EDIR_catch];
% % % % %         modified_store_ETRAN_catch = [modified_store_ETRAN_catch,modified_ETRAN_catch];
% % % % %         modified_store_ET_catch = [modified_store_ET_catch,modified_ET_catch];
% % % % %         modified_store_SWE_catch = [modified_store_SWE_catch,modified_SWE_catch];
% % % % %         modified_store_LAI_catch = [modified_store_LAI_catch,modified_LAI_catch];
% % % % %         modified_store_UGDRNOFF_catch = [modified_store_UGDRNOFF_catch,modified_UGDRNOFF_catch];
% % % % %         modified_store_SOIL_M1_catch = [modified_store_SOIL_M1_catch,modified_SOIL_M1_catch];
% % % % %         modified_store_SOIL_M2_catch = [modified_store_SOIL_M2_catch,modified_SOIL_M2_catch];
% % % % %         modified_store_SOIL_M3_catch = [modified_store_SOIL_M3_catch,modified_SOIL_M3_catch];
% % % % %         modified_store_SOIL_M4_catch = [modified_store_SOIL_M4_catch,modified_SOIL_M4_catch];
% % % % %     end
% % % % %     mean_catch_Albedo = nanmean(modified_store_Albedo_catch')';
% % % % %     mean_catch_ECAN = nanmean(modified_store_ECAN_catch')';
% % % % %     mean_catch_EDIR = nanmean(modified_store_EDIR_catch')';
% % % % %     mean_catch_ETRAN = nanmean(modified_store_ETRAN_catch')';
% % % % %     mean_catch_ET= nanmean(modified_store_ET_catch')';
% % % % %     mean_catch_SWE = nanmean(modified_store_SWE_catch')';
% % % % %     mean_catch_LAI = nanmean(modified_store_LAI_catch')';
% % % % %     mean_catch_UGDRNOFF = nanmean(modified_store_UGDRNOFF_catch')';
% % % % %     mean_catch_SOIL_M1 = nanmean(modified_store_SOIL_M1_catch')';
% % % % %     mean_catch_SOIL_M2 = nanmean(modified_store_SOIL_M2_catch')';
% % % % %     mean_catch_SOIL_M3 = nanmean(modified_store_SOIL_M3_catch')';
% % % % %     mean_catch_SOIL_M4 = nanmean(modified_store_SOIL_M4_catch')';
% % % % %     %write outputs to table:
% % % % %     out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
% % % % %     colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
% % % % %     outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/modified_param_GVF_LSM_outputs_BurnedPixels_Catch_%d.mat',i);
% % % % %     Catchment_outputs = array2table(out_data,'VariableNames',colnames);
% % % % %     Catchment_outputs.Properties.VariableNames=colnames;
% % % % %     save(outfilename,'Catchment_outputs','-v7.3');
% % % % %     
% % % % %     %% mod param +gvf + veg class
% % % % %     modified_store_Albedo_catch=[];
% % % % %     modified_store_ECAN_catch = [];
% % % % %     modified_store_EDIR_catch = [];
% % % % %     modified_store_ETRAN_catch = [];
% % % % %     modified_store_ET_catch = [];
% % % % %     modified_store_SWE_catch = [];
% % % % %     modified_store_LAI_catch = [];
% % % % %     modified_store_UGDRNOFF_catch = [];
% % % % %     modified_store_SOIL_M1_catch = [];
% % % % %     modified_store_SOIL_M2_catch = [];
% % % % %     modified_store_SOIL_M3_catch = [];
% % % % %     modified_store_SOIL_M4_catch = [];
% % % % %     for j=1:length(burned_catch_idx)
% % % % %         current_loc = burned_catch_idx(j,:);
% % % % %         modified_Albedo_catch = squeeze(modified_param_GVF_VegClass_Albedo_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ECAN_catch = squeeze(modified_param_GVF_VegClass_ECAN_all(current_loc(1),current_loc(2),:));
% % % % %         modified_EDIR_catch = squeeze(modified_param_GVF_VegClass_EDIR_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ETRAN_catch = squeeze(modified_param_GVF_VegClass_ETRAN_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ET_catch = squeeze(modified_param_GVF_VegClass_ET_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SWE_catch = squeeze(modified_param_GVF_VegClass_SWE_all(current_loc(1),current_loc(2),:));
% % % % %         modified_LAI_catch = squeeze(modified_param_GVF_VegClass_LAI_all(current_loc(1),current_loc(2),:));
% % % % %         modified_UGDRNOFF_catch = squeeze(modified_param_GVF_VegClass_UGDRNOFF_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M1_catch = squeeze(modified_param_GVF_VegClass_SOIL_M1_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M2_catch = squeeze(modified_param_GVF_VegClass_SOIL_M2_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M3_catch = squeeze(modified_param_GVF_VegClass_SOIL_M3_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M4_catch = squeeze(modified_param_GVF_VegClass_SOIL_M4_all(current_loc(1),current_loc(2),:));
% % % % % 
% % % % %         modified_store_Albedo_catch = [modified_store_Albedo_catch,modified_Albedo_catch];
% % % % %         modified_store_ECAN_catch = [modified_store_ECAN_catch,modified_ECAN_catch];
% % % % %         modified_store_EDIR_catch = [modified_store_EDIR_catch,modified_EDIR_catch];
% % % % %         modified_store_ETRAN_catch = [modified_store_ETRAN_catch,modified_ETRAN_catch];
% % % % %         modified_store_ET_catch = [modified_store_ET_catch,modified_ET_catch];
% % % % %         modified_store_SWE_catch = [modified_store_SWE_catch,modified_SWE_catch];
% % % % %         modified_store_LAI_catch = [modified_store_LAI_catch,modified_LAI_catch];
% % % % %         modified_store_UGDRNOFF_catch = [modified_store_UGDRNOFF_catch,modified_UGDRNOFF_catch];
% % % % %         modified_store_SOIL_M1_catch = [modified_store_SOIL_M1_catch,modified_SOIL_M1_catch];
% % % % %         modified_store_SOIL_M2_catch = [modified_store_SOIL_M2_catch,modified_SOIL_M2_catch];
% % % % %         modified_store_SOIL_M3_catch = [modified_store_SOIL_M3_catch,modified_SOIL_M3_catch];
% % % % %         modified_store_SOIL_M4_catch = [modified_store_SOIL_M4_catch,modified_SOIL_M4_catch];
% % % % %     end
% % % % %     mean_catch_Albedo = nanmean(modified_store_Albedo_catch')';
% % % % %     mean_catch_ECAN = nanmean(modified_store_ECAN_catch')';
% % % % %     mean_catch_EDIR = nanmean(modified_store_EDIR_catch')';
% % % % %     mean_catch_ETRAN = nanmean(modified_store_ETRAN_catch')';
% % % % %     mean_catch_ET= nanmean(modified_store_ET_catch')';
% % % % %     mean_catch_SWE = nanmean(modified_store_SWE_catch')';
% % % % %     mean_catch_LAI = nanmean(modified_store_LAI_catch')';
% % % % %     mean_catch_UGDRNOFF = nanmean(modified_store_UGDRNOFF_catch')';
% % % % %     mean_catch_SOIL_M1 = nanmean(modified_store_SOIL_M1_catch')';
% % % % %     mean_catch_SOIL_M2 = nanmean(modified_store_SOIL_M2_catch')';
% % % % %     mean_catch_SOIL_M3 = nanmean(modified_store_SOIL_M3_catch')';
% % % % %     mean_catch_SOIL_M4 = nanmean(modified_store_SOIL_M4_catch')';
% % % % %     %write outputs to table:
% % % % %     out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
% % % % %     colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
% % % % %     outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/modified_param_GVF_VegClass_LSM_outputs_BurnedPixels_Catch_%d.mat',i);
% % % % %     Catchment_outputs = array2table(out_data,'VariableNames',colnames);
% % % % %     Catchment_outputs.Properties.VariableNames=colnames;
% % % % %     save(outfilename,'Catchment_outputs','-v7.3');
% % % % %     
% % % % %     %% realistic adj
% % % % %     modified_store_Albedo_catch=[];
% % % % %     modified_store_ECAN_catch = [];
% % % % %     modified_store_EDIR_catch = [];
% % % % %     modified_store_ETRAN_catch = [];
% % % % %     modified_store_ET_catch = [];
% % % % %     modified_store_SWE_catch = [];
% % % % %     modified_store_LAI_catch = [];
% % % % %     modified_store_UGDRNOFF_catch = [];
% % % % %     modified_store_SOIL_M1_catch = [];
% % % % %     modified_store_SOIL_M2_catch = [];
% % % % %     modified_store_SOIL_M3_catch = [];
% % % % %     modified_store_SOIL_M4_catch = [];
% % % % %     for j=1:length(burned_catch_idx)
% % % % %         current_loc = burned_catch_idx(j,:);
% % % % %         modified_Albedo_catch = squeeze(realistic_Albedo_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ECAN_catch = squeeze(realistic_ECAN_all(current_loc(1),current_loc(2),:));
% % % % %         modified_EDIR_catch = squeeze(realistic_EDIR_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ETRAN_catch = squeeze(realistic_ETRAN_all(current_loc(1),current_loc(2),:));
% % % % %         modified_ET_catch = squeeze(realistic_ET_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SWE_catch = squeeze(realistic_SWE_all(current_loc(1),current_loc(2),:));
% % % % %         modified_LAI_catch = squeeze(realistic_LAI_all(current_loc(1),current_loc(2),:));
% % % % %         modified_UGDRNOFF_catch = squeeze(realistic_UGDRNOFF_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M1_catch = squeeze(realistic_SOIL_M1_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M2_catch = squeeze(realistic_SOIL_M2_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M3_catch = squeeze(realistic_SOIL_M3_all(current_loc(1),current_loc(2),:));
% % % % %         modified_SOIL_M4_catch = squeeze(realistic_SOIL_M4_all(current_loc(1),current_loc(2),:));
% % % % % 
% % % % %         modified_store_Albedo_catch = [modified_store_Albedo_catch,modified_Albedo_catch];
% % % % %         modified_store_ECAN_catch = [modified_store_ECAN_catch,modified_ECAN_catch];
% % % % %         modified_store_EDIR_catch = [modified_store_EDIR_catch,modified_EDIR_catch];
% % % % %         modified_store_ETRAN_catch = [modified_store_ETRAN_catch,modified_ETRAN_catch];
% % % % %         modified_store_ET_catch = [modified_store_ET_catch,modified_ET_catch];
% % % % %         modified_store_SWE_catch = [modified_store_SWE_catch,modified_SWE_catch];
% % % % %         modified_store_LAI_catch = [modified_store_LAI_catch,modified_LAI_catch];
% % % % %         modified_store_UGDRNOFF_catch = [modified_store_UGDRNOFF_catch,modified_UGDRNOFF_catch];
% % % % %         modified_store_SOIL_M1_catch = [modified_store_SOIL_M1_catch,modified_SOIL_M1_catch];
% % % % %         modified_store_SOIL_M2_catch = [modified_store_SOIL_M2_catch,modified_SOIL_M2_catch];
% % % % %         modified_store_SOIL_M3_catch = [modified_store_SOIL_M3_catch,modified_SOIL_M3_catch];
% % % % %         modified_store_SOIL_M4_catch = [modified_store_SOIL_M4_catch,modified_SOIL_M4_catch];
% % % % %     end
% % % % %     mean_catch_Albedo = nanmean(modified_store_Albedo_catch')';
% % % % %     mean_catch_ECAN = nanmean(modified_store_ECAN_catch')';
% % % % %     mean_catch_EDIR = nanmean(modified_store_EDIR_catch')';
% % % % %     mean_catch_ETRAN = nanmean(modified_store_ETRAN_catch')';
% % % % %     mean_catch_ET= nanmean(modified_store_ET_catch')';
% % % % %     mean_catch_SWE = nanmean(modified_store_SWE_catch')';
% % % % %     mean_catch_LAI = nanmean(modified_store_LAI_catch')';
% % % % %     mean_catch_UGDRNOFF = nanmean(modified_store_UGDRNOFF_catch')';
% % % % %     mean_catch_SOIL_M1 = nanmean(modified_store_SOIL_M1_catch')';
% % % % %     mean_catch_SOIL_M2 = nanmean(modified_store_SOIL_M2_catch')';
% % % % %     mean_catch_SOIL_M3 = nanmean(modified_store_SOIL_M3_catch')';
% % % % %     mean_catch_SOIL_M4 = nanmean(modified_store_SOIL_M4_catch')';
% % % % %     %write outputs to table:
% % % % %     out_data = [dates,mean_catch_Albedo,mean_catch_ECAN,mean_catch_ETRAN,mean_catch_EDIR , mean_catch_ET , mean_catch_SWE , mean_catch_LAI,mean_catch_UGDRNOFF , mean_catch_SOIL_M1,mean_catch_SOIL_M2,mean_catch_SOIL_M3,mean_catch_SOIL_M4];
% % % % %     colnames = {'date','albedo','ecan','etran','edir','et','swe','lai','ugdrnoff','soilm1','soilm2','soilm3','soilm4'};
% % % % %     outfilename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/realistic_LSM_outputs_BurnedPixels_Catch_%d.mat',i);
% % % % %     Catchment_outputs = array2table(out_data,'VariableNames',colnames);
% % % % %     Catchment_outputs.Properties.VariableNames=colnames;
% % % % %     save(outfilename,'Catchment_outputs','-v7.3');
% % % % % end

%% this 2nd part of the code compares the time series data created in part 1 (for catchment averages)
outdir = '/Users/abolafia/ASO_Fire/Plots/PaperPlots/Catchment_WB_comparissons/';
if exist(outdir,'dir') == 0
    system(['mkdir -p ',outdir])
end
data_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/';
catch_names = {'East_Branch_NF_Feather','Mid_Fork_Feather','NF_Feather','Upper_Yuba'};
for c=2:2
    
    %% baseline data
    baseline_data_catch = load([data_dir,sprintf('baseline_LSM_outputs_Catch_%d.mat',c)]);
    
    baseline_data_catch = baseline_data_catch.Catchment_outputs;
    dates = double(baseline_data_catch.date);
    datevecs = datevec(dates);
    baseline_catch_albedo = baseline_data_catch.albedo;
    baseline_catch_ecan = baseline_data_catch.ecan;
    baseline_catch_etran = baseline_data_catch.etran;
    baseline_catch_edir = baseline_data_catch.edir;
    baseline_catch_et = baseline_data_catch.et;
    baseline_catch_swe = baseline_data_catch.swe;
    baseline_catch_lai = baseline_data_catch.lai;
    baseline_catch_ugdrnoff = baseline_data_catch.ugdrnoff;
    baseline_catch_soilm1 = baseline_data_catch.soilm1;
    baseline_catch_soilm2 = baseline_data_catch.soilm2;
    baseline_catch_soilm3 = baseline_data_catch.soilm3;
    baseline_catch_soilm4 = baseline_data_catch.soilm4;
    baseline_catch_soilm_total = baseline_catch_soilm1.*(0.1/2) + baseline_catch_soilm2.*(0.3/2) + baseline_catch_soilm3.*(0.6/2) + baseline_catch_soilm4.*(1.0/2); 
    
    %convert accumulated variables to daily (mm/day)
    baseline_catch_ecan = diff(baseline_catch_ecan);
    baseline_catch_etran = diff(baseline_catch_etran);
    baseline_catch_edir = diff(baseline_catch_edir);
    baseline_catch_et = diff(baseline_catch_et);
    baseline_catch_ugdrnoff = diff(baseline_catch_ugdrnoff);
    %remove day 1 from other var:
    dates(1) = [];
    datevecs(1,:)=[];
    baseline_catch_albedo(1)=[];
    baseline_catch_swe(1)=[];
    baseline_catch_lai(1)=[];
    baseline_catch_soilm1(1)=[];
    baseline_catch_soilm2(1)=[];
    baseline_catch_soilm3(1)=[];
    baseline_catch_soilm4(1)=[];
    baseline_catch_soilm_total(1)=[];

    %aggregate to monthly average:
    [u,~,j] = unique(datevecs(:,1:2),'rows','stable');
    monthly_dates = datenum([u,ones(length(u),1)]);
    baseline_catch_monthly_albedo = accumarray(j,baseline_catch_albedo,[],@nanmean);
    baseline_catch_monthly_ecan = accumarray(j,baseline_catch_ecan,[],@nanmean);
    baseline_catch_monthly_etran = accumarray(j,baseline_catch_etran,[],@nanmean);
    baseline_catch_monthly_edir = accumarray(j,baseline_catch_edir,[],@nanmean);
    baseline_catch_monthly_et = accumarray(j,baseline_catch_et,[],@nanmean);
    baseline_catch_monthly_swe = accumarray(j,baseline_catch_swe,[],@nanmean);
    baseline_catch_monthly_lai = accumarray(j,baseline_catch_lai,[],@nanmean);
    baseline_catch_monthly_ugdrnoff = accumarray(j,baseline_catch_ugdrnoff,[],@nanmean);
    baseline_catch_monthly_soilm1 = accumarray(j,baseline_catch_soilm1,[],@nanmean);
    baseline_catch_monthly_soilm2 = accumarray(j,baseline_catch_soilm2,[],@nanmean);
    baseline_catch_monthly_soilm3 = accumarray(j,baseline_catch_soilm3,[],@nanmean);
    baseline_catch_monthly_soilm4 = accumarray(j,baseline_catch_soilm4,[],@nanmean);
    baseline_catch_monthly_soilm_total = accumarray(j,baseline_catch_soilm_total,[],@nanmean);
    
    %% mod param
    modified_param_data_catch = load([data_dir,sprintf('modified_param_LSM_outputs_Catch_%d.mat',c)]);
    modified_param_data_catch = modified_param_data_catch.Catchment_outputs;
    modified_param_catch_albedo = modified_param_data_catch.albedo;
    modified_param_catch_ecan = modified_param_data_catch.ecan;
    modified_param_catch_etran = modified_param_data_catch.etran;
    modified_param_catch_edir = modified_param_data_catch.edir;
    modified_param_catch_et = modified_param_data_catch.et;
    modified_param_catch_swe = modified_param_data_catch.swe;
    modified_param_catch_lai = modified_param_data_catch.lai;
    modified_param_catch_ugdrnoff = modified_param_data_catch.ugdrnoff;
    modified_param_catch_soilm1 = modified_param_data_catch.soilm1;
    modified_param_catch_soilm2 = modified_param_data_catch.soilm2;
    modified_param_catch_soilm3 = modified_param_data_catch.soilm3;
    modified_param_catch_soilm4 = modified_param_data_catch.soilm4;
    modified_param_catch_soilm_total = modified_param_catch_soilm1.*(0.1/2) + modified_param_catch_soilm2.*(0.3/2) + modified_param_catch_soilm3.*(0.6/2) + modified_param_catch_soilm4.*(1.0/2); 
    %convert accumulated variables to daily (mm/day)
    modified_param_catch_ecan = diff(modified_param_catch_ecan);
    modified_param_catch_etran = diff(modified_param_catch_etran);
    modified_param_catch_edir = diff(modified_param_catch_edir);
    modified_param_catch_et = diff(modified_param_catch_et);
    modified_param_catch_ugdrnoff = diff(modified_param_catch_ugdrnoff);
    %remove day 1 from other var:
    modified_param_catch_albedo(1)=[];
    modified_param_catch_swe(1)=[];
    modified_param_catch_lai(1)=[];
    modified_param_catch_soilm1(1)=[];
    modified_param_catch_soilm2(1)=[];
    modified_param_catch_soilm3(1)=[];
    modified_param_catch_soilm4(1)=[];
    modified_param_catch_soilm_total(1)=[];
    
    %get monthly averages:
    modified_param_catch_monthly_albedo = accumarray(j,modified_param_catch_albedo,[],@nanmean);
    modified_param_catch_monthly_ecan = accumarray(j,modified_param_catch_ecan,[],@nanmean);
    modified_param_catch_monthly_etran = accumarray(j,modified_param_catch_etran,[],@nanmean);
    modified_param_catch_monthly_edir = accumarray(j,modified_param_catch_edir,[],@nanmean);
    modified_param_catch_monthly_et = accumarray(j,modified_param_catch_et,[],@nanmean);
    modified_param_catch_monthly_swe = accumarray(j,modified_param_catch_swe,[],@nanmean);
    modified_param_catch_monthly_lai = accumarray(j,modified_param_catch_lai,[],@nanmean);
    modified_param_catch_monthly_ugdrnoff = accumarray(j,modified_param_catch_ugdrnoff,[],@nanmean);
    modified_param_catch_monthly_soilm1 = accumarray(j,modified_param_catch_soilm1,[],@nanmean);
    modified_param_catch_monthly_soilm2 = accumarray(j,modified_param_catch_soilm2,[],@nanmean);
    modified_param_catch_monthly_soilm3 = accumarray(j,modified_param_catch_soilm3,[],@nanmean);
    modified_param_catch_monthly_soilm4 = accumarray(j,modified_param_catch_soilm4,[],@nanmean);
    modified_param_catch_monthly_soilm_total = accumarray(j,modified_param_catch_soilm_total,[],@nanmean);
    
    %% mod param & veg class - GRASS
    modified_param_GVF_data_catch = load([data_dir,sprintf('modified_param_GVF_LSM_outputs_Catch_%d.mat',c)]);
    modified_param_GVF_data_catch = modified_param_GVF_data_catch.Catchment_outputs;
    modified_param_GVF_catch_albedo = modified_param_GVF_data_catch.albedo;
    modified_param_GVF_catch_ecan = modified_param_GVF_data_catch.ecan;
    modified_param_GVF_catch_etran = modified_param_GVF_data_catch.etran;
    modified_param_GVF_catch_edir = modified_param_GVF_data_catch.edir;
    modified_param_GVF_catch_et = modified_param_GVF_data_catch.et;
    modified_param_GVF_catch_swe = modified_param_GVF_data_catch.swe;
    modified_param_GVF_catch_lai = modified_param_GVF_data_catch.lai;
    modified_param_GVF_catch_ugdrnoff = modified_param_GVF_data_catch.ugdrnoff;
    modified_param_GVF_catch_soilm1 = modified_param_GVF_data_catch.soilm1;
    modified_param_GVF_catch_soilm2 = modified_param_GVF_data_catch.soilm2;
    modified_param_GVF_catch_soilm3 = modified_param_GVF_data_catch.soilm3;
    modified_param_GVF_catch_soilm4 = modified_param_GVF_data_catch.soilm4;
    modified_param_GVF_catch_soilm_total = modified_param_GVF_catch_soilm1.*(0.1/2) + modified_param_GVF_catch_soilm2.*(0.3/2) + modified_param_GVF_catch_soilm3.*(0.6/2) + modified_param_GVF_catch_soilm4.*(1.0/2); 
    %convert accumulated variables to daily (mm/day)
    modified_param_GVF_catch_ecan = diff(modified_param_GVF_catch_ecan);
    modified_param_GVF_catch_etran = diff(modified_param_GVF_catch_etran);
    modified_param_GVF_catch_edir = diff(modified_param_GVF_catch_edir);
    modified_param_GVF_catch_et = diff(modified_param_GVF_catch_et);
    modified_param_GVF_catch_ugdrnoff = diff(modified_param_GVF_catch_ugdrnoff);
    %remove day 1 from other var:
    modified_param_GVF_catch_albedo(1)=[];
    modified_param_GVF_catch_swe(1)=[];
    modified_param_GVF_catch_lai(1)=[];
    modified_param_GVF_catch_soilm1(1)=[];
    modified_param_GVF_catch_soilm2(1)=[];
    modified_param_GVF_catch_soilm3(1)=[];
    modified_param_GVF_catch_soilm4(1)=[];
    modified_param_GVF_catch_soilm_total(1)=[];
    
    %get monthly averages:
    modified_param_GVF_catch_monthly_albedo = accumarray(j,modified_param_GVF_catch_albedo,[],@nanmean);
    modified_param_GVF_catch_monthly_ecan = accumarray(j,modified_param_GVF_catch_ecan,[],@nanmean);
    modified_param_GVF_catch_monthly_etran = accumarray(j,modified_param_GVF_catch_etran,[],@nanmean);
    modified_param_GVF_catch_monthly_edir = accumarray(j,modified_param_GVF_catch_edir,[],@nanmean);
    modified_param_GVF_catch_monthly_et = accumarray(j,modified_param_GVF_catch_et,[],@nanmean);
    modified_param_GVF_catch_monthly_swe = accumarray(j,modified_param_GVF_catch_swe,[],@nanmean);
    modified_param_GVF_catch_monthly_lai = accumarray(j,modified_param_GVF_catch_lai,[],@nanmean);
    modified_param_GVF_catch_monthly_ugdrnoff = accumarray(j,modified_param_GVF_catch_ugdrnoff,[],@nanmean);
    modified_param_GVF_catch_monthly_soilm1 = accumarray(j,modified_param_GVF_catch_soilm1,[],@nanmean);
    modified_param_GVF_catch_monthly_soilm2 = accumarray(j,modified_param_GVF_catch_soilm2,[],@nanmean);
    modified_param_GVF_catch_monthly_soilm3 = accumarray(j,modified_param_GVF_catch_soilm3,[],@nanmean);
    modified_param_GVF_catch_monthly_soilm4 = accumarray(j,modified_param_GVF_catch_soilm4,[],@nanmean);
    modified_param_GVF_catch_monthly_soilm_total = accumarray(j,modified_param_GVF_catch_soilm_total,[],@nanmean);
    
    %% mod param & veg class - BARE
    modified_param_GVF_VegClass_data_catch = load([data_dir,sprintf('modified_param_GVF_VegClass_LSM_outputs_Catch_%d.mat',c)]);
    modified_param_GVF_VegClass_data_catch = modified_param_GVF_VegClass_data_catch.Catchment_outputs;
    modified_param_GVF_VegClass_catch_albedo = modified_param_GVF_VegClass_data_catch.albedo;
    modified_param_GVF_VegClass_catch_ecan = modified_param_GVF_VegClass_data_catch.ecan;
    modified_param_GVF_VegClass_catch_etran = modified_param_GVF_VegClass_data_catch.etran;
    modified_param_GVF_VegClass_catch_edir = modified_param_GVF_VegClass_data_catch.edir;
    modified_param_GVF_VegClass_catch_et = modified_param_GVF_VegClass_data_catch.et;
    modified_param_GVF_VegClass_catch_swe = modified_param_GVF_VegClass_data_catch.swe;
    modified_param_GVF_VegClass_catch_lai = modified_param_GVF_VegClass_data_catch.lai;
    modified_param_GVF_VegClass_catch_ugdrnoff = modified_param_GVF_VegClass_data_catch.ugdrnoff;
    modified_param_GVF_VegClass_catch_soilm1 = modified_param_GVF_VegClass_data_catch.soilm1;
    modified_param_GVF_VegClass_catch_soilm2 = modified_param_GVF_VegClass_data_catch.soilm2;
    modified_param_GVF_VegClass_catch_soilm3 = modified_param_GVF_VegClass_data_catch.soilm3;
    modified_param_GVF_VegClass_catch_soilm4 = modified_param_GVF_VegClass_data_catch.soilm4;
    modified_param_GVF_VegClass_catch_soilm_total = modified_param_GVF_VegClass_catch_soilm1.*(0.1/2) + modified_param_GVF_VegClass_catch_soilm2.*(0.3/2) + modified_param_GVF_VegClass_catch_soilm3.*(0.6/2) + modified_param_GVF_VegClass_catch_soilm4.*(1.0/2); 
    %convert accumulated variables to daily (mm/day)
    modified_param_GVF_VegClass_catch_ecan = diff(modified_param_GVF_VegClass_catch_ecan);
    modified_param_GVF_VegClass_catch_etran = diff(modified_param_GVF_VegClass_catch_etran);
    modified_param_GVF_VegClass_catch_edir = diff(modified_param_GVF_VegClass_catch_edir);
    modified_param_GVF_VegClass_catch_et = diff(modified_param_GVF_VegClass_catch_et);
    modified_param_GVF_VegClass_catch_ugdrnoff = diff(modified_param_GVF_VegClass_catch_ugdrnoff);
    %remove day 1 from other var:
    modified_param_GVF_VegClass_catch_albedo(1)=[];
    modified_param_GVF_VegClass_catch_swe(1)=[];
    modified_param_GVF_VegClass_catch_lai(1)=[];
    modified_param_GVF_VegClass_catch_soilm1(1)=[];
    modified_param_GVF_VegClass_catch_soilm2(1)=[];
    modified_param_GVF_VegClass_catch_soilm3(1)=[];
    modified_param_GVF_VegClass_catch_soilm4(1)=[];
    modified_param_GVF_VegClass_catch_soilm_total(1)=[];
    
    %get monthly averages:
    modified_param_GVF_VegClass_catch_monthly_albedo = accumarray(j,modified_param_GVF_VegClass_catch_albedo,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_ecan = accumarray(j,modified_param_GVF_VegClass_catch_ecan,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_etran = accumarray(j,modified_param_GVF_VegClass_catch_etran,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_edir = accumarray(j,modified_param_GVF_VegClass_catch_edir,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_et = accumarray(j,modified_param_GVF_VegClass_catch_et,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_swe = accumarray(j,modified_param_GVF_VegClass_catch_swe,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_lai = accumarray(j,modified_param_GVF_VegClass_catch_lai,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_ugdrnoff = accumarray(j,modified_param_GVF_VegClass_catch_ugdrnoff,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_soilm1 = accumarray(j,modified_param_GVF_VegClass_catch_soilm1,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_soilm2 = accumarray(j,modified_param_GVF_VegClass_catch_soilm2,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_soilm3 = accumarray(j,modified_param_GVF_VegClass_catch_soilm3,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_soilm4 = accumarray(j,modified_param_GVF_VegClass_catch_soilm4,[],@nanmean);
    modified_param_GVF_VegClass_catch_monthly_soilm_total = accumarray(j,modified_param_GVF_VegClass_catch_soilm_total,[],@nanmean);
    
    
    %% realistic adj
    realistic_data_catch = load([data_dir,sprintf('realistic_LSM_outputs_Catch_%d.mat',c)]);
    realistic_data_catch = realistic_data_catch.Catchment_outputs;
    realistic_catch_albedo = realistic_data_catch.albedo;
    realistic_catch_ecan = realistic_data_catch.ecan;
    realistic_catch_etran = realistic_data_catch.etran;
    realistic_catch_edir = realistic_data_catch.edir;
    realistic_catch_et = realistic_data_catch.et;
    realistic_catch_swe = realistic_data_catch.swe;
    realistic_catch_lai = realistic_data_catch.lai;
    realistic_catch_ugdrnoff = realistic_data_catch.ugdrnoff;
    realistic_catch_soilm1 = realistic_data_catch.soilm1;
    realistic_catch_soilm2 = realistic_data_catch.soilm2;
    realistic_catch_soilm3 = realistic_data_catch.soilm3;
    realistic_catch_soilm4 = realistic_data_catch.soilm4;
    realistic_catch_soilm_total = realistic_catch_soilm1.*(0.1/2) + realistic_catch_soilm2.*(0.3/2) + realistic_catch_soilm3.*(0.6/2) + realistic_catch_soilm4.*(1.0/2); 
    %convert accumulated variables to daily (mm/day)
    realistic_catch_ecan = diff(realistic_catch_ecan);
    realistic_catch_etran = diff(realistic_catch_etran);
    realistic_catch_edir = diff(realistic_catch_edir);
    realistic_catch_et = diff(realistic_catch_et);
    realistic_catch_ugdrnoff = diff(realistic_catch_ugdrnoff);
    %remove day 1 from other var:
    realistic_catch_albedo(1)=[];
    realistic_catch_swe(1)=[];
    realistic_catch_lai(1)=[];
    realistic_catch_soilm1(1)=[];
    realistic_catch_soilm2(1)=[];
    realistic_catch_soilm3(1)=[];
    realistic_catch_soilm4(1)=[];
    realistic_catch_soilm_total(1)=[];
    
    %get monthly averages:
    realistic_catch_monthly_albedo = accumarray(j,realistic_catch_albedo,[],@nanmean);
    realistic_catch_monthly_ecan = accumarray(j,realistic_catch_ecan,[],@nanmean);
    realistic_catch_monthly_etran = accumarray(j,realistic_catch_etran,[],@nanmean);
    realistic_catch_monthly_edir = accumarray(j,realistic_catch_edir,[],@nanmean);
    realistic_catch_monthly_et = accumarray(j,realistic_catch_et,[],@nanmean);
    realistic_catch_monthly_swe = accumarray(j,realistic_catch_swe,[],@nanmean);
    realistic_catch_monthly_lai = accumarray(j,realistic_catch_lai,[],@nanmean);
    realistic_catch_monthly_ugdrnoff = accumarray(j,realistic_catch_ugdrnoff,[],@nanmean);
    realistic_catch_monthly_soilm1 = accumarray(j,realistic_catch_soilm1,[],@nanmean);
    realistic_catch_monthly_soilm2 = accumarray(j,realistic_catch_soilm2,[],@nanmean);
    realistic_catch_monthly_soilm3 = accumarray(j,realistic_catch_soilm3,[],@nanmean);
    realistic_catch_monthly_soilm4 = accumarray(j,realistic_catch_soilm4,[],@nanmean);
    realistic_catch_monthly_soilm_total = accumarray(j,realistic_catch_soilm_total,[],@nanmean);
    
    
    %% plot data
    %plot comparisson of ET and its paritions:
    f=figure;
    f.Position = [-1919        -171        1920         976];
    subplot(2,2,1)
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_et,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_et,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_et,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_et,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_et,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('ET (mm/day)','fontsize',24)
    ylim([0 4])
    xlim([min(monthly_dates) max(monthly_dates)])
    legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
    
    subplot(2,2,2)
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_etran,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_etran,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_etran,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_etran,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_etran,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('E_{tran} (mm/day)','fontsize',24)
    ylim([0 4])
    xlim([min(monthly_dates) max(monthly_dates)])
    
    subplot(2,2,3)
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_edir,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_edir,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_edir,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_edir,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_edir,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('E_{soil} (mm/day)','fontsize',24)
    ylim([0 4])
    xlim([min(monthly_dates) max(monthly_dates)])
    
    subplot(2,2,4)
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_ecan,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_ecan,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_ecan,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_ecan,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_ecan,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('E_{can} (mm/day)','fontsize',24)
    ylim([0 4])
    xlim([min(monthly_dates) max(monthly_dates)])
    outfilename = sprintf('Catch%d_ET_comparisson.png',c);
    saveas(f,[outdir,outfilename])
    
    %plot comparisson of SM total and its layers:
    f=figure;
    f.Position = [-1919        -171        1920         976];
    subplot(4,2,[1:4])
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_soilm_total,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_soilm_total,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm_total,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm_total,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_soilm_total,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('total column SM (mm/mm)','fontsize',24)
    ylim([0.1 0.5])
    xlim([min(monthly_dates) max(monthly_dates)])
    legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
    
    subplot(4,2,5)
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_soilm1,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_soilm1,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm1,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm1,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_soilm1,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel({'soilm1'; '(0-10cm)'},'fontsize',24)
    ylim([0.1 0.5])
    xlim([min(monthly_dates) max(monthly_dates)])
    
    subplot(4,2,6)
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_soilm2,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_soilm2,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm2,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm2,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_soilm2,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel({'soilm2'; '(10-40cm)'},'fontsize',24)
    ylim([0.1 0.5])
    xlim([min(monthly_dates) max(monthly_dates)])
    
    subplot(4,2,7)
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_soilm3,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_soilm3,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm3,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm3,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_soilm3,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel({'soilm3'; '(40-100cm)'},'fontsize',24)
    ylim([0.1 0.5])
    xlim([min(monthly_dates) max(monthly_dates)])
    
    subplot(4,2,8)
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_soilm4,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_soilm4,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm4,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm4,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_soilm4,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel({'soilm4'; '(100-200cm)'},'fontsize',24)
    ylim([0.05 0.5])
    xlim([min(monthly_dates) max(monthly_dates)])
    outfilename = sprintf('Catch%d_SOILM_comparisson.png',c);
    saveas(f,[outdir,outfilename])
    
    %plot comparisson of SWE
    f=figure;
    f.Position = [-1919         418        1920         387];
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_swe,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_swe,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_swe,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_swe,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_swe,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('monthly SWE (mm)','fontsize',24)
    xlim([min(monthly_dates) max(monthly_dates)])
    legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
    outfilename = sprintf('Catch%d_SWE_comparisson.png',c);
    saveas(f,[outdir,outfilename])
    
    %plot comparisson of lai
    f=figure;
    f.Position = [-1919         418        1920         387];
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_lai,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_lai,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_lai,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_lai,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_lai,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('monthly LAI (mm/mm)','fontsize',24)
    xlim([min(monthly_dates) max(monthly_dates)])
    legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
    outfilename = sprintf('Catch%d_lai_comparisson.png',c);
    saveas(f,[outdir,outfilename])
    
    %plot comparisson of ugdrnoff
    f=figure;
    f.Position = [-1919         418        1920         387];
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_ugdrnoff,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_ugdrnoff,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_ugdrnoff,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_ugdrnoff,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_ugdrnoff,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('monthly baseflow (mm/day)','fontsize',24)
    xlim([min(monthly_dates) max(monthly_dates)])
    legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
    outfilename = sprintf('Catch%d_baseflow_comparisson.png',c);
    saveas(f,[outdir,outfilename])
    
    %plot comparisson of albedo
    f=figure;
    f.Position = [-1919         418        1920         387];
    hold on
    p1 = plot(monthly_dates,baseline_catch_monthly_albedo,'-','linewidth',2,'color','k');
    p2 = plot(monthly_dates, modified_param_catch_monthly_albedo,'-','linewidth',2,'color','b');
    p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_albedo,'-','linewidth',2,'color',[0 0.5 0]);
    p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_albedo,'--','linewidth',2,'color',[0 0.5 0]);
    p5 = plot(monthly_dates, realistic_catch_monthly_albedo,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    grid on
    box on
    set(gca,'fontsize',24)
    datetick('x','mm/yyyy','keepticks','keeplimits')
    ylabel('Albedo','fontsize',24)
    xlim([min(monthly_dates) max(monthly_dates)])
    legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
    outfilename = sprintf('Catch%d_albedo_comparisson.png',c);
    saveas(f,[outdir,outfilename])
end

% % %% this 3rd part of the code compares the time series data created in part 1 (for burn area averages)
% % outdir = '/Users/abolafia/ASO_Fire/Plots/PaperPlots/Catchment_WB_comparissons/';
% % if exist(outdir,'dir') == 0
% %     system(['mkdir -p ',outdir])
% % end
% % data_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs/';
% % catch_names = {'East_Branch_NF_Feather','Mid_Fork_Feather','NF_Feather'};
% % for c=1:3
% %     
% %     %% baseline data
% %     baseline_data_catch = load([data_dir,sprintf('baseline_LSM_outputs_BurnedPixels_Catch_%d.mat',c)]);
% %     
% %     baseline_data_catch = baseline_data_catch.Catchment_outputs;
% %     dates = double(baseline_data_catch.date);
% %     datevecs = datevec(dates);
% %     baseline_catch_albedo = baseline_data_catch.albedo;
% %     baseline_catch_ecan = baseline_data_catch.ecan;
% %     baseline_catch_etran = baseline_data_catch.etran;
% %     baseline_catch_edir = baseline_data_catch.edir;
% %     baseline_catch_et = baseline_data_catch.et;
% %     baseline_catch_swe = baseline_data_catch.swe;
% %     baseline_catch_lai = baseline_data_catch.lai;
% %     baseline_catch_ugdrnoff = baseline_data_catch.ugdrnoff;
% %     baseline_catch_soilm1 = baseline_data_catch.soilm1;
% %     baseline_catch_soilm2 = baseline_data_catch.soilm2;
% %     baseline_catch_soilm3 = baseline_data_catch.soilm3;
% %     baseline_catch_soilm4 = baseline_data_catch.soilm4;
% %     baseline_catch_soilm_total = baseline_catch_soilm1.*(0.1/2) + baseline_catch_soilm2.*(0.3/2) + baseline_catch_soilm3.*(0.6/2) + baseline_catch_soilm4.*(1.0/2); 
% %     
% %     %convert accumulated variables to daily (mm/day)
% %     baseline_catch_ecan = diff(baseline_catch_ecan);
% %     baseline_catch_etran = diff(baseline_catch_etran);
% %     baseline_catch_edir = diff(baseline_catch_edir);
% %     baseline_catch_et = diff(baseline_catch_et);
% %     baseline_catch_ugdrnoff = diff(baseline_catch_ugdrnoff);
% %     %remove day 1 from other var:
% %     dates(1) = [];
% %     datevecs(1,:)=[];
% %     baseline_catch_albedo(1)=[];
% %     baseline_catch_swe(1)=[];
% %     baseline_catch_lai(1)=[];
% %     baseline_catch_soilm1(1)=[];
% %     baseline_catch_soilm2(1)=[];
% %     baseline_catch_soilm3(1)=[];
% %     baseline_catch_soilm4(1)=[];
% %     baseline_catch_soilm_total(1)=[];
% % 
% %     %aggregate to monthly average:
% %     [u,~,j] = unique(datevecs(:,1:2),'rows','stable');
% %     monthly_dates = datenum([u,ones(length(u),1)]);
% %     baseline_catch_monthly_albedo = accumarray(j,baseline_catch_albedo,[],@nanmean);
% %     baseline_catch_monthly_ecan = accumarray(j,baseline_catch_ecan,[],@nanmean);
% %     baseline_catch_monthly_etran = accumarray(j,baseline_catch_etran,[],@nanmean);
% %     baseline_catch_monthly_edir = accumarray(j,baseline_catch_edir,[],@nanmean);
% %     baseline_catch_monthly_et = accumarray(j,baseline_catch_et,[],@nanmean);
% %     baseline_catch_monthly_swe = accumarray(j,baseline_catch_swe,[],@nanmean);
% %     baseline_catch_monthly_lai = accumarray(j,baseline_catch_lai,[],@nanmean);
% %     baseline_catch_monthly_ugdrnoff = accumarray(j,baseline_catch_ugdrnoff,[],@nanmean);
% %     baseline_catch_monthly_soilm1 = accumarray(j,baseline_catch_soilm1,[],@nanmean);
% %     baseline_catch_monthly_soilm2 = accumarray(j,baseline_catch_soilm2,[],@nanmean);
% %     baseline_catch_monthly_soilm3 = accumarray(j,baseline_catch_soilm3,[],@nanmean);
% %     baseline_catch_monthly_soilm4 = accumarray(j,baseline_catch_soilm4,[],@nanmean);
% %     baseline_catch_monthly_soilm_total = accumarray(j,baseline_catch_soilm_total,[],@nanmean);
% %     
% %     %% mod param
% %     modified_param_data_catch = load([data_dir,sprintf('modified_param_LSM_outputs_BurnedPixels_Catch_%d.mat',c)]);
% %     modified_param_data_catch = modified_param_data_catch.Catchment_outputs;
% %     modified_param_catch_albedo = modified_param_data_catch.albedo;
% %     modified_param_catch_ecan = modified_param_data_catch.ecan;
% %     modified_param_catch_etran = modified_param_data_catch.etran;
% %     modified_param_catch_edir = modified_param_data_catch.edir;
% %     modified_param_catch_et = modified_param_data_catch.et;
% %     modified_param_catch_swe = modified_param_data_catch.swe;
% %     modified_param_catch_lai = modified_param_data_catch.lai;
% %     modified_param_catch_ugdrnoff = modified_param_data_catch.ugdrnoff;
% %     modified_param_catch_soilm1 = modified_param_data_catch.soilm1;
% %     modified_param_catch_soilm2 = modified_param_data_catch.soilm2;
% %     modified_param_catch_soilm3 = modified_param_data_catch.soilm3;
% %     modified_param_catch_soilm4 = modified_param_data_catch.soilm4;
% %     modified_param_catch_soilm_total = modified_param_catch_soilm1.*(0.1/2) + modified_param_catch_soilm2.*(0.3/2) + modified_param_catch_soilm3.*(0.6/2) + modified_param_catch_soilm4.*(1.0/2); 
% %     %convert accumulated variables to daily (mm/day)
% %     modified_param_catch_ecan = diff(modified_param_catch_ecan);
% %     modified_param_catch_etran = diff(modified_param_catch_etran);
% %     modified_param_catch_edir = diff(modified_param_catch_edir);
% %     modified_param_catch_et = diff(modified_param_catch_et);
% %     modified_param_catch_ugdrnoff = diff(modified_param_catch_ugdrnoff);
% %     %remove day 1 from other var:
% %     modified_param_catch_albedo(1)=[];
% %     modified_param_catch_swe(1)=[];
% %     modified_param_catch_lai(1)=[];
% %     modified_param_catch_soilm1(1)=[];
% %     modified_param_catch_soilm2(1)=[];
% %     modified_param_catch_soilm3(1)=[];
% %     modified_param_catch_soilm4(1)=[];
% %     modified_param_catch_soilm_total(1)=[];
% %     
% %     %get monthly averages:
% %     modified_param_catch_monthly_albedo = accumarray(j,modified_param_catch_albedo,[],@nanmean);
% %     modified_param_catch_monthly_ecan = accumarray(j,modified_param_catch_ecan,[],@nanmean);
% %     modified_param_catch_monthly_etran = accumarray(j,modified_param_catch_etran,[],@nanmean);
% %     modified_param_catch_monthly_edir = accumarray(j,modified_param_catch_edir,[],@nanmean);
% %     modified_param_catch_monthly_et = accumarray(j,modified_param_catch_et,[],@nanmean);
% %     modified_param_catch_monthly_swe = accumarray(j,modified_param_catch_swe,[],@nanmean);
% %     modified_param_catch_monthly_lai = accumarray(j,modified_param_catch_lai,[],@nanmean);
% %     modified_param_catch_monthly_ugdrnoff = accumarray(j,modified_param_catch_ugdrnoff,[],@nanmean);
% %     modified_param_catch_monthly_soilm1 = accumarray(j,modified_param_catch_soilm1,[],@nanmean);
% %     modified_param_catch_monthly_soilm2 = accumarray(j,modified_param_catch_soilm2,[],@nanmean);
% %     modified_param_catch_monthly_soilm3 = accumarray(j,modified_param_catch_soilm3,[],@nanmean);
% %     modified_param_catch_monthly_soilm4 = accumarray(j,modified_param_catch_soilm4,[],@nanmean);
% %     modified_param_catch_monthly_soilm_total = accumarray(j,modified_param_catch_soilm_total,[],@nanmean);
% %     
% %     %% mod param & veg class - GRASS
% %     modified_param_GVF_data_catch = load([data_dir,sprintf('modified_param_GVF_LSM_outputs_BurnedPixels_Catch_%d.mat',c)]);
% %     modified_param_GVF_data_catch = modified_param_GVF_data_catch.Catchment_outputs;
% %     modified_param_GVF_catch_albedo = modified_param_GVF_data_catch.albedo;
% %     modified_param_GVF_catch_ecan = modified_param_GVF_data_catch.ecan;
% %     modified_param_GVF_catch_etran = modified_param_GVF_data_catch.etran;
% %     modified_param_GVF_catch_edir = modified_param_GVF_data_catch.edir;
% %     modified_param_GVF_catch_et = modified_param_GVF_data_catch.et;
% %     modified_param_GVF_catch_swe = modified_param_GVF_data_catch.swe;
% %     modified_param_GVF_catch_lai = modified_param_GVF_data_catch.lai;
% %     modified_param_GVF_catch_ugdrnoff = modified_param_GVF_data_catch.ugdrnoff;
% %     modified_param_GVF_catch_soilm1 = modified_param_GVF_data_catch.soilm1;
% %     modified_param_GVF_catch_soilm2 = modified_param_GVF_data_catch.soilm2;
% %     modified_param_GVF_catch_soilm3 = modified_param_GVF_data_catch.soilm3;
% %     modified_param_GVF_catch_soilm4 = modified_param_GVF_data_catch.soilm4;
% %     modified_param_GVF_catch_soilm_total = modified_param_GVF_catch_soilm1.*(0.1/2) + modified_param_GVF_catch_soilm2.*(0.3/2) + modified_param_GVF_catch_soilm3.*(0.6/2) + modified_param_GVF_catch_soilm4.*(1.0/2); 
% %     %convert accumulated variables to daily (mm/day)
% %     modified_param_GVF_catch_ecan = diff(modified_param_GVF_catch_ecan);
% %     modified_param_GVF_catch_etran = diff(modified_param_GVF_catch_etran);
% %     modified_param_GVF_catch_edir = diff(modified_param_GVF_catch_edir);
% %     modified_param_GVF_catch_et = diff(modified_param_GVF_catch_et);
% %     modified_param_GVF_catch_ugdrnoff = diff(modified_param_GVF_catch_ugdrnoff);
% %     %remove day 1 from other var:
% %     modified_param_GVF_catch_albedo(1)=[];
% %     modified_param_GVF_catch_swe(1)=[];
% %     modified_param_GVF_catch_lai(1)=[];
% %     modified_param_GVF_catch_soilm1(1)=[];
% %     modified_param_GVF_catch_soilm2(1)=[];
% %     modified_param_GVF_catch_soilm3(1)=[];
% %     modified_param_GVF_catch_soilm4(1)=[];
% %     modified_param_GVF_catch_soilm_total(1)=[];
% %     
% %     %get monthly averages:
% %     modified_param_GVF_catch_monthly_albedo = accumarray(j,modified_param_GVF_catch_albedo,[],@nanmean);
% %     modified_param_GVF_catch_monthly_ecan = accumarray(j,modified_param_GVF_catch_ecan,[],@nanmean);
% %     modified_param_GVF_catch_monthly_etran = accumarray(j,modified_param_GVF_catch_etran,[],@nanmean);
% %     modified_param_GVF_catch_monthly_edir = accumarray(j,modified_param_GVF_catch_edir,[],@nanmean);
% %     modified_param_GVF_catch_monthly_et = accumarray(j,modified_param_GVF_catch_et,[],@nanmean);
% %     modified_param_GVF_catch_monthly_swe = accumarray(j,modified_param_GVF_catch_swe,[],@nanmean);
% %     modified_param_GVF_catch_monthly_lai = accumarray(j,modified_param_GVF_catch_lai,[],@nanmean);
% %     modified_param_GVF_catch_monthly_ugdrnoff = accumarray(j,modified_param_GVF_catch_ugdrnoff,[],@nanmean);
% %     modified_param_GVF_catch_monthly_soilm1 = accumarray(j,modified_param_GVF_catch_soilm1,[],@nanmean);
% %     modified_param_GVF_catch_monthly_soilm2 = accumarray(j,modified_param_GVF_catch_soilm2,[],@nanmean);
% %     modified_param_GVF_catch_monthly_soilm3 = accumarray(j,modified_param_GVF_catch_soilm3,[],@nanmean);
% %     modified_param_GVF_catch_monthly_soilm4 = accumarray(j,modified_param_GVF_catch_soilm4,[],@nanmean);
% %     modified_param_GVF_catch_monthly_soilm_total = accumarray(j,modified_param_GVF_catch_soilm_total,[],@nanmean);
% %     
% %     %% mod param & veg class - BARE
% %     modified_param_GVF_VegClass_data_catch = load([data_dir,sprintf('modified_param_GVF_VegClass_LSM_outputs_BurnedPixels_Catch_%d.mat',c)]);
% %     modified_param_GVF_VegClass_data_catch = modified_param_GVF_VegClass_data_catch.Catchment_outputs;
% %     modified_param_GVF_VegClass_catch_albedo = modified_param_GVF_VegClass_data_catch.albedo;
% %     modified_param_GVF_VegClass_catch_ecan = modified_param_GVF_VegClass_data_catch.ecan;
% %     modified_param_GVF_VegClass_catch_etran = modified_param_GVF_VegClass_data_catch.etran;
% %     modified_param_GVF_VegClass_catch_edir = modified_param_GVF_VegClass_data_catch.edir;
% %     modified_param_GVF_VegClass_catch_et = modified_param_GVF_VegClass_data_catch.et;
% %     modified_param_GVF_VegClass_catch_swe = modified_param_GVF_VegClass_data_catch.swe;
% %     modified_param_GVF_VegClass_catch_lai = modified_param_GVF_VegClass_data_catch.lai;
% %     modified_param_GVF_VegClass_catch_ugdrnoff = modified_param_GVF_VegClass_data_catch.ugdrnoff;
% %     modified_param_GVF_VegClass_catch_soilm1 = modified_param_GVF_VegClass_data_catch.soilm1;
% %     modified_param_GVF_VegClass_catch_soilm2 = modified_param_GVF_VegClass_data_catch.soilm2;
% %     modified_param_GVF_VegClass_catch_soilm3 = modified_param_GVF_VegClass_data_catch.soilm3;
% %     modified_param_GVF_VegClass_catch_soilm4 = modified_param_GVF_VegClass_data_catch.soilm4;
% %     modified_param_GVF_VegClass_catch_soilm_total = modified_param_GVF_VegClass_catch_soilm1.*(0.1/2) + modified_param_GVF_VegClass_catch_soilm2.*(0.3/2) + modified_param_GVF_VegClass_catch_soilm3.*(0.6/2) + modified_param_GVF_VegClass_catch_soilm4.*(1.0/2); 
% %     %convert accumulated variables to daily (mm/day)
% %     modified_param_GVF_VegClass_catch_ecan = diff(modified_param_GVF_VegClass_catch_ecan);
% %     modified_param_GVF_VegClass_catch_etran = diff(modified_param_GVF_VegClass_catch_etran);
% %     modified_param_GVF_VegClass_catch_edir = diff(modified_param_GVF_VegClass_catch_edir);
% %     modified_param_GVF_VegClass_catch_et = diff(modified_param_GVF_VegClass_catch_et);
% %     modified_param_GVF_VegClass_catch_ugdrnoff = diff(modified_param_GVF_VegClass_catch_ugdrnoff);
% %     %remove day 1 from other var:
% %     modified_param_GVF_VegClass_catch_albedo(1)=[];
% %     modified_param_GVF_VegClass_catch_swe(1)=[];
% %     modified_param_GVF_VegClass_catch_lai(1)=[];
% %     modified_param_GVF_VegClass_catch_soilm1(1)=[];
% %     modified_param_GVF_VegClass_catch_soilm2(1)=[];
% %     modified_param_GVF_VegClass_catch_soilm3(1)=[];
% %     modified_param_GVF_VegClass_catch_soilm4(1)=[];
% %     modified_param_GVF_VegClass_catch_soilm_total(1)=[];
% %     
% %     %get monthly averages:
% %     modified_param_GVF_VegClass_catch_monthly_albedo = accumarray(j,modified_param_GVF_VegClass_catch_albedo,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_ecan = accumarray(j,modified_param_GVF_VegClass_catch_ecan,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_etran = accumarray(j,modified_param_GVF_VegClass_catch_etran,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_edir = accumarray(j,modified_param_GVF_VegClass_catch_edir,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_et = accumarray(j,modified_param_GVF_VegClass_catch_et,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_swe = accumarray(j,modified_param_GVF_VegClass_catch_swe,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_lai = accumarray(j,modified_param_GVF_VegClass_catch_lai,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_ugdrnoff = accumarray(j,modified_param_GVF_VegClass_catch_ugdrnoff,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_soilm1 = accumarray(j,modified_param_GVF_VegClass_catch_soilm1,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_soilm2 = accumarray(j,modified_param_GVF_VegClass_catch_soilm2,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_soilm3 = accumarray(j,modified_param_GVF_VegClass_catch_soilm3,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_soilm4 = accumarray(j,modified_param_GVF_VegClass_catch_soilm4,[],@nanmean);
% %     modified_param_GVF_VegClass_catch_monthly_soilm_total = accumarray(j,modified_param_GVF_VegClass_catch_soilm_total,[],@nanmean);
% %     
% %     %% mod param & veg class & albedo - GRASS
% %     modified_param_albedo_GRASS_data_catch = load([data_dir,sprintf('modified_param_albedo_GRASS_LSM_outputs_BurnedPixels_Catch_%d.mat',c)]);
% %     modified_param_albedo_GRASS_data_catch = modified_param_albedo_GRASS_data_catch.Catchment_outputs;
% %     modified_param_albedo_GRASS_catch_albedo = modified_param_albedo_GRASS_data_catch.albedo;
% %     modified_param_albedo_GRASS_catch_ecan = modified_param_albedo_GRASS_data_catch.ecan;
% %     modified_param_albedo_GRASS_catch_etran = modified_param_albedo_GRASS_data_catch.etran;
% %     modified_param_albedo_GRASS_catch_edir = modified_param_albedo_GRASS_data_catch.edir;
% %     modified_param_albedo_GRASS_catch_et = modified_param_albedo_GRASS_data_catch.et;
% %     modified_param_albedo_GRASS_catch_swe = modified_param_albedo_GRASS_data_catch.swe;
% %     modified_param_albedo_GRASS_catch_lai = modified_param_albedo_GRASS_data_catch.lai;
% %     modified_param_albedo_GRASS_catch_ugdrnoff = modified_param_albedo_GRASS_data_catch.ugdrnoff;
% %     modified_param_albedo_GRASS_catch_soilm1 = modified_param_albedo_GRASS_data_catch.soilm1;
% %     modified_param_albedo_GRASS_catch_soilm2 = modified_param_albedo_GRASS_data_catch.soilm2;
% %     modified_param_albedo_GRASS_catch_soilm3 = modified_param_albedo_GRASS_data_catch.soilm3;
% %     modified_param_albedo_GRASS_catch_soilm4 = modified_param_albedo_GRASS_data_catch.soilm4;
% %     modified_param_albedo_GRASS_catch_soilm_total = modified_param_albedo_GRASS_catch_soilm1.*(0.1/2) + modified_param_albedo_GRASS_catch_soilm2.*(0.3/2) + modified_param_albedo_GRASS_catch_soilm3.*(0.6/2) + modified_param_albedo_GRASS_catch_soilm4.*(1.0/2); 
% %     %convert accumulated variables to daily (mm/day)
% %     modified_param_albedo_GRASS_catch_ecan = diff(modified_param_albedo_GRASS_catch_ecan);
% %     modified_param_albedo_GRASS_catch_etran = diff(modified_param_albedo_GRASS_catch_etran);
% %     modified_param_albedo_GRASS_catch_edir = diff(modified_param_albedo_GRASS_catch_edir);
% %     modified_param_albedo_GRASS_catch_et = diff(modified_param_albedo_GRASS_catch_et);
% %     modified_param_albedo_GRASS_catch_ugdrnoff = diff(modified_param_albedo_GRASS_catch_ugdrnoff);
% %     %remove day 1 from other var:
% %     modified_param_albedo_GRASS_catch_albedo(1)=[];
% %     modified_param_albedo_GRASS_catch_swe(1)=[];
% %     modified_param_albedo_GRASS_catch_lai(1)=[];
% %     modified_param_albedo_GRASS_catch_soilm1(1)=[];
% %     modified_param_albedo_GRASS_catch_soilm2(1)=[];
% %     modified_param_albedo_GRASS_catch_soilm3(1)=[];
% %     modified_param_albedo_GRASS_catch_soilm4(1)=[];
% %     modified_param_albedo_GRASS_catch_soilm_total(1)=[];
% %     
% %     %get monthly averages:
% %     modified_param_albedo_GRASS_catch_monthly_albedo = accumarray(j,modified_param_albedo_GRASS_catch_albedo,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_ecan = accumarray(j,modified_param_albedo_GRASS_catch_ecan,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_etran = accumarray(j,modified_param_albedo_GRASS_catch_etran,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_edir = accumarray(j,modified_param_albedo_GRASS_catch_edir,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_et = accumarray(j,modified_param_albedo_GRASS_catch_et,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_swe = accumarray(j,modified_param_albedo_GRASS_catch_swe,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_lai = accumarray(j,modified_param_albedo_GRASS_catch_lai,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_ugdrnoff = accumarray(j,modified_param_albedo_GRASS_catch_ugdrnoff,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_soilm1 = accumarray(j,modified_param_albedo_GRASS_catch_soilm1,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_soilm2 = accumarray(j,modified_param_albedo_GRASS_catch_soilm2,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_soilm3 = accumarray(j,modified_param_albedo_GRASS_catch_soilm3,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_soilm4 = accumarray(j,modified_param_albedo_GRASS_catch_soilm4,[],@nanmean);
% %     modified_param_albedo_GRASS_catch_monthly_soilm_total = accumarray(j,modified_param_albedo_GRASS_catch_soilm_total,[],@nanmean);
% %     
% %     %% mod param & veg class & albedo - BARE
% %     modified_param_albedo_BARE_data_catch = load([data_dir,sprintf('modified_param_albedo_BARE_LSM_outputs_BurnedPixels_Catch_%d.mat',c)]);
% %     modified_param_albedo_BARE_data_catch = modified_param_albedo_BARE_data_catch.Catchment_outputs;
% %     modified_param_albedo_BARE_catch_albedo = modified_param_albedo_BARE_data_catch.albedo;
% %     modified_param_albedo_BARE_catch_ecan = modified_param_albedo_BARE_data_catch.ecan;
% %     modified_param_albedo_BARE_catch_etran = modified_param_albedo_BARE_data_catch.etran;
% %     modified_param_albedo_BARE_catch_edir = modified_param_albedo_BARE_data_catch.edir;
% %     modified_param_albedo_BARE_catch_et = modified_param_albedo_BARE_data_catch.et;
% %     modified_param_albedo_BARE_catch_swe = modified_param_albedo_BARE_data_catch.swe;
% %     modified_param_albedo_BARE_catch_lai = modified_param_albedo_BARE_data_catch.lai;
% %     modified_param_albedo_BARE_catch_ugdrnoff = modified_param_albedo_BARE_data_catch.ugdrnoff;
% %     modified_param_albedo_BARE_catch_soilm1 = modified_param_albedo_BARE_data_catch.soilm1;
% %     modified_param_albedo_BARE_catch_soilm2 = modified_param_albedo_BARE_data_catch.soilm2;
% %     modified_param_albedo_BARE_catch_soilm3 = modified_param_albedo_BARE_data_catch.soilm3;
% %     modified_param_albedo_BARE_catch_soilm4 = modified_param_albedo_BARE_data_catch.soilm4;
% %     modified_param_albedo_BARE_catch_soilm_total = modified_param_albedo_BARE_catch_soilm1.*(0.1/2) + modified_param_albedo_BARE_catch_soilm2.*(0.3/2) + modified_param_albedo_BARE_catch_soilm3.*(0.6/2) + modified_param_albedo_BARE_catch_soilm4.*(1.0/2); 
% %     %convert accumulated variables to daily (mm/day)
% %     modified_param_albedo_BARE_catch_ecan = diff(modified_param_albedo_BARE_catch_ecan);
% %     modified_param_albedo_BARE_catch_etran = diff(modified_param_albedo_BARE_catch_etran);
% %     modified_param_albedo_BARE_catch_edir = diff(modified_param_albedo_BARE_catch_edir);
% %     modified_param_albedo_BARE_catch_et = diff(modified_param_albedo_BARE_catch_et);
% %     modified_param_albedo_BARE_catch_ugdrnoff = diff(modified_param_albedo_BARE_catch_ugdrnoff);
% %     %remove day 1 from other var:
% %     modified_param_albedo_BARE_catch_albedo(1)=[];
% %     modified_param_albedo_BARE_catch_swe(1)=[];
% %     modified_param_albedo_BARE_catch_lai(1)=[];
% %     modified_param_albedo_BARE_catch_soilm1(1)=[];
% %     modified_param_albedo_BARE_catch_soilm2(1)=[];
% %     modified_param_albedo_BARE_catch_soilm3(1)=[];
% %     modified_param_albedo_BARE_catch_soilm4(1)=[];
% %     modified_param_albedo_BARE_catch_soilm_total(1)=[];
% %     
% %     %get monthly averages:
% %     modified_param_albedo_BARE_catch_monthly_albedo = accumarray(j,modified_param_albedo_BARE_catch_albedo,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_ecan = accumarray(j,modified_param_albedo_BARE_catch_ecan,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_etran = accumarray(j,modified_param_albedo_BARE_catch_etran,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_edir = accumarray(j,modified_param_albedo_BARE_catch_edir,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_et = accumarray(j,modified_param_albedo_BARE_catch_et,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_swe = accumarray(j,modified_param_albedo_BARE_catch_swe,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_lai = accumarray(j,modified_param_albedo_BARE_catch_lai,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_ugdrnoff = accumarray(j,modified_param_albedo_BARE_catch_ugdrnoff,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_soilm1 = accumarray(j,modified_param_albedo_BARE_catch_soilm1,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_soilm2 = accumarray(j,modified_param_albedo_BARE_catch_soilm2,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_soilm3 = accumarray(j,modified_param_albedo_BARE_catch_soilm3,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_soilm4 = accumarray(j,modified_param_albedo_BARE_catch_soilm4,[],@nanmean);
% %     modified_param_albedo_BARE_catch_monthly_soilm_total = accumarray(j,modified_param_albedo_BARE_catch_soilm_total,[],@nanmean);
% %     
% %     %% realistic adj
% %     realistic_data_catch = load([data_dir,sprintf('realistic_LSM_outputs_BurnedPixels_Catch_%d.mat',c)]);
% %     realistic_data_catch = realistic_data_catch.Catchment_outputs;
% %     realistic_catch_albedo = realistic_data_catch.albedo;
% %     realistic_catch_ecan = realistic_data_catch.ecan;
% %     realistic_catch_etran = realistic_data_catch.etran;
% %     realistic_catch_edir = realistic_data_catch.edir;
% %     realistic_catch_et = realistic_data_catch.et;
% %     realistic_catch_swe = realistic_data_catch.swe;
% %     realistic_catch_lai = realistic_data_catch.lai;
% %     realistic_catch_ugdrnoff = realistic_data_catch.ugdrnoff;
% %     realistic_catch_soilm1 = realistic_data_catch.soilm1;
% %     realistic_catch_soilm2 = realistic_data_catch.soilm2;
% %     realistic_catch_soilm3 = realistic_data_catch.soilm3;
% %     realistic_catch_soilm4 = realistic_data_catch.soilm4;
% %     realistic_catch_soilm_total = realistic_catch_soilm1.*(0.1/2) + realistic_catch_soilm2.*(0.3/2) + realistic_catch_soilm3.*(0.6/2) + realistic_catch_soilm4.*(1.0/2); 
% %     %convert accumulated variables to daily (mm/day)
% %     realistic_catch_ecan = diff(realistic_catch_ecan);
% %     realistic_catch_etran = diff(realistic_catch_etran);
% %     realistic_catch_edir = diff(realistic_catch_edir);
% %     realistic_catch_et = diff(realistic_catch_et);
% %     realistic_catch_ugdrnoff = diff(realistic_catch_ugdrnoff);
% %     %remove day 1 from other var:
% %     realistic_catch_albedo(1)=[];
% %     realistic_catch_swe(1)=[];
% %     realistic_catch_lai(1)=[];
% %     realistic_catch_soilm1(1)=[];
% %     realistic_catch_soilm2(1)=[];
% %     realistic_catch_soilm3(1)=[];
% %     realistic_catch_soilm4(1)=[];
% %     realistic_catch_soilm_total(1)=[];
% %     
% %     %get monthly averages:
% %     realistic_catch_monthly_albedo = accumarray(j,realistic_catch_albedo,[],@nanmean);
% %     realistic_catch_monthly_ecan = accumarray(j,realistic_catch_ecan,[],@nanmean);
% %     realistic_catch_monthly_etran = accumarray(j,realistic_catch_etran,[],@nanmean);
% %     realistic_catch_monthly_edir = accumarray(j,realistic_catch_edir,[],@nanmean);
% %     realistic_catch_monthly_et = accumarray(j,realistic_catch_et,[],@nanmean);
% %     realistic_catch_monthly_swe = accumarray(j,realistic_catch_swe,[],@nanmean);
% %     realistic_catch_monthly_lai = accumarray(j,realistic_catch_lai,[],@nanmean);
% %     realistic_catch_monthly_ugdrnoff = accumarray(j,realistic_catch_ugdrnoff,[],@nanmean);
% %     realistic_catch_monthly_soilm1 = accumarray(j,realistic_catch_soilm1,[],@nanmean);
% %     realistic_catch_monthly_soilm2 = accumarray(j,realistic_catch_soilm2,[],@nanmean);
% %     realistic_catch_monthly_soilm3 = accumarray(j,realistic_catch_soilm3,[],@nanmean);
% %     realistic_catch_monthly_soilm4 = accumarray(j,realistic_catch_soilm4,[],@nanmean);
% %     realistic_catch_monthly_soilm_total = accumarray(j,realistic_catch_soilm_total,[],@nanmean);
% %     
% %     %% realistic_cleaner adj
% %     realistic_cleaner_data_catch = load([data_dir,sprintf('realistic_cleaner_LSM_outputs_BurnedPixels_Catch_%d.mat',c)]);
% %     realistic_cleaner_data_catch = realistic_cleaner_data_catch.Catchment_outputs;
% %     realistic_cleaner_catch_albedo = realistic_cleaner_data_catch.albedo;
% %     realistic_cleaner_catch_ecan = realistic_cleaner_data_catch.ecan;
% %     realistic_cleaner_catch_etran = realistic_cleaner_data_catch.etran;
% %     realistic_cleaner_catch_edir = realistic_cleaner_data_catch.edir;
% %     realistic_cleaner_catch_et = realistic_cleaner_data_catch.et;
% %     realistic_cleaner_catch_swe = realistic_cleaner_data_catch.swe;
% %     realistic_cleaner_catch_lai = realistic_cleaner_data_catch.lai;
% %     realistic_cleaner_catch_ugdrnoff = realistic_cleaner_data_catch.ugdrnoff;
% %     realistic_cleaner_catch_soilm1 = realistic_cleaner_data_catch.soilm1;
% %     realistic_cleaner_catch_soilm2 = realistic_cleaner_data_catch.soilm2;
% %     realistic_cleaner_catch_soilm3 = realistic_cleaner_data_catch.soilm3;
% %     realistic_cleaner_catch_soilm4 = realistic_cleaner_data_catch.soilm4;
% %     realistic_cleaner_catch_soilm_total = realistic_cleaner_catch_soilm1.*(0.1/2) + realistic_cleaner_catch_soilm2.*(0.3/2) + realistic_cleaner_catch_soilm3.*(0.6/2) + realistic_cleaner_catch_soilm4.*(1.0/2); 
% %     %convert accumulated variables to daily (mm/day)
% %     realistic_cleaner_catch_ecan = diff(realistic_cleaner_catch_ecan);
% %     realistic_cleaner_catch_etran = diff(realistic_cleaner_catch_etran);
% %     realistic_cleaner_catch_edir = diff(realistic_cleaner_catch_edir);
% %     realistic_cleaner_catch_et = diff(realistic_cleaner_catch_et);
% %     realistic_cleaner_catch_ugdrnoff = diff(realistic_cleaner_catch_ugdrnoff);
% %     %remove day 1 from other var:
% %     realistic_cleaner_catch_albedo(1)=[];
% %     realistic_cleaner_catch_swe(1)=[];
% %     realistic_cleaner_catch_lai(1)=[];
% %     realistic_cleaner_catch_soilm1(1)=[];
% %     realistic_cleaner_catch_soilm2(1)=[];
% %     realistic_cleaner_catch_soilm3(1)=[];
% %     realistic_cleaner_catch_soilm4(1)=[];
% %     realistic_cleaner_catch_soilm_total(1)=[];
% %     
% %     %get monthly averages:
% %     realistic_cleaner_catch_monthly_albedo = accumarray(j,realistic_cleaner_catch_albedo,[],@nanmean);
% %     realistic_cleaner_catch_monthly_ecan = accumarray(j,realistic_cleaner_catch_ecan,[],@nanmean);
% %     realistic_cleaner_catch_monthly_etran = accumarray(j,realistic_cleaner_catch_etran,[],@nanmean);
% %     realistic_cleaner_catch_monthly_edir = accumarray(j,realistic_cleaner_catch_edir,[],@nanmean);
% %     realistic_cleaner_catch_monthly_et = accumarray(j,realistic_cleaner_catch_et,[],@nanmean);
% %     realistic_cleaner_catch_monthly_swe = accumarray(j,realistic_cleaner_catch_swe,[],@nanmean);
% %     realistic_cleaner_catch_monthly_lai = accumarray(j,realistic_cleaner_catch_lai,[],@nanmean);
% %     realistic_cleaner_catch_monthly_ugdrnoff = accumarray(j,realistic_cleaner_catch_ugdrnoff,[],@nanmean);
% %     realistic_cleaner_catch_monthly_soilm1 = accumarray(j,realistic_cleaner_catch_soilm1,[],@nanmean);
% %     realistic_cleaner_catch_monthly_soilm2 = accumarray(j,realistic_cleaner_catch_soilm2,[],@nanmean);
% %     realistic_cleaner_catch_monthly_soilm3 = accumarray(j,realistic_cleaner_catch_soilm3,[],@nanmean);
% %     realistic_cleaner_catch_monthly_soilm4 = accumarray(j,realistic_cleaner_catch_soilm4,[],@nanmean);
% %     realistic_cleaner_catch_monthly_soilm_total = accumarray(j,realistic_cleaner_catch_soilm_total,[],@nanmean);
% %     
% %     %% plot data
% %     %plot comparisson of ET and its paritions:
% %     f=figure;
% %     f.Position = [-1919        -171        1920         976];
% %     subplot(2,2,1)
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_et,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_et,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_et,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_et,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_et,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_et,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_et,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_et,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('ET (mm/day)','fontsize',24)
% %     ylim([0 4])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
% %     
% %     subplot(2,2,2)
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_etran,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_etran,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_etran,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_etran,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_etran,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_etran,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_etran,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_etran,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('E_{tran} (mm/day)','fontsize',24)
% %     ylim([0 4])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     
% %     subplot(2,2,3)
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_edir,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_edir,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_edir,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_edir,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_edir,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_edir,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_edir,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_edir,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('E_{soil} (mm/day)','fontsize',24)
% %     ylim([0 4])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     
% %     subplot(2,2,4)
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_ecan,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_ecan,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_ecan,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_ecan,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_ecan,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_ecan,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_ecan,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_ecan,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('E_{can} (mm/day)','fontsize',24)
% %     ylim([0 4])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     outfilename = sprintf('Catch%d_ET_comparisson_BurnedPixels.png',c);
% %     saveas(f,[outdir,outfilename])
% %     
% %     %plot comparisson of SM total and its layers:
% %     f=figure;
% %     f.Position = [-1919        -171        1920         976];
% %     subplot(4,2,[1:4])
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_soilm_total,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_soilm_total,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm_total,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm_total,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_soilm_total,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_soilm_total,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_soilm_total,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_soilm_total,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('total column SM (mm/mm)','fontsize',24)
% %     ylim([0.1 0.5])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
% %     
% %     subplot(4,2,5)
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_soilm1,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_soilm1,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm1,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm1,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_soilm1,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_soilm1,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_soilm1,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_soilm1,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel({'soilm1'; '(0-10cm)'},'fontsize',24)
% %     ylim([0.1 0.5])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     
% %     subplot(4,2,6)
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_soilm2,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_soilm2,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm2,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm2,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_soilm2,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_soilm2,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_soilm2,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_soilm2,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel({'soilm2'; '(10-40cm)'},'fontsize',24)
% %     ylim([0.1 0.5])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     
% %     subplot(4,2,7)
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_soilm3,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_soilm3,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm3,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm3,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_soilm3,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_soilm3,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_soilm3,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_soilm3,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel({'soilm3'; '(40-100cm)'},'fontsize',24)
% %     ylim([0.1 0.5])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     
% %     subplot(4,2,8)
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_soilm4,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_soilm4,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_soilm4,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_soilm4,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_soilm4,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_soilm4,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_soilm4,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_soilm4,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel({'soilm4'; '(100-200cm)'},'fontsize',24)
% %     ylim([0.05 0.5])
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     outfilename = sprintf('Catch%d_SOILM_comparisson_BurnedPixels.png',c);
% %     saveas(f,[outdir,outfilename])
% %     
% %     %plot comparisson of SWE
% %     f=figure;
% %     f.Position = [-1919         418        1920         387];
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_swe,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_swe,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_swe,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_swe,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_swe,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_swe,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_swe,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_swe,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('monthly SWE (mm)','fontsize',24)
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
% %     outfilename = sprintf('Catch%d_SWE_comparisson_BurnedPixels.png',c);
% %     saveas(f,[outdir,outfilename])
% %     
% %     %plot comparisson of lai
% %     f=figure;
% %     f.Position = [-1919         418        1920         387];
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_lai,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_lai,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_lai,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_lai,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_lai,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_lai,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_lai,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_lai,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('monthly LAI (mm/mm)','fontsize',24)
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
% %     outfilename = sprintf('Catch%d_lai_comparisson_BurnedPixels.png',c);
% %     saveas(f,[outdir,outfilename])
% %     
% %     %plot comparisson of ugdrnoff
% %     f=figure;
% %     f.Position = [-1919         418        1920         387];
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_ugdrnoff,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_ugdrnoff,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_ugdrnoff,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_ugdrnoff,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_ugdrnoff,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_ugdrnoff,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_ugdrnoff,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_ugdrnoff,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('monthly baseflow (mm/day)','fontsize',24)
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
% %     outfilename = sprintf('Catch%d_baseflow_comparisson_BurnedPixels.png',c);
% %     saveas(f,[outdir,outfilename])
% %     
% %     %plot comparisson of albedo
% %     f=figure;
% %     f.Position = [-1919         418        1920         387];
% %     hold on
% %     p1 = plot(monthly_dates,baseline_catch_monthly_albedo,'-','linewidth',2,'color','k');
% %     p2 = plot(monthly_dates, modified_param_catch_monthly_albedo,'-','linewidth',2,'color','b');
% %     p3 = plot(monthly_dates, modified_param_GVF_catch_monthly_albedo,'-','linewidth',2,'color',[0 0.5 0]);
% %     p4 = plot(monthly_dates, modified_param_GVF_VegClass_catch_monthly_albedo,'--','linewidth',2,'color',[0 0.5 0]);
% %     p5 = plot(monthly_dates, modified_param_albedo_GRASS_catch_monthly_albedo,'-','linewidth',2,'color','r');
% %     p6 = plot(monthly_dates, modified_param_albedo_BARE_catch_monthly_albedo,'--','linewidth',2,'color','r');
% %     p7 = plot(monthly_dates, realistic_catch_monthly_albedo,'-','linewidth',2,'color',[0.5 0.5 0.5]);
% %     p8 = plot(monthly_dates, realistic_cleaner_catch_monthly_albedo,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% %     grid on
% %     box on
% %     set(gca,'fontsize',24)
% %     datetick('x','mm/yyyy','keepticks','keeplimits')
% %     ylabel('Albedo','fontsize',24)
% %     xlim([min(monthly_dates) max(monthly_dates)])
% %     legend([p1 p2 p3 p4 p5],{'baseline','param-adjusted','param + GVF','param + GVF + veg-class','realistic'},'fontsize',16,'location','northoutside')
% %     outfilename = sprintf('Catch%d_albedo_comparisson_BurnedPixels.png',c);
% %     saveas(f,[outdir,outfilename])
% % end
