clc;clear all;close all;

%get burn severity data:
BSI = ncread('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile//geogrid_bsi.nc','bsi_merged');

%% load in shapefiles for the catchments:
%load in study shapes:
Catchments = shaperead('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
Catchments = Catchments(1:3);
shape_info_catchments = shapeinfo('/Users/abolafia/ASO_Fire/Data/Shapefiles_From_Aubrey/analysis/Catchments_of_interest/HUC8_CA_Simplified.shp');
p1_catchments = shape_info_catchments.CoordinateReferenceSystem;
ncatch = length(Catchments);

%% get domian grid:
LAT = ncread('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/wrfinput_d0x.nc','XLAT');
LONG = ncread('/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Inputs/geofile/wrfinput_d0x.nc','XLONG');
latvec = LAT(:);
lonvec =LONG(:);

%% define indices for AORC domain in each catchment:
store_IN_idx=[];
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
    [in] = inpolygon(latvec,lonvec,catch_lat,catch_lon);
    IDX=in;
    idx=find(IDX==1);
    idx = unique(idx);
    assert( max(idx) <= length(latvec), 'bad idx')

    store_IN_idx = [store_IN_idx;idx];
end
Total_IDX = 1:length(latvec);
IDX_out = Total_IDX;
IDX_out(store_IN_idx)=[];

%% load in Noah-MP
Baseline_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/Feather_Baseline/';
Realistic_dir = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/NoahMP/Outputs/For_paper/Realistic/';

store_Baseline_fall_median=[];
store_Baseline_spring_median=[];
store_Baseline_summer_median=[];
store_Baseline_winter_median=[];

store_Realistic_fall_median=[];
store_Realistic_spring_median=[];
store_Realistic_summer_median=[];
store_Realistic_winter_median=[];
for WY = 2000:2022
    infilename = sprintf('daily_WY%04d.nc',WY)
    WY_dates = [datenum([WY-1 10 1]):datenum([WY 9 30])];
    WY_dates= WY_dates(2:end);
    WY_datevecs = datevec(WY_dates);
    
    idx_fall = find(WY_datevecs(:,2) == 9 | WY_datevecs(:,2) ==10 | WY_datevecs(:,2) ==11);
    idx_winter = find(WY_datevecs(:,2) == 12 | WY_datevecs(:,2) ==1 | WY_datevecs(:,2) ==2);
    idx_spring = find(WY_datevecs(:,2) == 3 | WY_datevecs(:,2) ==4 | WY_datevecs(:,2) ==5);
    idx_summer = find(WY_datevecs(:,2) == 6 | WY_datevecs(:,2) ==7 | WY_datevecs(:,2) ==8);
    
    Baseline_data = ncread([Baseline_dir,infilename],'ET');
    Realistic_data = ncread([Realistic_dir,infilename],'ET');
    %change accumulated ET into just ET:
    Baseline_data = diff(Baseline_data,[],3);
    Realistic_data = diff(Realistic_data,[],3);
    
    %get seasonal averages:
    %Baseline
    Baseline_fall = Baseline_data(:,:,idx_fall);
    Baseline_spring = Baseline_data(:,:,idx_spring);
    Baseline_summer = Baseline_data(:,:,idx_summer);
    Baseline_winter = Baseline_data(:,:,idx_winter);
    
    Baseline_fall_median = squeeze(nanmedian(Baseline_fall,3));
    Baseline_spring_median = squeeze(nanmedian(Baseline_spring,3));
    Baseline_summer_median = squeeze(nanmedian(Baseline_summer,3));
    Baseline_winter_median = squeeze(nanmedian(Baseline_winter,3));
    
    store_Baseline_fall_median = cat(3,store_Baseline_fall_median,Baseline_fall_median);
    store_Baseline_spring_median = cat(3,store_Baseline_spring_median,Baseline_spring_median);
    store_Baseline_summer_median = cat(3,store_Baseline_summer_median,Baseline_summer_median);
    store_Baseline_winter_median = cat(3,store_Baseline_winter_median,Baseline_winter_median);
    
    %Realistic
    Realistic_fall = Realistic_data(:,:,idx_fall);
    Realistic_spring = Realistic_data(:,:,idx_spring);
    Realistic_summer = Realistic_data(:,:,idx_summer);
    Realistic_winter = Realistic_data(:,:,idx_winter);
    
    Realistic_fall_median = squeeze(nanmedian(Realistic_fall,3));
    Realistic_spring_median = squeeze(nanmedian(Realistic_spring,3));
    Realistic_summer_median = squeeze(nanmedian(Realistic_summer,3));
    Realistic_winter_median = squeeze(nanmedian(Realistic_winter,3));
    
    store_Realistic_fall_median = cat(3,store_Realistic_fall_median,Realistic_fall_median);
    store_Realistic_spring_median = cat(3,store_Realistic_spring_median,Realistic_spring_median);
    store_Realistic_summer_median = cat(3,store_Realistic_summer_median,Realistic_summer_median);
    store_Realistic_winter_median = cat(3,store_Realistic_winter_median,Realistic_winter_median);
end

%% take multiyear averages:
Baseline_fall_median = squeeze(nanmedian(store_Baseline_fall_median,3));
Baseline_winter_median = squeeze(nanmedian(store_Baseline_winter_median,3));
Baseline_spring_median = squeeze(nanmedian(store_Baseline_spring_median,3));
Baseline_summer_median = squeeze(nanmedian(store_Baseline_summer_median,3));

Realistic_fall_median = squeeze(nanmedian(store_Realistic_fall_median,3));
Realistic_winter_median = squeeze(nanmedian(store_Realistic_winter_median,3));
Realistic_spring_median = squeeze(nanmedian(store_Realistic_spring_median,3));
Realistic_summer_median = squeeze(nanmedian(store_Realistic_summer_median,3));

%compute differences:
delta_fall = Realistic_fall_median - Baseline_fall_median;
delta_winter = Realistic_winter_median - Baseline_winter_median;
delta_spring = Realistic_spring_median - Baseline_spring_median;
delta_summer = Realistic_summer_median - Baseline_summer_median;


%% create spatial plot for baseline fall data
PLOT_DATA = Baseline_fall_median;
PLOT_DATA(IDX_out) = NaN;

valid_min=0.01;
valid_max=3;

levels = linspace(valid_min, valid_max, 24);
levels = round(levels,2);
cmap1 = flipud(copper(12));
cmap2 = flipud(summer(13));

cmap = [cmap2([13:-1:1],:);cmap1([1:12],:)];
cmap = flipud(cmap);
Z = PLOT_DATA;

% Clamp the min and max values to the level index.
IDX_large = find(Z > levels(end));
if length(IDX_large) > 0
    Z(IDX_large) = length(levels)+1;
    iter=0;
else
    iter = 1;
     Z(iter) = length(levels)+1;
end
Z(Z < levels(1)) = 1;

% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 2:length(levels)
    IDX = find(PLOT_DATA >= levels(k-1) & PLOT_DATA < levels(k));
    if length(IDX) > 0
        Z(IDX) = double(k) ;
    else
        iter = iter+1;
        Z(iter) = double(k);
    end
end

f=figure;
f.Position = [-1854         -14         952         811];
lonlim = [-122 -120];
latlim = [39.2 40.6];

xlim(lonlim);
ylim(latlim);

hold on
% %tightmap
colormap(cmap)
geoshow(LAT,LONG,Z,'DisplayType','surface')

caxis auto
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
clevels{1} = ['<',num2str(levels(1))]; 
clevels{2} = [num2str((levels(1) + levels(2))/2 )]; 
clevels{3} = [num2str((levels(2) + levels(3))/2 )]; 
clevels{4} = [num2str((levels(3) + levels(4))/2 )];  
clevels{5} = [num2str((levels(4) + levels(5))/2 )]; 
clevels{6} = [num2str((levels(5) + levels(6))/2 )];  
clevels{7} = [num2str((levels(6) + levels(7))/2 )]; 
clevels{8} = [num2str((levels(7) + levels(8))/2 )]; 
clevels{9} = [num2str((levels(8) + levels(9))/2 )]; 
clevels{10} = [num2str((levels(9) + levels(10))/2 )];  
clevels{11} = [num2str((levels(10) + levels(11))/2 )]; 
clevels{12} = [num2str((levels(11) + levels(12))/2 )]; 
clevels{13} = [num2str((levels(12) + levels(13))/2 )]; 
clevels{14} = [num2str((levels(13) + levels(14))/2 )]; 
clevels{15} = [num2str((levels(14) + levels(15))/2 )];  
clevels{16} = [num2str((levels(15) + levels(16))/2 )]; 
clevels{17} = [num2str((levels(16) + levels(17))/2 )]; 
clevels{18} = [num2str((levels(17) + levels(18))/2 )]; 
clevels{19} = [num2str((levels(18) + levels(19))/2 )]; 
clevels{20} = [num2str((levels(19) + levels(20))/2 )];  
clevels{21} = [num2str((levels(20) + levels(21))/2 )]; 
clevels{22} = [num2str((levels(21) + levels(22))/2 )]; 
clevels{23} = [num2str((levels(22) + levels(23))/2 )];  
clevels{24} = [num2str((levels(23) + levels(24))/2 )];  
clevels{25} = ['>',num2str(levels(24))]; 

h = lcolorbar(clevels, 'Location', 'vertical','fontsize',30,'title','mm/day');
set(gca,'fontsize',30)
h.Position(3)=0.02;

tick_spots=1:4:25;
for i=1:length(h.YTickLabel)
    idx=find(tick_spots==i);
    if length(idx)==0
        h.YTickLabel{i}='';
    elseif length(idx)==1 && i~=1 && i~=length(h.YTickLabel)
        h.YTickLabel{i}=num2str(round(str2num(h.YTickLabel{i}),2));
    end
end

%overlay catchment shapes:
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
   p1=plot(catch_lon,catch_lat,'-k','linewidth',3);
   p1.ZData = 1000.*ones(length(p1.YData),1);
end

xlabel('longitude')
ylabel('latitude')
%%title('Baseline Fall ET')

outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_Fall_ET.png');
saveas(f,outfilename)

%% create spatial plot for baseline winter data
PLOT_DATA = Baseline_winter_median;
PLOT_DATA(IDX_out) = NaN;

valid_min=0.01;
valid_max=3;

levels = linspace(valid_min, valid_max, 24);
levels = round(levels,2);
cmap1 = flipud(copper(12));
cmap2 = flipud(summer(13));

cmap = [cmap2([13:-1:1],:);cmap1([1:12],:)];
cmap = flipud(cmap);
Z = PLOT_DATA;

% Clamp the min and max values to the level index.
IDX_large = find(Z > levels(end));
if length(IDX_large) > 0
    Z(IDX_large) = length(levels)+1;
    iter=0;
else
    iter = 1;
     Z(iter) = length(levels)+1;
end
Z(Z < levels(1)) = 1;

% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 2:length(levels)
    IDX = find(PLOT_DATA >= levels(k-1) & PLOT_DATA < levels(k));
    if length(IDX) > 0
        Z(IDX) = double(k) ;
    else
        iter = iter+1;
        Z(iter) = double(k);
    end
end

f=figure;
f.Position = [-1854         -14         952         811];
lonlim = [-122 -120];
latlim = [39.2 40.6];

xlim(lonlim);
ylim(latlim);

hold on
% %tightmap
colormap(cmap)
geoshow(LAT,LONG,Z,'DisplayType','surface')

caxis auto
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
clevels{1} = ['<',num2str(levels(1))]; 
clevels{2} = [num2str((levels(1) + levels(2))/2 )]; 
clevels{3} = [num2str((levels(2) + levels(3))/2 )]; 
clevels{4} = [num2str((levels(3) + levels(4))/2 )];  
clevels{5} = [num2str((levels(4) + levels(5))/2 )]; 
clevels{6} = [num2str((levels(5) + levels(6))/2 )];  
clevels{7} = [num2str((levels(6) + levels(7))/2 )]; 
clevels{8} = [num2str((levels(7) + levels(8))/2 )]; 
clevels{9} = [num2str((levels(8) + levels(9))/2 )]; 
clevels{10} = [num2str((levels(9) + levels(10))/2 )];  
clevels{11} = [num2str((levels(10) + levels(11))/2 )]; 
clevels{12} = [num2str((levels(11) + levels(12))/2 )]; 
clevels{13} = [num2str((levels(12) + levels(13))/2 )]; 
clevels{14} = [num2str((levels(13) + levels(14))/2 )]; 
clevels{15} = [num2str((levels(14) + levels(15))/2 )];  
clevels{16} = [num2str((levels(15) + levels(16))/2 )]; 
clevels{17} = [num2str((levels(16) + levels(17))/2 )]; 
clevels{18} = [num2str((levels(17) + levels(18))/2 )]; 
clevels{19} = [num2str((levels(18) + levels(19))/2 )]; 
clevels{20} = [num2str((levels(19) + levels(20))/2 )];  
clevels{21} = [num2str((levels(20) + levels(21))/2 )]; 
clevels{22} = [num2str((levels(21) + levels(22))/2 )]; 
clevels{23} = [num2str((levels(22) + levels(23))/2 )];  
clevels{24} = [num2str((levels(23) + levels(24))/2 )];  
clevels{25} = ['>',num2str(levels(24))]; 

h = lcolorbar(clevels, 'Location', 'vertical','fontsize',30,'title','mm/day');
set(gca,'fontsize',30)
h.Position(3)=0.02;

tick_spots=1:4:25;
for i=1:length(h.YTickLabel)
    idx=find(tick_spots==i);
    if length(idx)==0
        h.YTickLabel{i}='';
    elseif length(idx)==1 && i~=1 && i~=length(h.YTickLabel)
        h.YTickLabel{i}=num2str(round(str2num(h.YTickLabel{i}),2));
    end
end

%overlay catchment shapes:
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
   p1=plot(catch_lon,catch_lat,'-k','linewidth',3);
   p1.ZData = 1000.*ones(length(p1.YData),1);
end

xlabel('longitude')
ylabel('latitude')
%%title('Baseline Winter ET')
outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_winter_ET.png');
saveas(f,outfilename)

%% create spatial plot for baseline spring data
PLOT_DATA = Baseline_spring_median;
PLOT_DATA(IDX_out) = NaN;

valid_min=0.01;
valid_max=3;

levels = linspace(valid_min, valid_max, 24);
levels = round(levels,2);
cmap1 = flipud(copper(12));
cmap2 = flipud(summer(13));

cmap = [cmap2([13:-1:1],:);cmap1([1:12],:)];
cmap = flipud(cmap);
Z = PLOT_DATA;

% Clamp the min and max values to the level index.
IDX_large = find(Z > levels(end));
if length(IDX_large) > 0
    Z(IDX_large) = length(levels)+1;
    iter=0;
else
    iter = 1;
     Z(iter) = length(levels)+1;
end
Z(Z < levels(1)) = 1;

% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 2:length(levels)
    IDX = find(PLOT_DATA >= levels(k-1) & PLOT_DATA < levels(k));
    if length(IDX) > 0
        Z(IDX) = double(k) ;
    else
        iter = iter+1;
        Z(iter) = double(k);
    end
end

f=figure;
f.Position = [-1854         -14         952         811];
lonlim = [-122 -120];
latlim = [39.2 40.6];

xlim(lonlim);
ylim(latlim);

hold on
% %tightmap
colormap(cmap)
geoshow(LAT,LONG,Z,'DisplayType','surface')

caxis auto
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
clevels{1} = ['<',num2str(levels(1))]; 
clevels{2} = [num2str((levels(1) + levels(2))/2 )]; 
clevels{3} = [num2str((levels(2) + levels(3))/2 )]; 
clevels{4} = [num2str((levels(3) + levels(4))/2 )];  
clevels{5} = [num2str((levels(4) + levels(5))/2 )]; 
clevels{6} = [num2str((levels(5) + levels(6))/2 )];  
clevels{7} = [num2str((levels(6) + levels(7))/2 )]; 
clevels{8} = [num2str((levels(7) + levels(8))/2 )]; 
clevels{9} = [num2str((levels(8) + levels(9))/2 )]; 
clevels{10} = [num2str((levels(9) + levels(10))/2 )];  
clevels{11} = [num2str((levels(10) + levels(11))/2 )]; 
clevels{12} = [num2str((levels(11) + levels(12))/2 )]; 
clevels{13} = [num2str((levels(12) + levels(13))/2 )]; 
clevels{14} = [num2str((levels(13) + levels(14))/2 )]; 
clevels{15} = [num2str((levels(14) + levels(15))/2 )];  
clevels{16} = [num2str((levels(15) + levels(16))/2 )]; 
clevels{17} = [num2str((levels(16) + levels(17))/2 )]; 
clevels{18} = [num2str((levels(17) + levels(18))/2 )]; 
clevels{19} = [num2str((levels(18) + levels(19))/2 )]; 
clevels{20} = [num2str((levels(19) + levels(20))/2 )];  
clevels{21} = [num2str((levels(20) + levels(21))/2 )]; 
clevels{22} = [num2str((levels(21) + levels(22))/2 )]; 
clevels{23} = [num2str((levels(22) + levels(23))/2 )];  
clevels{24} = [num2str((levels(23) + levels(24))/2 )];  
clevels{25} = ['>',num2str(levels(24))]; 

h = lcolorbar(clevels, 'Location', 'vertical','fontsize',30,'title','mm/day');
set(gca,'fontsize',30)
h.Position(3)=0.02;

tick_spots=1:4:25;
for i=1:length(h.YTickLabel)
    idx=find(tick_spots==i);
    if length(idx)==0
        h.YTickLabel{i}='';
    elseif length(idx)==1 && i~=1 && i~=length(h.YTickLabel)
        h.YTickLabel{i}=num2str(round(str2num(h.YTickLabel{i}),2));
    end
end

%overlay catchment shapes:
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
   p1=plot(catch_lon,catch_lat,'-k','linewidth',3);
   p1.ZData = 1000.*ones(length(p1.YData),1);
end

xlabel('longitude')
ylabel('latitude')
%%title('Baseline Spring ET')
outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_spring_ET.png');
saveas(f,outfilename)

%% create spatial plot for baseline summer data
PLOT_DATA = Baseline_summer_median;
PLOT_DATA(IDX_out) = NaN;

valid_min=0.01;
valid_max=3;

levels = linspace(valid_min, valid_max, 24);
levels = round(levels,2);
cmap1 = flipud(copper(12));
cmap2 = flipud(summer(13));

cmap = [cmap2([13:-1:1],:);cmap1([1:12],:)];
cmap = flipud(cmap);
Z = PLOT_DATA;

% Clamp the min and max values to the level index.
IDX_large = find(Z > levels(end));
if length(IDX_large) > 0
    Z(IDX_large) = length(levels)+1;
    iter=0;
else
    iter = 1;
     Z(iter) = length(levels)+1;
end
Z(Z < levels(1)) = 1;

% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 2:length(levels)
    IDX = find(PLOT_DATA >= levels(k-1) & PLOT_DATA < levels(k));
    if length(IDX) > 0
        Z(IDX) = double(k) ;
    else
        iter = iter+1;
        Z(iter) = double(k);
    end
end

f=figure;
f.Position = [-1854         -14         952         811];
lonlim = [-122 -120];
latlim = [39.2 40.6];

xlim(lonlim);
ylim(latlim);

hold on
% %tightmap
colormap(cmap)
geoshow(LAT,LONG,Z,'DisplayType','surface')

caxis auto
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
clevels{1} = ['<',num2str(levels(1))]; 
clevels{2} = [num2str((levels(1) + levels(2))/2 )]; 
clevels{3} = [num2str((levels(2) + levels(3))/2 )]; 
clevels{4} = [num2str((levels(3) + levels(4))/2 )];  
clevels{5} = [num2str((levels(4) + levels(5))/2 )]; 
clevels{6} = [num2str((levels(5) + levels(6))/2 )];  
clevels{7} = [num2str((levels(6) + levels(7))/2 )]; 
clevels{8} = [num2str((levels(7) + levels(8))/2 )]; 
clevels{9} = [num2str((levels(8) + levels(9))/2 )]; 
clevels{10} = [num2str((levels(9) + levels(10))/2 )];  
clevels{11} = [num2str((levels(10) + levels(11))/2 )]; 
clevels{12} = [num2str((levels(11) + levels(12))/2 )]; 
clevels{13} = [num2str((levels(12) + levels(13))/2 )]; 
clevels{14} = [num2str((levels(13) + levels(14))/2 )]; 
clevels{15} = [num2str((levels(14) + levels(15))/2 )];  
clevels{16} = [num2str((levels(15) + levels(16))/2 )]; 
clevels{17} = [num2str((levels(16) + levels(17))/2 )]; 
clevels{18} = [num2str((levels(17) + levels(18))/2 )]; 
clevels{19} = [num2str((levels(18) + levels(19))/2 )]; 
clevels{20} = [num2str((levels(19) + levels(20))/2 )];  
clevels{21} = [num2str((levels(20) + levels(21))/2 )]; 
clevels{22} = [num2str((levels(21) + levels(22))/2 )]; 
clevels{23} = [num2str((levels(22) + levels(23))/2 )];  
clevels{24} = [num2str((levels(23) + levels(24))/2 )];  
clevels{25} = ['>',num2str(levels(24))]; 

h = lcolorbar(clevels, 'Location', 'vertical','fontsize',30,'title','mm/day');
set(gca,'fontsize',30)
h.Position(3)=0.02;

tick_spots=1:4:25;
for i=1:length(h.YTickLabel)
    idx=find(tick_spots==i);
    if length(idx)==0
        h.YTickLabel{i}='';
    elseif length(idx)==1 && i~=1 && i~=length(h.YTickLabel)
        h.YTickLabel{i}=num2str(round(str2num(h.YTickLabel{i}),2));
    end
end

%overlay catchment shapes:
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
   p1=plot(catch_lon,catch_lat,'-k','linewidth',3);
   p1.ZData = 1000.*ones(length(p1.YData),1);
end

xlabel('longitude')
ylabel('latitude')
%%title('Baseline Summer ET')
outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Baseline_summer_ET.png');
saveas(f,outfilename)

%% create spatial plot for delta fall data
idx_burn_increase=find(BSI>=2 & delta_fall>0);
idx_burn=find(BSI>=2);
sprintf('%.4f%% of burn areas have Fall ET increases',length(idx_burn_increase)/length(idx_burn)*100)
PLOT_DATA = delta_fall;
PLOT_DATA(IDX_out) = NaN;

valid_min=-0.5;
valid_max=0.5;

levels = linspace(valid_min, valid_max, 24);
levels = round(levels,2);
cmap = flipud(polarmap(25));

Z = PLOT_DATA;

% Clamp the min and max values to the level index.
IDX_large = find(Z > levels(end));
if length(IDX_large) > 0
    Z(IDX_large) = length(levels)+1;
    iter=0;
else
    iter = 1;
     Z(iter) = length(levels)+1;
end

IDX_small = find(Z <= levels(1));
if length(IDX_small) > 0
    Z(IDX_small) = 1;
else
    if iter==0
        iter = 1;
        Z(iter) = 1;
    else
        iter = 2;
        Z(iter) = 1;
    end
end

% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 2:length(levels)
    IDX = find(PLOT_DATA >= levels(k-1) & PLOT_DATA < levels(k));
    if length(IDX) > 0
        Z(IDX) = double(k) ;
    else
        iter = iter+1;
        Z(iter) = double(k);
    end
end

f=figure;
f.Position = [-1854         -14         952         811];
lonlim = [-122 -120];
latlim = [39.2 40.6];

xlim(lonlim);
ylim(latlim);

hold on
% %tightmap
colormap(cmap)
geoshow(LAT,LONG,Z,'DisplayType','surface')

caxis auto
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
clevels{1} = ['<',num2str(levels(1))]; 
clevels{2} = [num2str((levels(1) + levels(2))/2 )]; 
clevels{3} = [num2str((levels(2) + levels(3))/2 )]; 
clevels{4} = [num2str((levels(3) + levels(4))/2 )];  
clevels{5} = [num2str((levels(4) + levels(5))/2 )]; 
clevels{6} = [num2str((levels(5) + levels(6))/2 )];  
clevels{7} = [num2str((levels(6) + levels(7))/2 )]; 
clevels{8} = [num2str((levels(7) + levels(8))/2 )]; 
clevels{9} = [num2str((levels(8) + levels(9))/2 )]; 
clevels{10} = [num2str((levels(9) + levels(10))/2 )];  
clevels{11} = [num2str((levels(10) + levels(11))/2 )]; 
clevels{12} = [num2str((levels(11) + levels(12))/2 )]; 
clevels{13} = [num2str((levels(12) + levels(13))/2 )]; 
clevels{14} = [num2str((levels(13) + levels(14))/2 )]; 
clevels{15} = [num2str((levels(14) + levels(15))/2 )];  
clevels{16} = [num2str((levels(15) + levels(16))/2 )]; 
clevels{17} = [num2str((levels(16) + levels(17))/2 )]; 
clevels{18} = [num2str((levels(17) + levels(18))/2 )]; 
clevels{19} = [num2str((levels(18) + levels(19))/2 )]; 
clevels{20} = [num2str((levels(19) + levels(20))/2 )];  
clevels{21} = [num2str((levels(20) + levels(21))/2 )]; 
clevels{22} = [num2str((levels(21) + levels(22))/2 )]; 
clevels{23} = [num2str((levels(22) + levels(23))/2 )];  
clevels{24} = [num2str((levels(23) + levels(24))/2 )];  
clevels{25} = ['>',num2str(levels(24))]; 

h = lcolorbar(clevels, 'Location', 'vertical','fontsize',30,'title','mm/day');
set(gca,'fontsize',30)
h.Position(3)=0.02;

tick_spots=1:4:25;
for i=1:length(h.YTickLabel)
    idx=find(tick_spots==i);
    if length(idx)==0
        h.YTickLabel{i}='';
    elseif length(idx)==1 && i~=1 && i~=length(h.YTickLabel)
        h.YTickLabel{i}=num2str(round(str2num(h.YTickLabel{i}),2));
    end
end
%overlay catchment shapes:
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
   p1=plot(catch_lon,catch_lat,'-k','linewidth',3);
   p1.ZData = 1000.*ones(length(p1.YData),1);
end

xlabel('longitude')
ylabel('latitude')
%%title('Realistic - Baseline Fall ET')
outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Delta_Fall_ET.png');
saveas(f,outfilename)

%% create spatial plot for delta winter data
idx_burn_increase=find(BSI>=2 & delta_winter>0);
idx_burn=find(BSI>=2);
sprintf('%.4f%% of burn areas have Winter ET increases',length(idx_burn_increase)/length(idx_burn)*100)
PLOT_DATA = delta_winter;
PLOT_DATA(IDX_out) = NaN;

valid_min=-0.5;
valid_max=0.5;

levels = linspace(valid_min, valid_max, 24);
levels = round(levels,2);
cmap = flipud(polarmap(25));

Z = PLOT_DATA;

% Clamp the min and max values to the level index.
IDX_large = find(Z > levels(end));
if length(IDX_large) > 0
    Z(IDX_large) = length(levels)+1;
    iter=0;
else
    iter = 1;
     Z(iter) = length(levels)+1;
end

IDX_small = find(Z <= levels(1));
if length(IDX_small) > 0
    Z(IDX_small) = 1;
else
    if iter==0
        iter = 1;
        Z(iter) = 1;
    else
        iter = 2;
        Z(iter) = 1;
    end
end
% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 2:length(levels)
    IDX = find(PLOT_DATA >= levels(k-1) & PLOT_DATA < levels(k));
    if length(IDX) > 0
        Z(IDX) = double(k) ;
    else
        iter = iter+1;
        Z(iter) = double(k);
    end
end

f=figure;
f.Position = [-1854         -14         952         811];
lonlim = [-122 -120];
latlim = [39.2 40.6];

xlim(lonlim);
ylim(latlim);

hold on
% %tightmap
colormap(cmap)
geoshow(LAT,LONG,Z,'DisplayType','surface')

caxis auto
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
clevels{1} = ['<',num2str(levels(1))]; 
clevels{2} = [num2str((levels(1) + levels(2))/2 )]; 
clevels{3} = [num2str((levels(2) + levels(3))/2 )]; 
clevels{4} = [num2str((levels(3) + levels(4))/2 )];  
clevels{5} = [num2str((levels(4) + levels(5))/2 )]; 
clevels{6} = [num2str((levels(5) + levels(6))/2 )];  
clevels{7} = [num2str((levels(6) + levels(7))/2 )]; 
clevels{8} = [num2str((levels(7) + levels(8))/2 )]; 
clevels{9} = [num2str((levels(8) + levels(9))/2 )]; 
clevels{10} = [num2str((levels(9) + levels(10))/2 )];  
clevels{11} = [num2str((levels(10) + levels(11))/2 )]; 
clevels{12} = [num2str((levels(11) + levels(12))/2 )]; 
clevels{13} = [num2str((levels(12) + levels(13))/2 )]; 
clevels{14} = [num2str((levels(13) + levels(14))/2 )]; 
clevels{15} = [num2str((levels(14) + levels(15))/2 )];  
clevels{16} = [num2str((levels(15) + levels(16))/2 )]; 
clevels{17} = [num2str((levels(16) + levels(17))/2 )]; 
clevels{18} = [num2str((levels(17) + levels(18))/2 )]; 
clevels{19} = [num2str((levels(18) + levels(19))/2 )]; 
clevels{20} = [num2str((levels(19) + levels(20))/2 )];  
clevels{21} = [num2str((levels(20) + levels(21))/2 )]; 
clevels{22} = [num2str((levels(21) + levels(22))/2 )]; 
clevels{23} = [num2str((levels(22) + levels(23))/2 )];  
clevels{24} = [num2str((levels(23) + levels(24))/2 )];  
clevels{25} = ['>',num2str(levels(24))]; 

h = lcolorbar(clevels, 'Location', 'vertical','fontsize',30,'title','mm/day');
set(gca,'fontsize',30)
h.Position(3)=0.02;

tick_spots=1:4:25;
for i=1:length(h.YTickLabel)
    idx=find(tick_spots==i);
    if length(idx)==0
        h.YTickLabel{i}='';
    elseif length(idx)==1 && i~=1 && i~=length(h.YTickLabel)
        h.YTickLabel{i}=num2str(round(str2num(h.YTickLabel{i}),2));
    end
end


%overlay catchment shapes:
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
   p1=plot(catch_lon,catch_lat,'-k','linewidth',3);
   p1.ZData = 1000.*ones(length(p1.YData),1);
end

xlabel('longitude')
ylabel('latitude')
%%title('Realistic - Baseline Winter ET')
outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Delta_winter_ET.png');
saveas(f,outfilename)

%% create spatial plot for delta spring data
idx_burn_increase=find(BSI>=2 & delta_spring>0);
idx_burn=find(BSI>=2);
sprintf('%.4f%% of burn areas have Spring ET increases',length(idx_burn_increase)/length(idx_burn)*100)
PLOT_DATA = delta_spring;
PLOT_DATA(IDX_out) = NaN;

valid_min=-0.5;
valid_max=0.5;

levels = linspace(valid_min, valid_max, 24);
levels = round(levels,2);
cmap = flipud(polarmap(25));

Z = PLOT_DATA;

% Clamp the min and max values to the level index.
IDX_large = find(Z > levels(end));
if length(IDX_large) > 0
    Z(IDX_large) = length(levels)+1;
    iter=0;
else
    iter = 1;
     Z(iter) = length(levels)+1;
end

IDX_small = find(Z <= levels(1));
if length(IDX_small) > 0
    Z(IDX_small) = 1;
else
    if iter==0
        iter = 1;
        Z(iter) = 1;
    else
        iter = 2;
        Z(iter) = 1;
    end
end

% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 2:length(levels)
    IDX = find(PLOT_DATA >= levels(k-1) & PLOT_DATA < levels(k));
    if length(IDX) > 0
        Z(IDX) = double(k) ;
    else
        iter = iter+1;
        Z(iter) = double(k);
    end
end

f=figure;
f.Position = [-1854         -14         952         811];
lonlim = [-122 -120];
latlim = [39.2 40.6];

xlim(lonlim);
ylim(latlim);

hold on
% %tightmap
colormap(cmap)
geoshow(LAT,LONG,Z,'DisplayType','surface')

caxis auto
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
clevels{1} = ['<',num2str(levels(1))]; 
clevels{2} = [num2str((levels(1) + levels(2))/2 )]; 
clevels{3} = [num2str((levels(2) + levels(3))/2 )]; 
clevels{4} = [num2str((levels(3) + levels(4))/2 )];  
clevels{5} = [num2str((levels(4) + levels(5))/2 )]; 
clevels{6} = [num2str((levels(5) + levels(6))/2 )];  
clevels{7} = [num2str((levels(6) + levels(7))/2 )]; 
clevels{8} = [num2str((levels(7) + levels(8))/2 )]; 
clevels{9} = [num2str((levels(8) + levels(9))/2 )]; 
clevels{10} = [num2str((levels(9) + levels(10))/2 )];  
clevels{11} = [num2str((levels(10) + levels(11))/2 )]; 
clevels{12} = [num2str((levels(11) + levels(12))/2 )]; 
clevels{13} = [num2str((levels(12) + levels(13))/2 )]; 
clevels{14} = [num2str((levels(13) + levels(14))/2 )]; 
clevels{15} = [num2str((levels(14) + levels(15))/2 )];  
clevels{16} = [num2str((levels(15) + levels(16))/2 )]; 
clevels{17} = [num2str((levels(16) + levels(17))/2 )]; 
clevels{18} = [num2str((levels(17) + levels(18))/2 )]; 
clevels{19} = [num2str((levels(18) + levels(19))/2 )]; 
clevels{20} = [num2str((levels(19) + levels(20))/2 )];  
clevels{21} = [num2str((levels(20) + levels(21))/2 )]; 
clevels{22} = [num2str((levels(21) + levels(22))/2 )]; 
clevels{23} = [num2str((levels(22) + levels(23))/2 )];  
clevels{24} = [num2str((levels(23) + levels(24))/2 )];  
clevels{25} = ['>',num2str(levels(24))]; 

h = lcolorbar(clevels, 'Location', 'vertical','fontsize',30,'title','mm/day');
set(gca,'fontsize',30)
h.Position(3)=0.02;

tick_spots=1:4:25;
for i=1:length(h.YTickLabel)
    idx=find(tick_spots==i);
    if length(idx)==0
        h.YTickLabel{i}='';
    elseif length(idx)==1 && i~=1 && i~=length(h.YTickLabel)
        h.YTickLabel{i}=num2str(round(str2num(h.YTickLabel{i}),2));
    end
end
%overlay catchment shapes:
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
   p1=plot(catch_lon,catch_lat,'-k','linewidth',3);
   p1.ZData = 1000.*ones(length(p1.YData),1);
end

xlabel('longitude')
ylabel('latitude')
%%title('Realistic - Baseline Spring ET')
outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Delta_spring_ET.png');
saveas(f,outfilename)
%save this output for deeper analysis:
outfilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/Catchment_averaged_outputs/Gridded_deltaET_Spring.mat';
save(outfilename,'delta_spring');

outfilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/Catchment_averaged_outputs/Gridded_Baseline_ET_Spring.mat';
save(outfilename,'store_Baseline_spring_median')

outfilename = '/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/Catchment_averaged_outputs/Gridded_Realistic_ET_Spring.mat';
save(outfilename,'store_Realistic_spring_median')
%% create spatial plot for delta summer data
idx_burn_increase=find(BSI>=2 & delta_summer>0);
idx_burn=find(BSI>=2);
sprintf('%.4f%% of burn areas have summer ET increases',length(idx_burn_increase)/length(idx_burn)*100)
PLOT_DATA = delta_summer;
PLOT_DATA(IDX_out) = NaN;

valid_min=-0.5;
valid_max=0.5;

levels = linspace(valid_min, valid_max, 24);
levels = round(levels,2);
cmap = flipud(polarmap(25));

Z = PLOT_DATA;

% Clamp the min and max values to the level index.
IDX_large = find(Z > levels(end));
if length(IDX_large) > 0
    Z(IDX_large) = length(levels)+1;
    iter=0;
else
    iter = 1;
     Z(iter) = length(levels)+1;
end

IDX_small = find(Z <= levels(1));
if length(IDX_small) > 0
    Z(IDX_small) = 1;
else
    if iter==0
        iter = 1;
        Z(iter) = 1;
    else
        iter = 2;
        Z(iter) = 1;
    end
end

% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 2:length(levels)
    IDX = find(PLOT_DATA >= levels(k-1) & PLOT_DATA < levels(k));
    if length(IDX) > 0
        Z(IDX) = double(k) ;
    else
        iter = iter+1;
        Z(iter) = double(k);
    end
end

f=figure;
f.Position = [-1854         -14         952         811];
lonlim = [-122 -120];
latlim = [39.2 40.6];

xlim(lonlim);
ylim(latlim);

hold on
% %tightmap
colormap(cmap)
geoshow(LAT,LONG,Z,'DisplayType','surface')

caxis auto
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
clevels{1} = ['<',num2str(levels(1))]; 
clevels{2} = [num2str((levels(1) + levels(2))/2 )]; 
clevels{3} = [num2str((levels(2) + levels(3))/2 )]; 
clevels{4} = [num2str((levels(3) + levels(4))/2 )];  
clevels{5} = [num2str((levels(4) + levels(5))/2 )]; 
clevels{6} = [num2str((levels(5) + levels(6))/2 )];  
clevels{7} = [num2str((levels(6) + levels(7))/2 )]; 
clevels{8} = [num2str((levels(7) + levels(8))/2 )]; 
clevels{9} = [num2str((levels(8) + levels(9))/2 )]; 
clevels{10} = [num2str((levels(9) + levels(10))/2 )];  
clevels{11} = [num2str((levels(10) + levels(11))/2 )]; 
clevels{12} = [num2str((levels(11) + levels(12))/2 )]; 
clevels{13} = [num2str((levels(12) + levels(13))/2 )]; 
clevels{14} = [num2str((levels(13) + levels(14))/2 )]; 
clevels{15} = [num2str((levels(14) + levels(15))/2 )];  
clevels{16} = [num2str((levels(15) + levels(16))/2 )]; 
clevels{17} = [num2str((levels(16) + levels(17))/2 )]; 
clevels{18} = [num2str((levels(17) + levels(18))/2 )]; 
clevels{19} = [num2str((levels(18) + levels(19))/2 )]; 
clevels{20} = [num2str((levels(19) + levels(20))/2 )];  
clevels{21} = [num2str((levels(20) + levels(21))/2 )]; 
clevels{22} = [num2str((levels(21) + levels(22))/2 )]; 
clevels{23} = [num2str((levels(22) + levels(23))/2 )];  
clevels{24} = [num2str((levels(23) + levels(24))/2 )];  
clevels{25} = ['>',num2str(levels(24))]; 

h = lcolorbar(clevels, 'Location', 'vertical','fontsize',30,'title','mm/day');
set(gca,'fontsize',30)
h.Position(3)=0.02;

tick_spots=1:4:25;
for i=1:length(h.YTickLabel)
    idx=find(tick_spots==i);
    if length(idx)==0
        h.YTickLabel{i}='';
    elseif length(idx)==1 && i~=1 && i~=length(h.YTickLabel)
        h.YTickLabel{i}=num2str(round(str2num(h.YTickLabel{i}),2));
    end
end
%overlay catchment shapes:
for c=1:ncatch
    current_shape =  Catchments(c);
    Shape_X=current_shape.X;
    Shape_Y = current_shape.Y;
    [catch_lat,catch_lon] = projinv(p1_catchments,Shape_X,Shape_Y);
    idx=find(catch_lat==90);
    catch_lat(idx)=[];
    catch_lon(idx)=[];
    
   p1=plot(catch_lon,catch_lat,'-k','linewidth',3);
   p1.ZData = 1000.*ones(length(p1.YData),1);
end

xlabel('longitude')
ylabel('latitude')
%%title('Realistic - Baseline Summer ET')
outfilename = sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/Delta_summer_ET.png');
saveas(f,outfilename)