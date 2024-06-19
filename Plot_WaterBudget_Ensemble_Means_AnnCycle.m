clc;clear all;close all;

%% load in outputs from:Compare_Baseline_and_ModParam_WaterBudget.m
catch_names = {'East_Branch_NF_Feather','Mid_Fork_Feather','NF_Feather'};

for c=1:3
    baseline_filename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs//baseline_LSM_outputs_Catch_%d.mat',c);
    mod_param_filename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs//modified_param_LSM_outputs_Catch_%d.mat',c);
    modified_param_GVF_filename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs//modified_param_GVF_LSM_outputs_Catch_%d.mat',c);
    modified_param_GVF_VegClass_filename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs//modified_param_GVF_VegClass_LSM_outputs_Catch_%d.mat',c);
    realistic_filename = sprintf('/Volumes/Pruina_External_Elements/ASO_Fire/Data/Analysis_Data/For_Paper//Catchment_averaged_outputs//realistic_LSM_outputs_Catch_%d.mat',c);

    %baseline data
    baseline_data = load(baseline_filename);
    baseline_data = baseline_data.Catchment_outputs;
    dates = baseline_data.date;
    dates = double(dates);
    baseline_ET = baseline_data.et;
    baseline_edir = baseline_data.edir;
    baseline_etran = baseline_data.etran;
    baseline_ecan = baseline_data.ecan;
    baseline_SWE = baseline_data.swe;
    baseline_ET = diff(baseline_ET);
    baseline_edir = diff(baseline_edir);
    baseline_etran = diff(baseline_etran);
    baseline_ecan = diff(baseline_ecan);

    baseline_soilm1 = baseline_data.soilm1;
    baseline_soilm2 = baseline_data.soilm2;
    baseline_soilm3 = baseline_data.soilm3;
    baseline_soilm4 = baseline_data.soilm4;
    baseline_soilm = baseline_soilm1.*(0.1/2) + baseline_soilm2.*(0.3/2) + baseline_soilm4.*(0.6/2) + baseline_soilm4.*(1/2);
    baseline_baseflow = baseline_data.ugdrnoff;
    baseline_baseflow = diff(baseline_baseflow);

    %get annual cycles:
    datevecs = datevec(dates);
    datevecs(1,:)=[];
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    baseline_ET(idx)=[];
    baseline_edir(idx)=[];
    baseline_etran(idx)=[];
    baseline_ecan(idx)=[];
    baseline_baseflow(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_ET_baseline = accumarray(j,baseline_ET,[],@nanmean);
    annual_ET_baseline = [annual_ET_baseline(end);annual_ET_baseline(1:end-1)];

    annual_edir_baseline = accumarray(j,baseline_edir,[],@nanmean);
    annual_edir_baseline = [annual_edir_baseline(end);annual_edir_baseline(1:end-1)];

    annual_etran_baseline = accumarray(j,baseline_etran,[],@nanmean);
    annual_etran_baseline = [annual_etran_baseline(end);annual_etran_baseline(1:end-1)];

    annual_ecan_baseline = accumarray(j,baseline_ecan,[],@nanmean);
    annual_ecan_baseline = [annual_ecan_baseline(end);annual_ecan_baseline(1:end-1)];

    annual_baseflow_baseline = accumarray(j,baseline_baseflow,[],@nanmean);
    annual_baseflow_baseline = [annual_baseflow_baseline(end);annual_baseflow_baseline(1:end-1)];


    datevecs = datevec(dates);
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    baseline_SWE(idx)=[];
    baseline_soilm(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_SWE_baseline = accumarray(j,baseline_SWE,[],@nanmean);
    annual_SM_baseline = accumarray(j,baseline_soilm,[],@nanmean);

    %mod_param data
    mod_param_data = load(mod_param_filename);
    mod_param_data = mod_param_data.Catchment_outputs;
    mod_param_ET = mod_param_data.et;
    mod_param_SWE = mod_param_data.swe;
    mod_param_ET = diff(mod_param_ET);
    mod_param_edir = mod_param_data.edir;
    mod_param_etran = mod_param_data.etran;
    mod_param_ecan = mod_param_data.ecan;
    mod_param_edir = diff(mod_param_edir);
    mod_param_etran = diff(mod_param_etran);
    mod_param_ecan = diff(mod_param_ecan);

    mod_param_soilm1 = mod_param_data.soilm1;
    mod_param_soilm2 = mod_param_data.soilm2;
    mod_param_soilm3 = mod_param_data.soilm3;
    mod_param_soilm4 = mod_param_data.soilm4;
    mod_param_soilm = mod_param_soilm1.*(0.1/2) + mod_param_soilm2.*(0.3/2) + mod_param_soilm4.*(0.6/2) + mod_param_soilm4.*(1/2);
    mod_param_baseflow = mod_param_data.ugdrnoff;
    mod_param_baseflow = diff(mod_param_baseflow);
    %get annual cycles:
    datevecs = datevec(dates);
    datevecs(1,:)=[];
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    mod_param_ET(idx)=[];
    mod_param_edir(idx)=[];
    mod_param_etran(idx)=[];
    mod_param_ecan(idx)=[];
    mod_param_baseflow(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_ET_mod_param = accumarray(j,mod_param_ET,[],@nanmean);
    annual_ET_mod_param = [annual_ET_mod_param(end);annual_ET_mod_param(1:end-1)];

    annual_edir_mod_param = accumarray(j,mod_param_edir,[],@nanmean);
    annual_edir_mod_param = [annual_edir_mod_param(end);annual_edir_mod_param(1:end-1)];

    annual_etran_mod_param = accumarray(j,mod_param_etran,[],@nanmean);
    annual_etran_mod_param = [annual_etran_mod_param(end);annual_etran_mod_param(1:end-1)];

    annual_ecan_mod_param = accumarray(j,mod_param_ecan,[],@nanmean);
    annual_ecan_mod_param = [annual_ecan_mod_param(end);annual_ecan_mod_param(1:end-1)];

    annual_baseflow_mod_param = accumarray(j,mod_param_baseflow,[],@nanmean);
    annual_baseflow_mod_param = [annual_baseflow_mod_param(end);annual_baseflow_mod_param(1:end-1)];

    datevecs = datevec(dates);
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    mod_param_SWE(idx)=[];
    mod_param_soilm(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_SWE_mod_param = accumarray(j,mod_param_SWE,[],@nanmean);
    annual_SM_mod_param = accumarray(j,mod_param_soilm,[],@nanmean);

    %modified_param_GVF data
    modified_param_GVF_data = load(modified_param_GVF_filename);
    modified_param_GVF_data = modified_param_GVF_data.Catchment_outputs;
    modified_param_GVF_ET = modified_param_GVF_data.et;
    modified_param_GVF_SWE = modified_param_GVF_data.swe;
    modified_param_GVF_ET = diff(modified_param_GVF_ET);
    modified_param_GVF_edir = modified_param_GVF_data.edir;
    modified_param_GVF_etran = modified_param_GVF_data.etran;
    modified_param_GVF_ecan = modified_param_GVF_data.ecan;

    modified_param_GVF_edir = diff(modified_param_GVF_edir);
    modified_param_GVF_etran = diff(modified_param_GVF_etran);
    modified_param_GVF_ecan = diff(modified_param_GVF_ecan);

    modified_param_GVF_soilm1 = modified_param_GVF_data.soilm1;
    modified_param_GVF_soilm2 = modified_param_GVF_data.soilm2;
    modified_param_GVF_soilm3 = modified_param_GVF_data.soilm3;
    modified_param_GVF_soilm4 = modified_param_GVF_data.soilm4;
    modified_param_GVF_soilm = modified_param_GVF_soilm1.*(0.1/2) + modified_param_GVF_soilm2.*(0.3/2) + modified_param_GVF_soilm4.*(0.6/2) + modified_param_GVF_soilm4.*(1/2);
    modified_param_GVF_baseflow = modified_param_GVF_data.ugdrnoff;
    modified_param_GVF_baseflow = diff(modified_param_GVF_baseflow);
    %get annual cycles:
    datevecs = datevec(dates);
    datevecs(1,:)=[];
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    modified_param_GVF_ET(idx)=[];
    modified_param_GVF_edir(idx)=[];
    modified_param_GVF_etran(idx)=[];
    modified_param_GVF_ecan(idx)=[];
    modified_param_GVF_baseflow(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_ET_modified_param_GVF = accumarray(j,modified_param_GVF_ET,[],@nanmean);
    annual_ET_modified_param_GVF = [annual_ET_modified_param_GVF(end);annual_ET_modified_param_GVF(1:end-1)];

    annual_edir_modified_param_GVF = accumarray(j,modified_param_GVF_edir,[],@nanmean);
    annual_edir_modified_param_GVF = [annual_edir_modified_param_GVF(end);annual_edir_modified_param_GVF(1:end-1)];

    annual_etran_modified_param_GVF = accumarray(j,modified_param_GVF_etran,[],@nanmean);
    annual_etran_modified_param_GVF = [annual_etran_modified_param_GVF(end);annual_etran_modified_param_GVF(1:end-1)];

    annual_ecan_modified_param_GVF = accumarray(j,modified_param_GVF_ecan,[],@nanmean);
    annual_ecan_modified_param_GVF = [annual_ecan_modified_param_GVF(end);annual_ecan_modified_param_GVF(1:end-1)];

    annual_baseflow_modified_param_GVF = accumarray(j,modified_param_GVF_baseflow,[],@nanmean);
    annual_baseflow_modified_param_GVF = [annual_baseflow_modified_param_GVF(end);annual_baseflow_modified_param_GVF(1:end-1)];

    datevecs = datevec(dates);
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    modified_param_GVF_SWE(idx)=[];
    modified_param_GVF_soilm(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_SWE_modified_param_GVF = accumarray(j,modified_param_GVF_SWE,[],@nanmean);
    annual_SM_modified_param_GVF = accumarray(j,modified_param_GVF_soilm,[],@nanmean);

    %modified_param_GVF_VegClass data
    modified_param_GVF_VegClass_data = load(modified_param_GVF_VegClass_filename);
    modified_param_GVF_VegClass_data = modified_param_GVF_VegClass_data.Catchment_outputs;
    modified_param_GVF_VegClass_ET = modified_param_GVF_VegClass_data.et;
    modified_param_GVF_VegClass_SWE = modified_param_GVF_VegClass_data.swe;
    modified_param_GVF_VegClass_ET = diff(modified_param_GVF_VegClass_ET);
    modified_param_GVF_VegClass_edir = modified_param_GVF_VegClass_data.edir;
    modified_param_GVF_VegClass_etran = modified_param_GVF_VegClass_data.etran;
    modified_param_GVF_VegClass_ecan = modified_param_GVF_VegClass_data.ecan;

    modified_param_GVF_VegClass_edir = diff(modified_param_GVF_VegClass_edir);
    modified_param_GVF_VegClass_etran = diff(modified_param_GVF_VegClass_etran);
    modified_param_GVF_VegClass_ecan = diff(modified_param_GVF_VegClass_ecan);

    modified_param_GVF_VegClass_soilm1 = modified_param_GVF_VegClass_data.soilm1;
    modified_param_GVF_VegClass_soilm2 = modified_param_GVF_VegClass_data.soilm2;
    modified_param_GVF_VegClass_soilm3 = modified_param_GVF_VegClass_data.soilm3;
    modified_param_GVF_VegClass_soilm4 = modified_param_GVF_VegClass_data.soilm4;
    modified_param_GVF_VegClass_soilm = modified_param_GVF_VegClass_soilm1.*(0.1/2) + modified_param_GVF_VegClass_soilm2.*(0.3/2) + modified_param_GVF_VegClass_soilm4.*(0.6/2) + modified_param_GVF_VegClass_soilm4.*(1/2);
    modified_param_GVF_VegClass_baseflow = modified_param_GVF_VegClass_data.ugdrnoff;
    modified_param_GVF_VegClass_baseflow = diff(modified_param_GVF_VegClass_baseflow);
    %get annual cycles:
    datevecs = datevec(dates);
    datevecs(1,:)=[];
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    modified_param_GVF_VegClass_ET(idx)=[];
    modified_param_GVF_VegClass_edir(idx)=[];
    modified_param_GVF_VegClass_etran(idx)=[];
    modified_param_GVF_VegClass_ecan(idx)=[];
    modified_param_GVF_VegClass_baseflow(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_ET_modified_param_GVF_VegClass = accumarray(j,modified_param_GVF_VegClass_ET,[],@nanmean);
    annual_ET_modified_param_GVF_VegClass = [annual_ET_modified_param_GVF_VegClass(end);annual_ET_modified_param_GVF_VegClass(1:end-1)];

    annual_edir_modified_param_GVF_VegClass = accumarray(j,modified_param_GVF_VegClass_edir,[],@nanmean);
    annual_edir_modified_param_GVF_VegClass = [annual_edir_modified_param_GVF_VegClass(end);annual_edir_modified_param_GVF_VegClass(1:end-1)];

    annual_etran_modified_param_GVF_VegClass = accumarray(j,modified_param_GVF_VegClass_etran,[],@nanmean);
    annual_etran_modified_param_GVF_VegClass = [annual_etran_modified_param_GVF_VegClass(end);annual_etran_modified_param_GVF_VegClass(1:end-1)];

    annual_ecan_modified_param_GVF_VegClass = accumarray(j,modified_param_GVF_VegClass_ecan,[],@nanmean);
    annual_ecan_modified_param_GVF_VegClass = [annual_ecan_modified_param_GVF_VegClass(end);annual_ecan_modified_param_GVF_VegClass(1:end-1)];

    annual_baseflow_modified_param_GVF_VegClass = accumarray(j,modified_param_GVF_VegClass_baseflow,[],@nanmean);
    annual_baseflow_modified_param_GVF_VegClass = [annual_baseflow_modified_param_GVF_VegClass(end);annual_baseflow_modified_param_GVF_VegClass(1:end-1)];

    datevecs = datevec(dates);
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    modified_param_GVF_VegClass_SWE(idx)=[];
    modified_param_GVF_VegClass_soilm(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_SWE_modified_param_GVF_VegClass = accumarray(j,modified_param_GVF_VegClass_SWE,[],@nanmean);
    annual_SM_modified_param_GVF_VegClass = accumarray(j,modified_param_GVF_VegClass_soilm,[],@nanmean);

    %realistic data
    realistic_data = load(realistic_filename);
    realistic_data = realistic_data.Catchment_outputs;
    realistic_ET = realistic_data.et;
    realistic_SWE = realistic_data.swe;
    realistic_ET = diff(realistic_ET);
    realistic_edir = realistic_data.edir;
    realistic_etran = realistic_data.etran;
    realistic_ecan = realistic_data.ecan;

    realistic_edir = diff(realistic_edir);
    realistic_etran = diff(realistic_etran);
    realistic_ecan = diff(realistic_ecan);

    realistic_soilm1 = realistic_data.soilm1;
    realistic_soilm2 = realistic_data.soilm2;
    realistic_soilm3 = realistic_data.soilm3;
    realistic_soilm4 = realistic_data.soilm4;
    realistic_soilm = realistic_soilm1.*(0.1/2) + realistic_soilm2.*(0.3/2) + realistic_soilm4.*(0.6/2) + realistic_soilm4.*(1/2);
    realistic_baseflow = realistic_data.ugdrnoff;
    realistic_baseflow = diff(realistic_baseflow);
    %get annual cycles:
    datevecs = datevec(dates);
    datevecs(1,:)=[];
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    realistic_ET(idx)=[];
    realistic_edir(idx)=[];
    realistic_etran(idx)=[];
    realistic_ecan(idx)=[];
    realistic_baseflow(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_ET_realistic = accumarray(j,realistic_ET,[],@nanmean);
    annual_ET_realistic = [annual_ET_realistic(end);annual_ET_realistic(1:end-1)];

    annual_edir_realistic = accumarray(j,realistic_edir,[],@nanmean);
    annual_edir_realistic = [annual_edir_realistic(end);annual_edir_realistic(1:end-1)];

    annual_etran_realistic = accumarray(j,realistic_etran,[],@nanmean);
    annual_etran_realistic = [annual_etran_realistic(end);annual_etran_realistic(1:end-1)];

    annual_ecan_realistic = accumarray(j,realistic_ecan,[],@nanmean);
    annual_ecan_realistic = [annual_ecan_realistic(end);annual_ecan_realistic(1:end-1)];

    annual_baseflow_realistic = accumarray(j,realistic_baseflow,[],@nanmean);
    annual_baseflow_realistic = [annual_baseflow_realistic(end);annual_baseflow_realistic(1:end-1)];
    datevecs = datevec(dates);
    idx=find(datevecs(:,2) == 2 & datevecs(:,3)==29);
    datevecs(idx,:)=[];
    realistic_SWE(idx)=[];
    realistic_soilm(idx)=[];
    [u,~,j] = unique(datevecs(:,2:3),'rows','stable');
    annual_SWE_realistic = accumarray(j,realistic_SWE,[],@nanmean);
    annual_SM_realistic = accumarray(j,realistic_soilm,[],@nanmean);

    if c==1
        f=figure;
        subplot(2,1,1)
        hold on
        plot_x=1:365;
        p01 = plot(plot_x,annual_ET_baseline,'-','linewidth',3.5,'color','k');
        p02 = plot(plot_x, annual_ET_mod_param,'-','linewidth',1,'color','b');
        p03 = plot(plot_x, annual_ET_modified_param_GVF,'-','linewidth',1,'color',[0 0.5 0]);
        p04 = plot(plot_x, annual_ET_modified_param_GVF_VegClass,'--','linewidth',1,'color',[0 0.5 0]);
        p05 = plot(plot_x, annual_ET_realistic,'-','linewidth',1,'color','r');
        grid on
        box on

        set(gca,'fontsize',44)
        ylabel('ET (mm/day)','fontsize',44)
        xlim([min(plot_x) max(plot_x)])
        

        f.Position = [-1919        -299        1165        1096];

        subplot(2,1,2)
        hold on
        p01 = plot(-1.*plot_x, annual_ET_mod_param-annual_ET_baseline,'-','linewidth',1,'color','k'); %for legend
        p02 = plot(plot_x, annual_ET_mod_param-annual_ET_baseline,'-','linewidth',1,'color','b');
        p03 = plot(plot_x, annual_ET_modified_param_GVF-annual_ET_baseline,'-','linewidth',1,'color',[0 0.5 0]);
        p04 = plot(plot_x, annual_ET_modified_param_GVF_VegClass-annual_ET_baseline,'--','linewidth',1,'color',[0 0.5 0]);
        p05 = plot(plot_x, annual_ET_realistic-annual_ET_baseline,'-','linewidth',1,'color','r');
        plot(plot_x,linspace(0,0,length(plot_x)),'--k','linewidth',2)
        grid on
        box on
        set(gca,'fontsize',44)
        ylabel({'\DeltaET relative to'; 'baseline (mm/day)'},'fontsize',44)
        xlim([min(plot_x) max(plot_x)])
        xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);

        leg = legend([p01 p02 p03 p04 p05],{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'},'fontsize',44,'location','best');
        leg.Position = [0.1202    0.1063    0.5313    0.1886];

        saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/AnnCycle_LEGENDS.png'))
    end
    
    f=figure;
    subplot(2,1,1)
    hold on
    plot_x=1:365;
    p01 = plot(plot_x,annual_ET_baseline,'-','linewidth',3.5,'color','k');
    p02 = plot(plot_x, annual_ET_mod_param,'-','linewidth',1,'color','b');
    p03 = plot(plot_x, annual_ET_modified_param_GVF,'-','linewidth',1,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_ET_modified_param_GVF_VegClass,'--','linewidth',1,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_ET_realistic,'-','linewidth',1,'color','r');
    grid on
    box on

    set(gca,'fontsize',44)
    ylabel('ET (mm/day)','fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);

    f.Position = [-1919        -299        1165        1096];

    subplot(2,1,2)
    hold on
    p01 = plot(-1.*plot_x, annual_ET_mod_param-annual_ET_baseline,'-','linewidth',1,'color','k'); %for legend
    p02 = plot(plot_x, annual_ET_mod_param-annual_ET_baseline,'-','linewidth',1,'color','b');
    p03 = plot(plot_x, annual_ET_modified_param_GVF-annual_ET_baseline,'-','linewidth',1,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_ET_modified_param_GVF_VegClass-annual_ET_baseline,'--','linewidth',1,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_ET_realistic-annual_ET_baseline,'-','linewidth',1,'color','r');
    plot(plot_x,linspace(0,0,length(plot_x)),'--k','linewidth',2)
    grid on
    box on
    set(gca,'fontsize',44)
    ylabel({'\DeltaET relative to'; 'baseline (mm/day)'},'fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
%     if c==1
%         leg = legend([p01 p02 p03 p04 p05],{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'},'fontsize',44,'location','best');
%         leg.Position = [0.1202    0.1063    0.5313    0.1886];
%     end
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/ET_ensemble_AnnCycle_%s.png',catch_names{c}))

    if c==1
        f=figure;
        subplot(2,1,1)
        hold on
        plot_x=1:365;
        p01 = plot(plot_x,annual_etran_baseline,'-','linewidth',1.5,'color',[0 0.5 0]);
        p02 = plot(plot_x,annual_ecan_baseline,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
        p03 = plot(plot_x,annual_edir_baseline,'-','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);

        p011 = plot(plot_x,annual_etran_realistic,'--','linewidth',1.5,'color',[0 0.5 0]);
        p022 = plot(plot_x,annual_ecan_realistic,'--','linewidth',1.5,'color',[0.5 0.5 0.5]);
        p033 = plot(plot_x,annual_edir_realistic,'--','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);

        grid on
        box on

        set(gca,'fontsize',44)
        ylabel('ET partition (mm/day)','fontsize',44)
        xlim([min(plot_x) max(plot_x)])
        xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);

        leg = legend([p01 p02 p03 p011 p022 p033],{'Baseline E_{tran}','Baseline E_{can}','Baseline E_{dir}','Mod-params+GVF+Veg-class+snow-alb E_{tran}','Mod-params+GVF+Veg-class+snow-alb E_{can}','Mod-params+GVF+Veg-class+snow-alb E_{dir}'},'fontsize',44,'location','best');
        f.Position = [-1919        -299        1165        1096];
        saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/ET_partition_ensemble_FORLEGEND.png'))
    end


    f=figure;
    subplot(2,1,1)
    hold on
    plot_x=1:365;
    p01 = plot(plot_x,annual_etran_baseline,'-','linewidth',1.5,'color',[0 0.5 0]);
    p02 = plot(plot_x,annual_ecan_baseline,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p03 = plot(plot_x,annual_edir_baseline,'-','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);

    p011 = plot(plot_x,annual_etran_realistic,'--','linewidth',1.5,'color',[0 0.5 0]);
    p022 = plot(plot_x,annual_ecan_realistic,'--','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p033 = plot(plot_x,annual_edir_realistic,'--','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);

    grid on
    box on

    set(gca,'fontsize',44)
    ylabel('ET partition (mm/day)','fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
    %     if c==1
    %         leg = legend([p01 p02 p03 p011 p022 p033],{'Baseline E_{tran}','Baseline E_{can}','Baseline E_{dir}','Mod-params+GVF+Veg-class+snow-alb E_{tran}','Mod-params+GVF+Veg-class+snow-alb E_{can}','Mod-params+GVF+Veg-class+snow-alb E_{dir}'},'fontsize',44,'location','best');
    %     end
    f.Position = [-1919        -299        1165        1096];

    subplot(2,1,2)
    hold on
    p01 = plot(plot_x,annual_etran_realistic-annual_etran_baseline,'-','linewidth',1.5,'color',[0 0.5 0]);
    p02 = plot(plot_x,annual_ecan_realistic-annual_ecan_baseline,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p03 = plot(plot_x,annual_edir_realistic-annual_edir_baseline,'-','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);
    plot(plot_x,linspace(0,0,length(plot_x)),'--k','linewidth',2)
    grid on
    box on
    set(gca,'fontsize',44)
    ylabel({'\DeltaET partition relative to'; 'baseline (mm/day)'},'fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);

    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/ET_partition_ensemble_AnnCycle_%s.png',catch_names{c}))
    %report when Etran outweighs EDIR:
    delta_edir = annual_edir_realistic-annual_edir_baseline;
    delta_etran = annual_etran_realistic-annual_etran_baseline;
    idx=find(abs(delta_etran) > abs(delta_edir));
    id = find(idx>175);
    ID = idx(id(1));
    ex_dates =datenum([2009 10 1]):datenum([2010 9 30]);
    ex_datevecs = datevec(ex_dates);
    sprintf('Etran outweighs Edir on %d/%d/%d at %s',ex_datevecs(ID,1:3),catch_names{c})

    f=figure;
    subplot(2,1,1)
    hold on
    plot_x=1:365;
    p01 = plot(plot_x,annual_etran_modified_param_GVF,'-','linewidth',1.5,'color',[0 0.5 0]);
    p02 = plot(plot_x,annual_ecan_modified_param_GVF,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p03 = plot(plot_x,annual_edir_modified_param_GVF,'-','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);

    p011 = plot(plot_x,annual_etran_modified_param_GVF_VegClass,'--','linewidth',1.5,'color',[0 0.5 0]);
    p022 = plot(plot_x,annual_ecan_modified_param_GVF_VegClass,'--','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p033 = plot(plot_x,annual_edir_modified_param_GVF_VegClass,'--','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);

    grid on
    box on

    set(gca,'fontsize',44)
    ylabel('ET partition (mm/day)','fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
    %     if c==1
    %         leg = legend([p01 p02 p03 p011 p022 p033],{'mod-params+GVF E_{tran}','mod-params+GVF E_{can}','mod-params+GVF E_{dir}','Mod-params+GVF+Veg-class E_{tran}','Mod-params+GVF+Veg-class E_{can}','Mod-params+GVF+Veg-class E_{dir}'},'fontsize',44,'location','best');
    %         leg.Position = [0.1371    0.7680    0.3562    0.1784];
    %     end
    f.Position = [-1919        -299        1165        1096];

    subplot(2,1,2)
    hold on
    p01 = plot(plot_x,annual_etran_modified_param_GVF_VegClass-annual_etran_modified_param_GVF,'-','linewidth',1.5,'color',[0 0.5 0]);
    p02 = plot(plot_x,annual_ecan_modified_param_GVF_VegClass-annual_ecan_modified_param_GVF,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p03 = plot(plot_x,annual_edir_modified_param_GVF_VegClass-annual_edir_modified_param_GVF,'-','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);
    plot(plot_x,linspace(0,0,length(plot_x)),'--k','linewidth',2)
    grid on
    box on
    set(gca,'fontsize',44)
    ylabel({'\DeltaET partition';'Mod-params+GVF+Veg-class'; '- mod-params+GVF'; '(mm/day)'},'fontsize',20)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
    if c==2
    ylim([-0.0225 0.0425])
    end
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/ET_partition_ensemble_AnnCycle_%s_VegClassImpact.png',catch_names{c}))

    f=figure;
    subplot(2,1,1)
    hold on
    plot_x=1:365;
    p01 = plot(plot_x,annual_etran_modified_param_GVF_VegClass,'-','linewidth',1.5,'color',[0 0.5 0]);
    p02 = plot(plot_x,annual_ecan_modified_param_GVF_VegClass,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p03 = plot(plot_x,annual_edir_modified_param_GVF_VegClass,'-','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);

    p011 = plot(plot_x,annual_etran_realistic,'--','linewidth',1.5,'color',[0 0.5 0]);
    p022 = plot(plot_x,annual_ecan_realistic,'--','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p033 = plot(plot_x,annual_edir_realistic,'--','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);

    grid on
    box on

    set(gca,'fontsize',44)
    ylabel('ET partition (mm/day)','fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
    % %     if c==1
    % %         leg = legend([p01 p02 p03 p011 p022 p033],{'mod-params+GVF+Veg-class E_{tran}','mod-params+GVF+Veg-class E_{can}','mod-params+GVF+Veg-class E_{dir}','Mod-params+GVF+Veg-class+snow-alb E_{tran}','Mod-params+GVF+Veg-class+snow-alb E_{can}','Mod-params+GVF+Veg-class+snow-alb E_{dir}'},'fontsize',23,'location','best');
    % %         leg.Position = [0.1371    0.7680    0.3562    0.1784];
    % %     end
    f.Position = [-1919        -299        1165        1096];

    subplot(2,1,2)
    hold on
    p01 = plot(plot_x,annual_etran_realistic-annual_etran_modified_param_GVF_VegClass,'-','linewidth',1.5,'color',[0 0.5 0]);
    p02 = plot(plot_x,annual_ecan_realistic-annual_ecan_modified_param_GVF_VegClass,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
    p03 = plot(plot_x,annual_edir_realistic-annual_edir_modified_param_GVF_VegClass,'-','linewidth',1.5,'color',[0.8500, 0.3250, 0.0980]);
    plot(plot_x,linspace(0,0,length(plot_x)),'--k','linewidth',2)
    grid on
    box on
    set(gca,'fontsize',44)
    ylabel({'\DeltaET partition';'Mod-params+GVF+Veg-class+snow-alb -'; 'mod-params+GVF+Veg-class'; '(mm/day)'},'fontsize',20)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/ET_partition_ensemble_AnnCycle_%s_AlbedoImpact.png',catch_names{c}))

    f=figure;
    subplot(2,1,1)
    hold on
    p01 = plot(plot_x,annual_SWE_baseline,'-','linewidth',2,'color','k');
    p02 = plot(plot_x, annual_SWE_mod_param,'-','linewidth',1.5,'color','b');
    p03 = plot(plot_x, annual_SWE_modified_param_GVF,'-','linewidth',1.5,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_SWE_modified_param_GVF_VegClass,'--','linewidth',1.5,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_SWE_realistic,'-','linewidth',1.5,'color','r');
    grid on
    box on

    set(gca,'fontsize',44)
    ylabel('SWE (mm)','fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
    f.Position = [-1919        -299        1165        1096];

    subplot(2,1,2)
    hold on
    p02 = plot(plot_x, annual_SWE_mod_param - annual_SWE_baseline,'-','linewidth',1.5,'color','b');
    p03 = plot(plot_x, annual_SWE_modified_param_GVF - annual_SWE_baseline,'-','linewidth',1.5,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_SWE_modified_param_GVF_VegClass - annual_SWE_baseline,'--','linewidth',1.5,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_SWE_realistic - annual_SWE_baseline,'-','linewidth',1.5,'color','r');
    grid on
    box on

    set(gca,'fontsize',44)
    ylabel({'\Delta SWE';'relative to baseline (mm)'},'fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/SWE_ensemble_AnnCycle_%s.png',catch_names{c}))

    % report stats for paper:
    delta_peak_SWE_GVF = max(annual_SWE_modified_param_GVF) - max(annual_SWE_baseline);
    delta_peak_SWE_GVF_pct = delta_peak_SWE_GVF/max(annual_SWE_baseline) * 100;
    sprintf('GVF reduction causes %.4f%% (%.4fmm) peak SWE change at %s',delta_peak_SWE_GVF_pct,delta_peak_SWE_GVF,catch_names{c})

    delta_peak_SWE_class = max(annual_SWE_modified_param_GVF_VegClass) - max(annual_SWE_modified_param_GVF);
    delta_peak_SWE_class_pct = delta_peak_SWE_class/max(annual_SWE_modified_param_GVF) * 100;
    sprintf('veg class reduction causes %.4f%% (%.4fmm) peak SWE change at %s',delta_peak_SWE_class_pct,delta_peak_SWE_class,catch_names{c})


    delta_peak_SWE_real= max(annual_SWE_realistic) - max(annual_SWE_modified_param_GVF_VegClass);
    delta_peak_SWE_real_pct = delta_peak_SWE_real/max(annual_SWE_modified_param_GVF_VegClass) * 100;
    sprintf('real reduction causes %.4f%% (%.4fmm) peak SWE change at %s',delta_peak_SWE_real_pct,delta_peak_SWE_real,catch_names{c})

    delta_peak_SWE_real= max(annual_SWE_realistic) - max(annual_SWE_baseline);
    delta_peak_SWE_real_pct = delta_peak_SWE_real/max(annual_SWE_baseline) * 100;
    sprintf('real reduction causes %.4f%% (%.4fmm) peak SWE change at %s relative to baseline',delta_peak_SWE_real_pct,delta_peak_SWE_real,catch_names{c})

    %get ablation rates:
    [baseline_peak,i] = max(annual_SWE_baseline);
    idx_80 = find(annual_SWE_baseline < 0.8* baseline_peak);
    ID=find(idx_80 > i);
    start_of_ablation = idx_80(ID(1));

    idx_20 = find(annual_SWE_baseline < 0.2* baseline_peak);
    ID=find(idx_20 > i);
    end_of_ablation = idx_20(ID(1));
    baseline_ablation = (annual_SWE_baseline(start_of_ablation) - annual_SWE_baseline(end_of_ablation))/length([start_of_ablation:end_of_ablation]);
    sprintf('baseline ablation rate - %.4f mm/day for %s',baseline_ablation,catch_names{c})
    ex_dates = datenum([2009 10 1]):datenum([2010 9 30]);
    ex_datevecs = datevec(ex_dates);
    sprintf('ablation rate start: %04d/%02d/%02d',ex_datevecs(start_of_ablation,1),ex_datevecs(start_of_ablation,2),ex_datevecs(start_of_ablation,3))
    sprintf('ablation rate end: %04d/%02d/%02d',ex_datevecs(end_of_ablation,1),ex_datevecs(end_of_ablation,2),ex_datevecs(end_of_ablation,3))

    [modified_param_GVF_peak,i] = max(annual_SWE_modified_param_GVF);
    idx_80 = find(annual_SWE_modified_param_GVF < 0.8* modified_param_GVF_peak);
    ID=find(idx_80 > i);
    start_of_ablation = idx_80(ID(1));

    idx_20 = find(annual_SWE_modified_param_GVF < 0.2* modified_param_GVF_peak);
    ID=find(idx_20 > i);
    end_of_ablation = idx_20(ID(1));
    modified_param_GVF_ablation = (annual_SWE_modified_param_GVF(start_of_ablation) - annual_SWE_modified_param_GVF(end_of_ablation))/length([start_of_ablation:end_of_ablation]);
    sprintf('modified_param_GVF ablation rate - %.4f mm/day for %s',modified_param_GVF_ablation,catch_names{c})

    [modified_param_GVF_VegClass_peak,i] = max(annual_SWE_modified_param_GVF_VegClass);
    idx_80 = find(annual_SWE_modified_param_GVF_VegClass < 0.8* modified_param_GVF_VegClass_peak);
    ID=find(idx_80 > i);
    start_of_ablation = idx_80(ID(1));

    idx_20 = find(annual_SWE_modified_param_GVF_VegClass < 0.2* modified_param_GVF_VegClass_peak);
    ID=find(idx_20 > i);
    end_of_ablation = idx_20(ID(1));
    modified_param_GVF_VegClass_ablation = (annual_SWE_modified_param_GVF_VegClass(start_of_ablation) - annual_SWE_modified_param_GVF_VegClass(end_of_ablation))/length([start_of_ablation:end_of_ablation]);
    sprintf('modified_param_GVF_VegClass ablation rate - %.4f mm/day for %s',modified_param_GVF_VegClass_ablation,catch_names{c})

    [realistic_peak,i] = max(annual_SWE_realistic);
    idx_80 = find(annual_SWE_realistic < 0.8* realistic_peak);
    ID=find(idx_80 > i);
    start_of_ablation = idx_80(ID(1));

    idx_20 = find(annual_SWE_realistic < 0.2* realistic_peak);
    ID=find(idx_20 > i);
    end_of_ablation = idx_20(ID(1));
    realistic_ablation = (annual_SWE_realistic(start_of_ablation) - annual_SWE_realistic(end_of_ablation))/length([start_of_ablation:end_of_ablation]);
    sprintf('realistic ablation rate - %.4f mm/day for %s',realistic_ablation,catch_names{c})

    p01 = plot(plot_x,annual_SWE_baseline,'-','linewidth',2,'color','k');
    p02 = plot(plot_x, annual_SWE_mod_param,'-','linewidth',1.5,'color','b');
    p03 = plot(plot_x, annual_SWE_modified_param_GVF,'-','linewidth',1.5,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_SWE_modified_param_GVF_VegClass,'--','linewidth',1.5,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_SWE_realistic,'-','linewidth',1.5,'color','r');

    f=figure;
    subplot(2,1,1)
    hold on
    p01 = plot(plot_x,annual_SM_baseline,'-','linewidth',3.5,'color','k');
    p02 = plot(plot_x, annual_SM_mod_param,'-','linewidth',1,'color','b');
    p03 = plot(plot_x, annual_SM_modified_param_GVF,'-','linewidth',1,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_SM_modified_param_GVF_VegClass,'--','linewidth',1,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_SM_realistic,'-','linewidth',3,'color','r');
    grid on
    box on
    set(gca,'fontsize',44)
    ylabel({'total column'; 'SM (mm/mm)'},'fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);


    subplot(2,1,2)
    hold on
    p01 = plot(plot_x, annual_SM_mod_param - annual_SM_baseline,'-','linewidth',1.5,'color','k');
    p02 = plot(plot_x, annual_SM_mod_param - annual_SM_baseline,'-','linewidth',1.5,'color','b');
    p03 = plot(plot_x, annual_SM_modified_param_GVF - annual_SM_baseline,'-','linewidth',1.5,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_SM_modified_param_GVF_VegClass - annual_SM_baseline,'--','linewidth',1.5,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_SM_realistic - annual_SM_baseline,'-','linewidth',1.5,'color','r');
    grid on
    box on

    set(gca,'fontsize',44)
    ylabel({'\Delta SM relative to'; 'baseline (mm/mm)'},'fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
% %     if c==1
% %         leg = legend([p01 p02 p03 p04 p05],{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'},'fontsize',44,'location','best');
% %         leg.Position = [0.2927    0.4114    0.5313    0.1886];
% %     end
    f.Position = [-1919        -299        1165        1096];
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/SM_ensemble_AnnCycle_%s.png',catch_names{c}))
    ex_dates = datenum([2005 10 1]):datenum([2006 9 30]);
    ex_datevecs =datevec(ex_dates);

    idx_fall = find(ex_datevecs(:,2)==9 | ex_datevecs(:,2)==10 | ex_datevecs(:,2)==11);
    idx_winter = find(ex_datevecs(:,2)==12 | ex_datevecs(:,2)==1 | ex_datevecs(:,2)==2);
    idx_spring = find(ex_datevecs(:,2)==3 | ex_datevecs(:,2)==4 | ex_datevecs(:,2)==5);
    idx_summer = find(ex_datevecs(:,2)==6 | ex_datevecs(:,2)==7 | ex_datevecs(:,2)==8);

    sprintf('mod-params results in %.4f%% change in SM at %s',mean((annual_SM_mod_param - annual_SM_baseline)./annual_SM_baseline) * 100,catch_names{c})
    sprintf('mod-params+gvf results in %.4f%% change in SM at %s',mean((annual_SM_modified_param_GVF - annual_SM_baseline)./annual_SM_baseline) * 100,catch_names{c})
    sprintf('mod-params+gvf results in %.4f%% change in fall SM at %s',mean((annual_SM_modified_param_GVF(idx_fall) - annual_SM_baseline(idx_fall))./annual_SM_baseline(idx_fall)) * 100,catch_names{c})
    sprintf('mod-params+gvf results in %.4f%% change in winter SM at %s',mean((annual_SM_modified_param_GVF(idx_winter) - annual_SM_baseline(idx_winter))./annual_SM_baseline(idx_winter)) * 100,catch_names{c})
    sprintf('mod-params+gvf results in %.4f%% change in spring SM at %s',mean((annual_SM_modified_param_GVF(idx_spring) - annual_SM_baseline(idx_spring))./annual_SM_baseline(idx_spring)) * 100,catch_names{c})
    sprintf('mod-params+gvf results in %.4f%% change in summer SM at %s',mean((annual_SM_modified_param_GVF(idx_summer) - annual_SM_baseline(idx_summer))./annual_SM_baseline(idx_summer)) * 100,catch_names{c})
    sprintf('mod-params+gvf+class+alb results in %.4f%% change in SM at %s relative to baseline',mean((annual_SM_realistic - annual_SM_baseline)./annual_SM_baseline) * 100,catch_names{c})
    f=figure;
    subplot(2,1,1)
    hold on
    p01 = plot(plot_x,annual_baseflow_baseline,'-','linewidth',3,'color','k');
    p02 = plot(plot_x, annual_baseflow_mod_param,'-','linewidth',1,'color','b');
    p03 = plot(plot_x, annual_baseflow_modified_param_GVF,'-','linewidth',1,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_baseflow_modified_param_GVF_VegClass,'--','linewidth',1,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_baseflow_realistic,'-','linewidth',1.5,'color','r');
    grid on
    box on
    set(gca,'fontsize',44)
    ylabel({'subsurface flow'; '(mm/day)'},'fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);


    subplot(2,1,2)
    hold on
    p01 = plot(plot_x, annual_baseflow_mod_param - annual_baseflow_baseline,'-','linewidth',1.5,'color','k');
    p02 = plot(plot_x, annual_baseflow_mod_param - annual_baseflow_baseline,'-','linewidth',1.5,'color','b');
    p03 = plot(plot_x, annual_baseflow_modified_param_GVF - annual_baseflow_baseline,'-','linewidth',1.5,'color',[0 0.5 0]);
    p04 = plot(plot_x, annual_baseflow_modified_param_GVF_VegClass - annual_baseflow_baseline,'--','linewidth',1.5,'color',[0 0.5 0]);
    p05 = plot(plot_x, annual_baseflow_realistic - annual_baseflow_baseline,'-','linewidth',1.5,'color','r');
    grid on
    box on

    set(gca,'fontsize',44)
    ylabel({'\Delta subsurface flow';'relative to baseline (mm)'},'fontsize',44)
    xlim([min(plot_x) max(plot_x)])
    xlabel('day of water year','fontsize',44);xticks([1:50:365]);xtickangle(90);
% %     if c==1
% %         leg = legend([p01 p02 p03 p04 p05],{'Baseline','Mod-params','Mod-params+GVF','Mod-params+GVF+Veg-class','Mod-params+GVF+Veg-class+snow-alb'},'fontsize',44,'location','best');
% %         leg.Position = [0.3869    0.3588    0.5313    0.1886];
% %     end
    f.Position = [-1919        -299        1165        1096];
    saveas(f,sprintf('/Users/abolafia/ASO_Fire/Plots/PaperPlots/baseflow_ensemble_AnnCycle_%s.png',catch_names{c}))

    delta_params = (annual_baseflow_mod_param - annual_baseflow_baseline)./annual_baseflow_baseline;
    delta_params_gvf = (annual_baseflow_modified_param_GVF - annual_baseflow_baseline)./annual_baseflow_baseline;
    delta_params_gvf_class = (annual_baseflow_modified_param_GVF_VegClass - annual_baseflow_baseline)./annual_baseflow_baseline;
    delta_params_gvf_class_albedo = (annual_baseflow_realistic - annual_baseflow_baseline)./annual_baseflow_baseline;
    delta_params_class = (delta_params_gvf_class - delta_params_gvf)./delta_params_gvf;

    sprintf('mean baseflow change from params = %.4f at %s', nanmean(delta_params) * 100 , catch_names{c})
    sprintf('mean baseflow change from params+gvf = %.4f at %s', nanmean(delta_params_gvf) * 100 , catch_names{c})
    sprintf('mean baseflow change from class = %.4f at %s', nanmean(annual_baseflow_modified_param_GVF_VegClass - annual_baseflow_modified_param_GVF)/nanmean(annual_baseflow_modified_param_GVF) * 100 , catch_names{c})
    sprintf('mean baseflow change from params+gvf+class+alb = %.4f at %s', nanmean(delta_params_gvf_class_albedo) * 100 , catch_names{c})

    idx_jan_mar = 93:182;
    sprintf('mean baseflow change from class in Jan-Mar = %.4f at %s', nanmean(annual_baseflow_modified_param_GVF_VegClass(idx_jan_mar) - annual_baseflow_modified_param_GVF(idx_jan_mar))/nanmean(annual_baseflow_modified_param_GVF(idx_jan_mar)) * 100 , catch_names{c})

    idx_apr_jun = 183:273;
    sprintf('mean baseflow change from class in Apr-Jun = %.4f at %s', nanmean(annual_baseflow_modified_param_GVF_VegClass(idx_apr_jun) - annual_baseflow_modified_param_GVF(idx_apr_jun))/nanmean(annual_baseflow_modified_param_GVF(idx_apr_jun)) * 100 , catch_names{c})

    idx_midApr_midJune = 197:258;
    sprintf('mean baseflow change from params+gvf+class+alb From apr 15 - Jun 15 = %.4f at %s', nanmean(delta_params_gvf_class_albedo(idx_midApr_midJune) - annual_baseflow_baseline(idx_midApr_midJune))/nanmean(annual_baseflow_baseline(idx_midApr_midJune)) * 100 , catch_names{c})

    %%report stats to go into paper:
    delta_ET_real = realistic_ET - baseline_ET;
    delta_ET_params = mod_param_ET - baseline_ET;
    ET_dates = datenum(datevecs);
    ET_dates(1)=[];
    store_delta_ET=[];
    store_delta_ET_params=[];
    %get range of ET differences in each year:
    for WY=2000:2022
        idx_WY = find(ET_dates>=datenum([WY-1 10 1]) & ET_dates<=datenum([WY 9 30]));
        current_delta = delta_ET_real(idx_WY);
        current_delta_params = delta_ET_params(idx_WY);
        current_baseline =baseline_ET(idx_WY);
        store_delta_ET = [store_delta_ET;nanmean(current_delta)/nanmean(current_baseline)*100];
        store_delta_ET_params = [store_delta_ET_params;nanmean(current_delta_params)/nanmean(current_baseline)*100];
    end
    sprintf('real - baseline ET varies from %.4f-%.4f%% (mean=%.4f%%) at %s',min(store_delta_ET),max(store_delta_ET),mean(store_delta_ET),catch_names{c})
    sprintf('modparams - baseline ET varies from %.4f-%.4f%% (mean=%.4f%%) at %s',min(store_delta_ET_params),max(store_delta_ET_params),mean(store_delta_ET_params),catch_names{c})

end
