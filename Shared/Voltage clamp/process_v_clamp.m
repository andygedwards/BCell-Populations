function [] = process_v_clamp(directory)
    
    cd(directory)
    files = dir(fullfile(directory, '*.mat'));
    bins = 100;
    Cm = 10;
    for i = 1:length(files)
        load(files(i).name);
        if exist('v_half_act', 'var')
            act_flag = 1;
            data(i, 1) = peak_INa;
            data(i, 2) = v_half;
            data(i, 3) = v_half_act;
            data(i, 4) = peak_ICa;
            data(i, 5) = late_ICa;
            data(i, 6) = EE;
            data(i, 7) = LE;
            data(i, 8) = TE;
            data(i, 9) = ramp_IK;
            data(i, 10) = Peak_IK_minus20;
            data(i, 11) = Peak_IK_max;
            data(i, 12) = Peak_IK_minus20/Peak_IK_max;
            data(i, 13) = v_half_act_IK;
            data(i, 14) = n_act_IK;
            data(i, 15) = vh_r2;
            data(i, 16) = vha_r2;
        else
            act_flag = 0;
            data(i, 1) = peak_INa;
            data(i, 2) = v_half;
            data(i, 3) = peak_ICa;
            data(i, 4) = late_ICa;
            data(i, 5) = EE;
            data(i, 6) = LE;
            data(i, 7) = TE;
            data(i, 8) = ramp_IK;
            data(i, 9) = Peak_IK_minus20;
            data(i, 10) = Peak_IK_max;
            data(i, 11) = Peak_IK_minus20/Peak_IK_max;
            data(i, 12) = v_half_act_IK;
            data(i, 13) = n_act_IK;
        end
    end

    disp(strcat('No. of failed peak I_{inward} analyses: ', num2str(sum(isnan(data(:,1))))));
    disp(strcat('No. of failed v_half analyses: ', num2str(sum(isnan(data(:,2))))));
    disp(strcat('No. of failed exocytosis analyses: ', num2str(sum(isnan(data(:,6))))));
    disp(strcat('No. of failed peak I_K analyses: ', num2str(sum(isnan(data(:,9))))));
    disp(strcat('No. of failed fraction A-type I_K analyses: ', num2str(sum(isnan(data(:,10))))));
    disp(strcat('No. of failed I_{K,tail} analyses: ', num2str(sum(isnan(data(:,11))))));
    
    if act_flag
        bad_INa_act = find(data(:,16) < 0.7);
        bad_INa_act2 = find(data(:,3) > -10);
        good = setdiff([1:1:length(files)],[bad_INa_act;bad_INa_act2]);
        data = data(good,:);

        [peak_INa_counts, INa_edges] = histcounts(-data(:,1), bins, 'Normalization','pdf');
        [vhalf_counts, vhalf_edges] = histcounts(data(:,2), bins, 'Normalization','pdf');
        [vhalf_act_counts, vhalf_act_edges] = histcounts(data(:,3), bins, 'Normalization','pdf');
        [peak_ICa_counts, ICa_edges] = histcounts(-data(:,4), bins, 'Normalization','pdf');
        [late_ICa_counts, ICa_late_edges] = histcounts(-data(:,5), bins, 'Normalization','pdf');
        [EE_counts, EE_edges] = histcounts(data(:,6), bins, 'Normalization','pdf');
        [LE_counts, LE_edges] = histcounts(data(:,7), bins, 'Normalization','pdf');
        [TE_counts, TE_edges] = histcounts(data(:,8), bins, 'Normalization','pdf');
        [rampIK_counts, rampIK_edges] = histcounts(data(:,9), bins, 'Normalization','pdf');
        [peakIK_counts, peakIK_edges] = histcounts(data(:,10), bins, 'Normalization','pdf');
        [peakIK_max_counts, peakIK_max_edges] = histcounts(data(:,11), bins, 'Normalization','pdf');
        [normAK_counts, normAK_edges] = histcounts(data(:,12), bins, 'Normalization','pdf');
        [vhalfIK_counts, vhalfIK_edges] = histcounts(data(:,13), bins, 'Normalization','pdf');
        [nIK_counts, nIK_edges] = histcounts(data(:,14), bins, 'Normalization','pdf');
    
        peak_INa_fit = fitdist(-data(:,1),'Normal'); 
        peak_INa_fit = pdf(peak_INa_fit,INa_edges);
        vhalf_INa_fit = fitdist(data(:,2),'Normal');
        vhalf_INa_fit = pdf(vhalf_INa_fit,vhalf_edges);
        vhalf_act_INa_fit = fitdist(data(:,3),'Normal');
        vhalf_act_INa_fit = pdf(vhalf_act_INa_fit,vhalf_edges);
        peak_ICa_fit = fitdist(-data(:,4),'Normal');
        peak_ICa_fit = pdf(peak_ICa_fit,ICa_edges);
        late_ICa_fit = fitdist(-data(:,5),'Normal');
        late_ICa_fit = pdf(late_ICa_fit,ICa_late_edges);
        EE_fit = fitdist(data(:,6),'Normal');
        EE_fit = pdf(EE_fit,EE_edges);
        LE_fit = fitdist(data(:,7),'Normal');
        LE_fit = pdf(LE_fit,LE_edges);
        TE_fit = fitdist(data(:,8),'Normal');
        TE_fit = pdf(TE_fit,TE_edges);
        rampIK_fit = fitdist(data(:,9),'Normal');
        rampIK_fit = pdf(rampIK_fit,rampIK_edges);
        peakIK_fit = fitdist(data(:,10),'Normal');
        peakIK_fit = pdf(peakIK_fit,peakIK_edges);
        peakIK_max_fit = fitdist(data(:,11),'Normal');
        peakIK_max_fit = pdf(peakIK_max_fit,peakIK_max_edges);
        normAK_fit = fitdist(data(:,12),'Normal');
        normAK_fit = pdf(normAK_fit,normAK_edges);
        vhalfIK_fit = fitdist(data(:,13),'Normal');
        vhalfIK_fit = pdf(vhalfIK_fit,vhalfIK_edges);
        nIK_fit = fitdist(data(:,14),'Normal');
        nIK_fit = pdf(nIK_fit,nIK_edges);
    else
        [peak_INa_counts, INa_edges] = histcounts(-data(:,1), bins, 'Normalization','pdf');
        [vhalf_counts, vhalf_edges] = histcounts(data(:,2), bins, 'Normalization','pdf');
        [peak_ICa_counts, ICa_edges] = histcounts(-data(:,3), bins, 'Normalization','pdf');
        [late_ICa_counts, ICa_late_edges] = histcounts(-data(:,4), bins, 'Normalization','pdf');
        [EE_counts, EE_edges] = histcounts(data(:,5), bins, 'Normalization','pdf');
        [LE_counts, LE_edges] = histcounts(data(:,6), bins, 'Normalization','pdf');
        [TE_counts, TE_edges] = histcounts(data(:,7), bins, 'Normalization','pdf');
        [rampIK_counts, rampIK_edges] = histcounts(data(:,8), bins, 'Normalization','pdf');
        [peakIK_counts, peakIK_edges] = histcounts(data(:,9), bins, 'Normalization','pdf');
        [peakIK_max_counts, peakIK_max_edges] = histcounts(data(:,10), bins, 'Normalization','pdf');
        [normAK_counts, normAK_edges] = histcounts(data(:,11), bins, 'Normalization','pdf');
        [vhalfIK_counts, vhalfIK_edges] = histcounts(data(:,12), bins, 'Normalization','pdf');
        [nIK_counts, nIK_edges] = histcounts(data(:,13), bins, 'Normalization','pdf');
    
        peak_INa_fit = fitdist(-data(:,1),'Normal'); 
        peak_INa_fit = pdf(peak_INa_fit,INa_edges);
        vhalf_INa_fit = fitdist(data(:,2),'Normal');
        vhalf_INa_fit = pdf(vhalf_INa_fit,vhalf_edges);
        peak_ICa_fit = fitdist(-data(:,3),'Normal');
        peak_ICa_fit = pdf(peak_ICa_fit,ICa_edges);
        late_ICa_fit = fitdist(-data(:,4),'Normal');
        late_ICa_fit = pdf(late_ICa_fit,ICa_late_edges);
        EE_fit = fitdist(data(:,5),'Normal');
        EE_fit = pdf(EE_fit,EE_edges);
        LE_fit = fitdist(data(:,6),'Normal');
        LE_fit = pdf(LE_fit,LE_edges);
        TE_fit = fitdist(data(:,7),'Normal');
        TE_fit = pdf(TE_fit,TE_edges);
        rampIK_fit = fitdist(data(:,8),'Normal');
        rampIK_fit = pdf(rampIK_fit,rampIK_edges);
        peakIK_fit = fitdist(data(:,9),'Normal');
        peakIK_fit = pdf(peakIK_fit,peakIK_edges);
        peakIK_max_fit = fitdist(data(:,10),'Normal');
        peakIK_max_fit = pdf(peakIK_max_fit,peakIK_max_edges);
        normAK_fit = fitdist(data(:,11),'Normal');
        normAK_fit = pdf(normAK_fit,normAK_edges);
        vhalfIK_fit = fitdist(data(:,12),'Normal');
        vhalfIK_fit = pdf(vhalfIK_fit,vhalfIK_edges);
        nIK_fit = fitdist(data(:,13),'Normal');
        nIK_fit = pdf(nIK_fit,nIK_edges);
    end

    cd ../..
    cd('Results')
    if ~isfolder(string(datetime('today')))
        mkdir(string(datetime('today')))
        cd(string(datetime('today')))
    else
        all_files = dir;
        num_dir = nnz(strfind({all_files.name},{char(datetime('today'))})&[all_files.isdir]);
        mkdir(strcat(string(datetime('today')),'_',num2str(num_dir+1)))
        cd(strcat(string(datetime('today')),'_',num2str(num_dir+1)))
    end

    f1 = figure;
    if act_flag
        subplot(3,1,1), histogram('BinEdges', INa_edges, 'BinCounts', peak_INa_counts), hold on, plot(INa_edges,peak_INa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('peak I_{Na} (pa/pF)'), xlim([0,120])
        subplot(3,1,2), histogram('BinEdges', vhalf_edges, 'BinCounts', vhalf_counts), hold on, plot(vhalf_edges,vhalf_INa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('I_{Na} inactivation V_{1/2} (mV)'), xlim([-85,-15])
        subplot(3,1,3), histogram('BinEdges', vhalf_act_edges, 'BinCounts', vhalf_act_counts), hold on, plot(vhalf_act_edges,vhalf_act_INa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('I_{Na} activation V_{1/2} (mV)'),xlim([-60,0])
        set(gcf,'position',[100,100,150,400])
        saveas(f1,strcat('INa PoM behaviour'));
        saveas(f1,strcat('INa PoM behaviour'),'png');
    else
        subplot(2,1,1), histogram('BinEdges', INa_edges, 'BinCounts', peak_INa_counts), hold on, plot(INa_edges,peak_INa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('peak I_{Na} (pa/pF)')
        subplot(2,1,2), histogram('BinEdges', vhalf_edges, 'BinCounts', vhalf_counts), hold on, plot(vhalf_edges,vhalf_INa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('I_{Na} inactivation V_{1/2} (mV)')
        set(gcf,'position',[100,100,150,400])
        saveas(f1,strcat('INa PoM (wide window) behaviour'));
        saveas(f1,strcat('INa PoM (wide window) behaviour'),'png');
    end

    f2 = figure;
    if act_flag
        subplot(2,1,1), histogram('BinEdges', ICa_edges, 'BinCounts', peak_ICa_counts), hold on, plot(ICa_edges,peak_ICa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('peak I_{Ca} (pa/pF)')
        subplot(2,1,2), histogram('BinEdges', ICa_late_edges, 'BinCounts', late_ICa_counts), hold on, plot(ICa_late_edges,late_ICa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('late I_{Ca} (pa/pF)')
        set(gcf,'position',[100,100,150,400])
        saveas(f2,strcat('ICa PoM behaviour'));
        saveas(f2,strcat('ICa PoM behaviour'),'png');
    else
        subplot(2,1,1), histogram('BinEdges', ICa_edges, 'BinCounts', peak_ICa_counts), hold on, plot(ICa_edges,peak_ICa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('peak I_{Ca} (pa/pF)')
        subplot(2,1,2), histogram('BinEdges', ICa_late_edges, 'BinCounts', late_ICa_counts), hold on, plot(ICa_late_edges,late_ICa_fit,'r','LineWidth',2), ylabel('Density'), xlabel('late I_{Ca} (pa/pF)')
        set(gcf,'position',[100,100,150,400])
        saveas(f2,strcat('ICa PoM (wide window) behaviour'));
        saveas(f2,strcat('ICa PoM (wide window) behaviour'),'png');
    end

    
    f3 = figure;
    if act_flag
        subplot(3,1,1), histogram('BinEdges', EE_edges, 'BinCounts', EE_counts), hold on, plot(EE_edges,EE_fit,'r','LineWidth',2), ylabel('Density'), xlabel('Early exocytosis (via Ca^{2+}_i, \muM)')
        subplot(3,1,2), histogram('BinEdges', LE_edges, 'BinCounts', LE_counts), hold on, plot(LE_edges,LE_fit,'r','LineWidth',2), ylabel('Density'), xlabel('Late exocytosis (via Ca^{2+}_i, \muM)')
        subplot(3,1,3), histogram('BinEdges', TE_edges, 'BinCounts', TE_counts), hold on, plot(TE_edges,TE_fit,'r','LineWidth',2), ylabel('Density'), xlabel('Total exocytosis (via Ca^{2+}_i, \muM)')
        set(gcf,'position',[100,100,150,400])
        saveas(f3,strcat('Exocytosis PoM behaviour'));
        saveas(f3,strcat('Exocytosis PoM behaviour'),'png');
    else
        subplot(3,1,1), histogram('BinEdges', EE_edges, 'BinCounts', EE_counts), hold on, plot(EE_edges,EE_fit,'r','LineWidth',2), ylabel('Density'), xlabel('Early exocytosis (via Ca^{2+}_i, \muM)')
        subplot(3,1,2), histogram('BinEdges', LE_edges, 'BinCounts', LE_counts), hold on, plot(LE_edges,LE_fit,'r','LineWidth',2), ylabel('Density'), xlabel('Late exocytosis (via Ca^{2+}_i, \muM)')
        subplot(3,1,3), histogram('BinEdges', TE_edges, 'BinCounts', TE_counts), hold on, plot(TE_edges,TE_fit,'r','LineWidth',2), ylabel('Density'), xlabel('Total exocytosis (via Ca^{2+}_i, \muM)')
        set(gcf,'position',[100,100,150,400])
        saveas(f3,strcat('Exocytosis PoM (wide window) behaviour'));
        saveas(f3,strcat('Exocytosis PoM (wide window) behaviour'),'png');
    end

    f4 = figure;
    if act_flag
        subplot(2,2,1), histogram('BinEdges', rampIK_edges, 'BinCounts', rampIK_counts), hold on, plot(rampIK_edges,rampIK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('peak total I_{K} - ramp protocol (pa/pF)')
        subplot(2,2,2), histogram('BinEdges', peakIK_edges, 'BinCounts', peakIK_counts), hold on, plot(peakIK_edges,peakIK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('total I_{K} at -20 mV (pa/pF)')
        subplot(2,2,3), histogram('BinEdges', peakIK_max_edges, 'BinCounts', peakIK_max_counts), hold on, plot(peakIK_max_edges,peakIK_max_fit,'r','LineWidth',2), ylabel('Density'), xlabel('peak total I_{K} - square pulse (pa/pF)')
        subplot(2,2,4), histogram('BinEdges', normAK_edges, 'BinCounts', normAK_counts), hold on, plot(normAK_edges,normAK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('I_{K} @ -20 mV / peak I_{K} (pa/pF)')
        set(gcf,'position',[100,100,150,400])
        saveas(f4,strcat('IK PoM behaviour'));
        saveas(f4,strcat('IK PoM behaviour'),'png');
    else
        subplot(2,2,1), histogram('BinEdges', rampIK_edges, 'BinCounts', rampIK_counts), hold on, plot(rampIK_edges,rampIK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('peak total I_{K} - ramp protocol (pa/pF)')
        subplot(2,2,2), histogram('BinEdges', peakIK_edges, 'BinCounts', peakIK_counts), hold on, plot(peakIK_edges,peakIK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('total I_{K} at -20 mV (pa/pF)')
        subplot(2,2,3), histogram('BinEdges', peakIK_max_edges, 'BinCounts', peakIK_max_counts), hold on, plot(peakIK_max_edges,peakIK_max_fit,'r','LineWidth',2), ylabel('Density'), xlabel('peak total I_{K} - square pulse (pa/pF)')
        subplot(2,2,4), histogram('BinEdges', normAK_edges, 'BinCounts', normAK_counts), hold on, plot(normAK_edges,normAK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('I_{K} @ -20 mV / peak I_{K} (pa/pF)')
        set(gcf,'position',[100,100,150,400])
        saveas(f4,strcat('IK PoM (wide window) behaviour'));
        saveas(f4,strcat('IK PoM (wide window) behaviour'),'png');
    end

    f5 = figure;
    if act_flag
        subplot(2,1,1), histogram('BinEdges', vhalfIK_edges, 'BinCounts', vhalfIK_counts), hold on, plot(vhalfIK_edges,vhalfIK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('total I_{K} V_{1/2} (mV)')
        subplot(2,1,2), histogram('BinEdges', nIK_edges, 'BinCounts', nIK_counts), hold on, plot(nIK_edges,nIK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('I_{K} activation slope factor (mV)')
        set(gcf,'position',[100,100,150,400])
        saveas(f4,strcat('IK v-dependence PoM behaviour'));
        saveas(f4,strcat('IK v-dependence PoM behaviour'),'png'); 
    else
        subplot(2,1,1), histogram('BinEdges', vhalfIK_edges, 'BinCounts', vhalfIK_counts), hold on, plot(vhalfIK_edges,vhalfIK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('total I_{K} V_{1/2} (mV)')
        subplot(2,1,2), histogram('BinEdges', nIK_edges, 'BinCounts', nIK_counts), hold on, plot(nIK_edges,nIK_fit,'r','LineWidth',2), ylabel('Density'), xlabel('I_{K} activation slope factor (mV)')
        set(gcf,'position',[100,100,150,400])
        saveas(f4,strcat('IK v-dependence (wide window) PoM behaviour'));
        saveas(f4,strcat('IK v-dependence (wide window) PoM behaviour'),'png');
    end

end
