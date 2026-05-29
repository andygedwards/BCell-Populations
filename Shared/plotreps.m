function [varargout] = plotreps(allvals, params, param_names, states, state_names, model_rhsFun, protocol, savedir)
mean_val = mean(allvals);
if mean_val < 0
    min_ind = find(allvals==max(allvals),1,'first'); 
    max_ind = find(allvals==min(allvals),1,'first'); 
else
    min_ind = find(allvals==min(allvals),1,'first'); 
    max_ind = find(allvals==max(allvals),1,'first'); 
end
if mod(length(allvals),2)==0
    mindiff = min(abs(allvals-median(allvals)));
    med_ind = find(abs(allvals-median(allvals))==mindiff,1,'first');
else 
    med_ind = find(allvals==median(allvals),1,'first');
end
inds = [min_ind,med_ind,max_ind];
disp(['min =', num2str(allvals(min_ind))])
disp(['med =', num2str(allvals(med_ind))])
disp(['max =', num2str(allvals(max_ind))])
lines = {'-','--',':'};
if strcmp(protocol,'INa')
    ts = figure;
    boltz = figure;
    V = (-120:1:0);
    for i = 1:length(inds)
        line = lines{i};
        figure(ts), hold on;
        [v_half,r_square,f,T,Y] = INa_inact_plot(params(inds(i),:), param_names, states, state_names, model_rhsFun, line);
        figure(boltz), hold on;
        plot(V,feval(f,V),[line,'k']);
        v_halfs{i} = ['v_{half} = ',num2str(v_half)];
    end
    saveas(ts, [savedir, filesep, 'INa_TS.fig'], 'fig')
    exportgraphics(ts, [savedir, filesep, 'INa_TS.eps'],  'ContentType', 'vector')

    legend(boltz.CurrentAxes, v_halfs, 'Location', 'NorthEast')
    saveas(boltz, [savedir, filesep, 'INa_boltz.fig'], 'fig')
    exportgraphics(boltz, [savedir, filesep, 'INa_boltz.eps'],  'ContentType', 'vector')
    
elseif strcmp(protocol,'peak currents')
    figure
    for i = 1:length(inds)
        line = lines(i);
        [peak_INa(i), v_half_act(i), vha_r2(i), peak_ICa(i), late_ICa(i)] = PeakCurrents_plot(params(inds(i),:), param_names, states, state_names, model_rhsFun, line);
    end
    saveas(gcf, [savedir, filesep, 'Peak_currents.fig'],  'fig')
    exportgraphics(gcf, [savedir, filesep, 'Peak_currents.eps'],  'ContentType', 'vector')
elseif strcmp(protocol,'Exocytosis')
    figure
    for i = 1:length(inds)
        line = lines(i);
        [EE(i), ~, TE(i),T,Y] = Exocytosis_plot(params(inds(i),:), param_names, states, state_names, model_rhsFun, line);
    end
    saveas(gcf, [savedir, filesep, 'Exocytosis.fig'], 'fig')
    exportgraphics(gcf, [savedir, filesep, 'Exocytosis.eps'],  'ContentType', 'vector')
end

varargout = {1}; %placeholder for return type if needed
end

