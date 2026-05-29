function NaHisto(classes_low,classes_high,class_names,VhNa,VmNa,n_bins,savedir)
    num_classes = length(class_names);
    for c = 1:num_classes
        current_class = class_names{c};
        low_inds = find(strcmp(classes_low, current_class));
        high_inds = find(strcmp(classes_high, current_class));
        low_non_inds = find(~strcmp(classes_low, current_class));
        high_non_inds = find(~strcmp(classes_high, current_class));
        
        Vh_class_low = VhNa(low_inds);
        Vm_class_low = VmNa(low_inds);
        Vh_class_high = VhNa(high_inds);
        Vm_class_high = VmNa(high_inds);

        Vh_non_class_low = VhNa(low_non_inds);
        Vm_non_class_low = VmNa(low_non_inds);
        Vh_non_class_high = VhNa(high_non_inds);
        Vm_non_class_high = VmNa(high_non_inds);

            % Colormap
        colors = {[125, 62, 119] / 255 ...  % violet
            [85, 138, 163] / 255 ...    % blue
            [246, 88, 95] / 255 ...  % red
            }; 

        if strcmp(current_class,'Bursting')
            nbins = n_bins(2);
        else
            nbins = n_bins(1);
        end

        figure
        sgtitle([current_class,' (Low glucose)'],'FontSize',24)
        
        % Compute common bin edges for the first subplot
        min_edge1 = min([Vh_class_low; Vh_non_class_low]);
        max_edge1 = max([Vh_class_low; Vh_non_class_low]);
        bin_edges1 = linspace(min_edge1, max_edge1, nbins+1); % Ensure same bin edges
        
        subplot(2,1,1)
        histogram(Vh_class_low, 'BinEdges', bin_edges1, 'FaceColor', colors{1}, 'FaceAlpha', 0.6, 'Normalization', 'probability'), hold on
        histogram(Vh_non_class_low, 'BinEdges', bin_edges1, 'FaceColor', colors{2}, 'FaceAlpha', 0.6, 'Normalization', 'probability')
        xlabel('V_{hNa} (mV)')
        
        % Compute common bin edges for the second subplot
        min_edge2 = min([Vm_class_low; Vm_non_class_low]);
        max_edge2 = max([Vm_class_low; Vm_non_class_low]);
        bin_edges2 = linspace(min_edge2, max_edge2, nbins+1); % Ensure same bin edges
        
        subplot(2,1,2)
        histogram(Vm_class_low, 'BinEdges', bin_edges2, 'FaceColor', colors{1}, 'FaceAlpha', 0.6, 'Normalization', 'probability'), hold on
        histogram(Vm_non_class_low, 'BinEdges', bin_edges2, 'FaceColor', colors{2}, 'FaceAlpha', 0.6, 'Normalization', 'probability')
        xlabel('V_{mNa} (mV)')
        exportgraphics(gcf,[savedir,filesep,current_class,'_low.pdf'], 'ContentType', 'vector')
    
        figure
        sgtitle([current_class, ' (High glucose)'], 'FontSize', 24)
        
        % Compute common bin edges for the first subplot
        min_edge1 = min([Vh_class_high; Vh_non_class_high]);
        max_edge1 = max([Vh_class_high; Vh_non_class_high]);
        bin_edges1 = linspace(min_edge1, max_edge1, nbins+1); % Ensure same bin edges
        
        subplot(2,1,1)
        histogram(Vh_class_high, 'BinEdges', bin_edges1, 'FaceColor', colors{1}, 'FaceAlpha', 0.6, 'Normalization', 'probability'), hold on
        histogram(Vh_non_class_high, 'BinEdges', bin_edges1, 'FaceColor', colors{2}, 'FaceAlpha', 0.6, 'Normalization', 'probability')
        xlabel('V_{hNa} (mV)')
        
        % Compute common bin edges for the second subplot
        min_edge2 = min([Vm_class_high; Vm_non_class_high]);
        max_edge2 = max([Vm_class_high; Vm_non_class_high]);
        bin_edges2 = linspace(min_edge2, max_edge2, nbins+1); % Ensure same bin edges
        
        subplot(2,1,2)
        histogram(Vm_class_high, 'BinEdges', bin_edges2, 'FaceColor', colors{1}, 'FaceAlpha', 0.6, 'Normalization', 'probability'), hold on
        histogram(Vm_non_class_high, 'BinEdges', bin_edges2, 'FaceColor', colors{2}, 'FaceAlpha', 0.6, 'Normalization', 'probability')
        xlabel('V_{mNa} (mV)')
        exportgraphics(gcf,[savedir,filesep,current_class,'_high.pdf'], 'ContentType', 'vector')

    end
end