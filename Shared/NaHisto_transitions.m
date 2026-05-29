function NaHisto_transitions(classes,class_names,VhNa,VmNa,nbins,savedir)
    num_classes = length(class_names);
    for c = 1:num_classes
        current_class = class_names{c};
        if strcmp(current_class,'Silent to Bursting')
            inds = find(strcmp(classes, current_class));
            non_inds = find(~strcmp(classes, current_class));
    
            Vh_class = VhNa(inds);
            Vm_class = VmNa(inds);
    
            Vh_non_class = VhNa(non_inds);
            Vm_non_class = VmNa(non_inds);
            
                % Colormap
            colors = {[125, 62, 119] / 255 ...  % violet
                [85, 138, 163] / 255 ...    % blue
                [246, 88, 95] / 255 ...  % red
                }; 
    
            figure
            sgtitle([current_class,' +'],'FontSize',24)
    
            % Compute common bin edges for the first subplot
            min_edge1 = min([Vh_class; Vh_non_class]);
            max_edge1 = max([Vh_class; Vh_non_class]);
            bin_edges1 = linspace(min_edge1, max_edge1, nbins+1); % Ensure same bin edges
    
            subplot(2,1,1)
            histogram(Vh_class, 'BinEdges', bin_edges1, 'FaceColor', colors{1}, 'FaceAlpha', 0.6, 'Normalization', 'probability'), hold on
            histogram(Vh_non_class, 'BinEdges', bin_edges1, 'FaceColor', colors{2}, 'FaceAlpha', 0.6, 'Normalization', 'probability')
            xlabel('V_{hNa} (mV)')
            
            % Compute common bin edges for the second subplot
            min_edge2 = min([Vm_class; Vm_non_class]);
            max_edge2 = max([Vm_class; Vm_non_class]);
            bin_edges2 = linspace(min_edge2, max_edge2, nbins+1); % Ensure same bin edges
            
            subplot(2,1,2)
            histogram(Vm_class, 'BinEdges', bin_edges2, 'FaceColor', colors{1}, 'FaceAlpha', 0.6, 'Normalization', 'probability'), hold on
            histogram(Vm_non_class, 'BinEdges', bin_edges2, 'FaceColor', colors{2}, 'FaceAlpha', 0.6, 'Normalization', 'probability')
            xlabel('V_{mNa} (mV)')
            exportgraphics(gcf,[savedir,filesep,current_class,'_transition.pdf'], 'ContentType', 'vector')

        end
    end
end