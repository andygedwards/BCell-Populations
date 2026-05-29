function KATPHisto(classes_low,classes_high,class_names,KATP,n_bins,savedir)
    num_classes = length(class_names);
    for c = 1:num_classes
        current_class = class_names{c};
        low_inds = find(strcmp(classes_low, current_class));
        high_inds = find(strcmp(classes_high, current_class));
        low_non_inds = find(~strcmp(classes_low, current_class));
        high_non_inds = find(~strcmp(classes_high, current_class));
        
        KATP_low = KATP(low_inds);
        KATP_high = KATP(high_inds);

        KATP_non_low = KATP(low_non_inds);
        KATP_non_high = KATP(high_non_inds);
        
            % Colormap
        colors = {[125, 62, 119] / 255 ...  % violet
            [85, 138, 163] / 255 ...    % blue
            [246, 88, 95] / 255 ...  % red
            }; 



        nbins = n_bins(1);

        figure
        sgtitle([current_class,' (Low glucose)'],'FontSize',24)
        
        % Compute common bin edges for the first subplot
        min_edge1 = min([KATP_low; KATP_non_low]);
        max_edge1 = max([KATP_low; KATP_non_low]);
        bin_edges1 = linspace(min_edge1, max_edge1, nbins+1); % Ensure same bin edges
        
        subplot(2,1,1)
        histogram(KATP_low, 'BinEdges', bin_edges1, 'FaceColor', colors{1}, 'FaceAlpha', 0.6, 'Normalization', 'probability'), hold on
        histogram(KATP_non_low, 'BinEdges', bin_edges1, 'FaceColor', colors{2}, 'FaceAlpha', 0.6, 'Normalization', 'probability')
        xlabel('g_{KATP} (nS/pF)')
        
        % Compute common bin edges for the first subplot
        min_edge1 = min([KATP_high; KATP_non_high]);
        max_edge1 = max([KATP_high; KATP_non_high]);
        bin_edges1 = linspace(min_edge1, max_edge1, nbins+1); % Ensure same bin edges
        
        subplot(2,1,2)
        histogram(KATP_high, 'BinEdges', bin_edges1, 'FaceColor', colors{1}, 'FaceAlpha', 0.6, 'Normalization', 'probability'), hold on
        histogram(KATP_non_high, 'BinEdges', bin_edges1, 'FaceColor', colors{2}, 'FaceAlpha', 0.6, 'Normalization', 'probability')
        xlabel('g_{KATP} (nS/pF)')
        
        exportgraphics(gcf,[savedir,filesep,'KATP_',current_class,'_high.pdf'], 'ContentType', 'vector')

    end
end