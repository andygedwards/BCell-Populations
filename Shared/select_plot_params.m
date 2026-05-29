function [modParam_names,modParam_inds,plot_names,plot_units] = select_plot_params(param_selection, param_names, param_inds)
      plot_names = {""};
      plot_units = {""};
      index_names = ["V_GK_max", "K_GK", "V_PFK_max", "K_PFK", "h_PFK", "V_GAPDH_max", "g_Kv", "g_BK", "g_CaL", "g_CaPQ", "g_CaT", "g_KATP_hat", "g_HERG", "V_hNa", "n_hNa", "g_Na", "V_mNa", "V_hNa_low", "n_hNa_low", "g_Na_low", "V_mNa_low"];
      index_plot_names = {'V_{GKmax}', 'K_{GK}', 'V_{PFKmax}', 'K_{PFK}', 'h_{PFK}', 'V_{GAPDHmax}', 'g_{Kv}', 'g_{BK}', 'g_{CaL}', 'g_{CaPQ}', 'g_{CaT}', 'g_{KATP}', 'g_{hERG}', 'V_{hNa}', 'n_{hNa}', 'g_{Na}', 'V_{mNa}', 'V_{hNa}', 'n_{hNa}', 'g_{Na}', 'V_{mNa}'}; 
      index_plot_units = {'M/s', 'mM', 'M/s', 'mM', 'dimensionless', 'M/s', 'nS/pF', 'nS/pF','nS/pF','nS/pF','nS/pF','nS/pF','nS/pF','mV','dimensionless','nS/pF','mV','mV','dimensionless','nS/pF','mV'};
      for ploop = 1:length(param_selection)
          if find(strcmp(param_names,param_selection(ploop)))
              modParam_names(ploop) = param_selection(ploop);
              modParam_inds(ploop) = param_inds(find(strcmp(param_names,param_selection(ploop))));
              target = index_plot_names{find(strcmp(index_names, param_selection(ploop)))};
              if ~isempty(find(cellfun(@(x) contains(x, target), plot_names)))
                disp('skipped')
                continue
              else
                disp(index_plot_names{find(strcmp(index_names,param_selection(ploop)))})
                plot_names{ploop} = index_plot_names{find(strcmp(index_names,param_selection(ploop)))};
                plot_units{ploop} = index_plot_units{find(strcmp(index_names,param_selection(ploop)))}; 
              end
          end
      end 
end