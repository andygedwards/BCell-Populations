function [p_matrix,p_mat_names,p_mat_inds,high_inds,low_inds] = convert_param(param,param_names,param_inds)
    p_matrix = zeros(length(param),length(param{1}));
    high_count = 0;
    low_count = 0;
    for p_num = 1:length(param)
        if param{p_num}(param_inds(find(strcmp(param_names,'V_hNa'))))==0 && param{p_num}(param_inds(find(strcmp(param_names,'g_Na'))))==0 && param{p_num}(param_inds(find(strcmp(param_names,'V_mNa'))))==0 && param{p_num}(param_inds(find(strcmp(param_names,'n_hNa'))))==0 
            param{p_num}(param_inds(find(strcmp(param_names,'V_hNa'))))=param{p_num}(param_inds(find(strcmp(param_names,'V_hNa_low'))));
            param{p_num}(param_inds(find(strcmp(param_names,'V_mNa'))))=param{p_num}(param_inds(find(strcmp(param_names,'V_mNa_low'))));
            param{p_num}(param_inds(find(strcmp(param_names,'n_hNa'))))=param{p_num}(param_inds(find(strcmp(param_names,'n_hNa_low'))));
            param{p_num}(param_inds(find(strcmp(param_names,'g_Na'))))=param{p_num}(param_inds(find(strcmp(param_names,'g_Na_low'))));
            param{p_num}(param_inds(find(strcmp(param_names,'V_hNa_low')))) = 0;
            param{p_num}(param_inds(find(strcmp(param_names,'V_mNa_low')))) = 0;
            param{p_num}(param_inds(find(strcmp(param_names,'n_hNa_low')))) = 0;
            param{p_num}(param_inds(find(strcmp(param_names,'g_Na_low')))) = 0;
            low_count = low_count+1;
            low_inds(low_count) = p_num;
        else
            high_count = high_count+1;
            high_inds(high_count) = p_num;
        end
        p_matrix(p_num,:) = param{p_num};
    end
    param_inds(find(strcmp(param_names,'V_hNa_low'))) = [];
    param_names(find(strcmp(param_names,'V_hNa_low'))) = [];
    param_inds(find(strcmp(param_names,'V_mNa_low'))) = [];
    param_names(find(strcmp(param_names,'V_mNa_low'))) = [];
    param_inds(find(strcmp(param_names,'n_hNa_low'))) = [];
    param_names(find(strcmp(param_names,'n_hNa_low'))) = [];
    param_inds(find(strcmp(param_names,'g_Na_low'))) = [];
    param_names(find(strcmp(param_names,'g_Na_low'))) = [];
    
    p_mat_inds = param_inds;
    p_mat_names = param_names;
    if low_count == 0
           low_inds = [];
    end
    if high_count == 0
           high_inds = [];
    end
end