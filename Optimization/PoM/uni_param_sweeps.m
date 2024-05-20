function[Quiescent, Stable, Bursting, FailedRepol] = uni_param_sweeps(scaled_params, param_names,mod_param_scaling,mod_param_names,nTrials,states,state_names)

    Bursting = zeros(nTrials);
    FailedRepol = zeros(nTrials);
    Quiescent = zeros(nTrials);
    Stable = zeros(nTrials);
    
    options = [];
    for i = 1:length(mod_param_scaling)
        [T, Y] = ode15s(@Riz2014_rhs,[0,420000], states, options, scaled_params(i,:));
        %Figure 6 Riz
        a_idx = find(strcmp(state_names, 'a'));
        G6PF6P_idx = find(strcmp(state_names, 'G6PF6P'));
        FBP_idx = find(strcmp(state_names, 'FBP'));
        DHAPG3P_idx = find(strcmp(state_names, 'DHAPG3P'));
        g_KATP_hat = scaled_params(i,find(strcmp(param_names, 'g_KATP_hat')));
        G6PF6P = Y(:,G6PF6P_idx);
        FBP = Y(:,FBP_idx);
        DHAPG3P = Y(:,DHAPG3P_idx);
        a = Y(:,a_idx);

        figure
        subplot(5,1,1), plot(T/(1000*60), G6PF6P)
        title(strcat(mod_param_names,'; scaling =',num2str(mod_param_scaling(i)),'; metabolic coupling'))
        subplot(5,1,2), plot(T/(1000*60), FBP)
        subplot(5,1,3), plot(T/(1000*60), DHAPG3P)
        subplot(5,1,4), plot(T/(1000*60), a)
        subplot(5,1,5), plot(T/(1000*60), g_KATP_hat./(1+a))
        subplot(5,1,1), ylabel('G6PF6P')
        subplot(5,1,2), ylabel('FBP')
        subplot(5,1,3), ylabel('DHAPG3P')
        subplot(5,1,4), ylabel('a')
        subplot(5,1,5), ylabel('g_KATP')
        subplot(5,1,5), xlabel('time (mins)')

        V_idx = find(strcmp(state_names, 'V'));
        Cac_idx = find(strcmp(state_names, 'Ca_c'));
        Cam_idx = find(strcmp(state_names, 'Ca_m'));
        V = Y(:, V_idx);
        Ca_c = Y(:, Cac_idx);
        Ca_m = Y(:, Cam_idx);

        figure
        subplot(3,1,1), plot(T/(1000*60), V)
        title(strcat(mod_param_names,'; scaling =',num2str(mod_param_scaling(i)),'; EP outcomes'))
        subplot(3,1,2), plot(T/(1000*60), Ca_c)
        subplot(3,1,3), plot(T/(1000*60), Ca_m)
        subplot(3,1,1), ylabel('V_m (mV)')
        subplot(3,1,2), ylabel('Cytosolic Ca^{2+} (\muM)')
        subplot(3,1,3), ylabel('Membrane Ca^{2+} (\muM)')
        subplot(3,1,3), xlabel('time (mins)')
        list = {'Bursting','Failed repolarisation','Quiescence','Stable'};
       [answers,tf] = listdlg('PromptString',{'Are there any noteworthy instabilities?'},'ListString',list);
        % Handle response
        if ~isempty(find(answers==1,1))
            Bursting(i) = 1;
        end
        if ~isempty(find(answers==2,1))
            FailedRepol(i) = 1;
        end
        if ~isempty(find(answers==3,1))
            Quiescent(i) = 1;
        end
        if ~isempty(find(answers==4,1))
            Stable(i) = 1;
        end
        close all
    end
    cd('/Users/Andy/Dropbox/Simula/Research/Current Projects/PancreaticIslets/Models/AGE/Riz/Long parameter sweeps') 
    save(strcat(char(mod_param_names),'_sweep'),'Bursting','FailedRepol','Quiescent','Stable')
    close all
end

