function [outs] = analyze_v_clamp(Y, T, parameters, state_names, CL, PP, EL, model, varargin)

% Analysis and plotting

V_idx = find(strcmp(state_names, 'V'));
Cac_idx = find(strcmp(state_names, 'Ca_c'));
Cam_idx = find(strcmp(state_names, 'Ca_m'));
V = Y(:, V_idx);
Ca_c = Y(:, Cac_idx);
Ca_m = Y(:, Cam_idx);

if ~isempty(varargin)
    v_clamp_amps = varargin{1};
    v_clamp_times = varargin{2};
    R_clamp = varargin{3};
    current = varargin{4};
    plotting = varargin{5};
    currents = get_currents(T, Y, parameters, model, v_clamp_amps, v_clamp_times, R_clamp);
    if func2str(model)=="Riz2014_rhs_INa_low"
        I_Na_high = currents{1}; 
        I_Na_low = currents{2};
        I_BK = currents{3};
        I_CaL = currents{4};
        I_CaPQ = currents{5};
        I_CaT = currents{6};
        I_GABAR = currents{7};
        I_HERG = currents{8};
        I_KATP = currents{9};
        I_Kv = currents{10};
        I_SK = currents{11};
        I_leak = currents{12};
        I_stim = currents{13};
        I_vclamp = currents{14};
    else
        I_Na = currents{1}; 
        I_BK = currents{2};
        I_CaL = currents{3};
        I_CaPQ = currents{4};
        I_CaT = currents{5};
        I_GABAR = currents{6};
        I_HERG = currents{7};
        I_KATP = currents{8};
        I_Kv = currents{9};
        I_SK = currents{10};
        I_leak = currents{11};
        I_stim = currents{12};
        I_vclamp = currents{13};
    end
else
    currents = get_currents(T, Y, parameters, model);
    if func2str(model)=="Riz2014_rhs_INa_low"
        I_Na_high = currents{1}; 
        I_Na_low = currents{2};
        I_BK = currents{3};
        I_CaL = currents{4};
        I_CaPQ = currents{5};
        I_CaT = currents{6};
        I_GABAR = currents{7};
        I_HERG = currents{8};
        I_KATP = currents{9};
        I_Kv = currents{10};
        I_SK = currents{11};
        I_leak = currents{12};
        I_stim = currents{13};
        I_vclamp = currents{14};
    else
        I_Na = currents{1}; 
        I_BK = currents{2};
        I_CaL = currents{3};
        I_CaPQ = currents{4};
        I_CaT = currents{5};
        I_GABAR = currents{6};
        I_HERG = currents{7};
        I_KATP = currents{8};
        I_Kv = currents{9};
        I_SK = currents{10};
        I_leak = currents{11};
        I_stim = currents{12};
        I_vclamp = currents{13};
    end
end

if func2str(model)=="Riz2014_rhs_INa_low"
    I = [I_Na_high, I_Na_low, I_BK, I_CaL, I_CaPQ, I_CaT, I_GABAR, I_HERG, I_KATP, I_Kv, I_SK, I_leak, I_stim, I_vclamp];
    I_in = [I_Na_high, I_Na_low, I_CaL, I_CaPQ, I_CaT];
else
    I = [I_Na, I_BK, I_CaL, I_CaPQ, I_CaT, I_GABAR, I_HERG, I_KATP, I_Kv, I_SK, I_leak, I_stim, I_vclamp];
    I_in = [I_Na, I_CaL, I_CaPQ, I_CaT];
end
I_tot = sum(I(:,1:end-2),2);
I_in_tot = sum(I_in,2);

% if plotting 
%     figure
%     subplot(2,2,1), plot(T, V)
%     subplot(2,2,2), plot(T, Ca_c)
%     subplot(2,2,3), plot(T, Ca_m)
%     subplot(2,2,4), plot(T, I_tot)
%     subplot(2,2,1), ylabel('V_m (mV)')
%     subplot(2,2,2), ylabel('Cytosolic Ca^{2+} (\muM)')
%     subplot(2,2,3), ylabel('Membrane Ca^{2+} (\muM)')
%     subplot(2,2,4), ylabel('Total current (pA/pF)')
%     subplot(2,2,3), xlabel('time (ms)')
%     subplot(2,2,4), xlabel('time (ms)')
% end

% Decompose sweeps
end_ind = find(T<=T(end)-EL,1,"last");
T = T(1:end_ind);
Y = Y(1:end_ind,:);
V = V(1:end_ind);
Ca_c = Ca_c(1:end_ind);
Ca_m = Ca_m(1:end_ind);
I_tot = I_tot(1:end_ind);
I_in_tot = I_in_tot(1:end_ind);
I_CaL = I_CaL(1:end_ind);
I_CaPQ = I_CaPQ(1:end_ind);
I_CaT = I_CaT(1:end_ind);
if func2str(model)=="Riz2014_rhs_INa_low"
    I_Na_high = I_Na_high(1:end_ind);
    I_Na_low = I_Na_low(1:end_ind);
else
    I_Na = I_Na(1:end_ind);
end

start_ind = find(T>=PP,1,"first");
T = T(start_ind:end);
Y = Y(start_ind:length(Y(:,1)),:);
V = V(start_ind:end);
Ca_c = Ca_c(start_ind:end);
Ca_m = Ca_m(start_ind:end);
I_tot = I_tot(start_ind:end);
I_in_tot = I_in_tot(start_ind:end);
I_CaL = I_CaL(start_ind:end);
I_CaPQ = I_CaPQ(start_ind:end);
I_CaT = I_CaT(start_ind:end);
if func2str(model)=="Riz2014_rhs_INa_low"
    I_Na_high = I_Na_high(start_ind:end);
    I_Na_low = I_Na_low(start_ind:end);
else
    I_Na = I_Na(start_ind:end);
end

T = T-PP;

if ~strcmp(current,'ramp') && ~strcmp(current,'K_inact_2')
    num_pulses = round(max(T),0)/CL;
%     figure
    for i = 1:num_pulses
        if i == 1
            sweep_start = 1;
        else
            sweep_start = find(T>=(i-1)*CL,1,"first");
        end
        sweep_end = find(T>=i*CL,1,"first");
        Sweep_Y = Y(sweep_start:sweep_end,:);
        Sweep_T = T(sweep_start:sweep_end)-T(sweep_start);
        Sweep_I = I_tot(sweep_start:sweep_end);
        Sweep_I_in = I_in_tot(sweep_start:sweep_end);
        Sweep_I_CaL = I_CaL(sweep_start:sweep_end);
        Sweep_I_CaPQ = I_CaPQ(sweep_start:sweep_end);
        Sweep_I_CaT = I_CaT(sweep_start:sweep_end);
        if func2str(model)=="Riz2014_rhs_INa_low"
            Sweep_I_Na_high = I_Na_high(sweep_start:sweep_end);
            Sweep_I_Na_low = I_Na_low(sweep_start:sweep_end);
        else
            Sweep_I_Na = I_Na(sweep_start:sweep_end);
        end
        Sweep_V = Sweep_Y (:, V_idx);
        Sweep_Ca_c = Sweep_Y(:, Cac_idx);
        Sweep_Ca_m = Sweep_Y(:, Cam_idx);
        
%         if plotting
%             subplot(2,2,1), plot(Sweep_T, Sweep_V), hold on
%             subplot(2,2,2), plot(Sweep_T, Sweep_Ca_c), hold on
%             subplot(2,2,3), plot(Sweep_T, Sweep_Ca_m), hold on
%             subplot(2,2,4), plot(Sweep_T, Sweep_I), hold on
%             subplot(2,2,1), ylabel('V_m (mV)')
%             subplot(2,2,2), ylabel('Cytosolic Ca^{2+} (\muM)')
%             subplot(2,2,3), ylabel('Membrane Ca^{2+} (\muM)')
%             subplot(2,2,4), ylabel('total current (pA/pF)')
%             subplot(2,2,3), xlabel('time (ms)'), subplot(2,2,4), xlabel('time (ms)')
%         end
    
        Sweeps_Y{i,:,:} = Sweep_Y;
        Sweeps_T{i,:} = Sweep_T;
        Sweeps_V{i,:} = Sweep_V;
        Sweeps_Ca_c{i,:} = Sweep_Ca_c;
        Sweeps_Ca_m{i,:} = Sweep_Ca_m;
        Sweeps_I{i,:} = Sweep_I;
        Sweeps_I_in{i,:} = Sweep_I_in;
        Sweeps_I_CaL{i,:} = Sweep_I_CaL;
        Sweeps_I_CaPQ{i,:} = Sweep_I_CaPQ;
        Sweeps_I_CaT{i,:} = Sweep_I_CaT;
        if func2str(model)=="Riz2014_rhs_INa_low"
            Sweeps_I_Na_high{i,:} = Sweep_I_Na_high;
            Sweeps_I_Na_low{i,:} = Sweep_I_Na_low;
        else
            Sweeps_I_Na{i,:} = Sweep_I_Na;
        end
    end
end

if ~isempty(varargin)
    switch current
        case 'exocytosis'
            [early_peak,early_ind] = max(Sweeps_Ca_c{1,:});
            [late_peak, late_ind] = max(cell2mat(Sweeps_Ca_c));

            early_exocytosis = max(Sweeps_Ca_c{1,:})-min(Sweeps_Ca_c{1,:});
            late_exocytosis = max(cell2mat(Sweeps_Ca_c))-early_exocytosis;
            total_exocytosis = early_exocytosis + late_exocytosis;

%             figure
%             plot(T,Ca_c), hold on
% %             plot(T(early_ind),early_peak,'r*',T(late_ind),late_peak,'b*');
%             xlabel('time (ms)')
%             ylabel('cytosolic [Ca^{2+}] (\mu{}M)')
            outs(1) = {early_exocytosis};
            outs(2) = {late_exocytosis};
            outs(3) = {total_exocytosis};
        case 'I_Na'
            try
                V = (-120:10:-20);
%                 figure
                for i = 1:length(V)
                    start_pulse = find(Sweeps_T{i}>502,1,'first');
                    end_plat = find(Sweeps_T{i}<500,1,'last');
                    start_plat = find(Sweeps_T{i}>400,1,'first');
                    [peak_INa(i),peak_INa_ind(i)] = min(Sweeps_I{i}(start_pulse:end));
                    [plat_INa(i)] = mean(Sweeps_I{i}(start_plat:end_plat));
%                     if plotting
%                         subplot(2,1,1), plot(Sweeps_T{i}, Sweeps_V{i}), hold on
%                         subplot(2,1,2), plot(Sweeps_T{i}, Sweeps_I{i}), hold on 
% %                         plot(Sweeps_T{i}(peak_INa_ind(i)+start_pulse), peak_INa(i), 'r*')
% %                         plot(Sweeps_T{i}(end_plat), plat_INa(i), 'k*')
% %                         pause
%                     end
                end
%                 subplot(2,1,1), ylabel('V_{m} (mV)')
%                 subplot(2,1,2), ylabel('I_{Na} (A/F)')
%                 subplot(2,1,2), xlabel('time (ms)')
                peak_INa = peak_INa-plat_INa;
                peak_INa = peak_INa-(max(peak_INa));
%                 max_peak_INa = max(abs(peak_INa))-min(abs(peak_INa));
                peak_INa_norm = peak_INa./min(peak_INa);
                boltzmann = '1/(1+exp((x-a)/b))';
                start = [-40,9];
                [f1,gof] = fit(V',peak_INa_norm',boltzmann,'Start',start);
                % if plotting
                %     figure
                %     plot(f1,V,peak_INa_norm,'r*')
                %     disp(strcat('v_half',num2str(f1.a)))
                % end
%                 pause
                v_half = f1.a;
                outs(1) = {v_half};
                outs(2) = {gof.rsquare};
            catch
                outs(1) = {NaN};
            end
        case 'peak currents'
            try
                V = (-60:10:30);
                %%%  CHECK THIS!!!! %%%
                base_start = find(Sweeps_T{i}>900,1,'first');
                base_end = find(Sweeps_T{i}>900,1,'first');
                base_I(i) = mean(Sweeps_I{i}(base_start:base_end));
                %%%  CHECK THIS!!!! %%%
%                 figure
                for i = 1:length(V)
                    %Late I_Ca
                    late_ICa(i) = Sweeps_I{i}(end)-base_I(i);%Sweeps_I_CaT{i}(end)+Sweeps_I_CaL{i}(end)+Sweeps_I_CaPQ{i}(end);

                    %Peak I_Ca
                    I_Ca_peak_start = find(Sweeps_T{i}>=1010,1,"first");
                    I_Ca_peak_end = find(Sweeps_T{i}>=1100,1,"first");
                    [peak_ICa(i), peak_ICa_ind(i)] = min(Sweeps_I{i}(I_Ca_peak_start:I_Ca_peak_end));%Sweeps_I_CaT{i}(I_Ca_peak_start:I_Ca_peak_end)+Sweeps_I_CaL{i}(I_Ca_peak_start:I_Ca_peak_end)+Sweeps_I_CaPQ{i}(I_Ca_peak_start:I_Ca_peak_end));
                    peak_ICa_ind(i) = peak_ICa_ind(i)+I_Ca_peak_start;
                    peak_ICa(i) = peak_ICa(i)-base_I(i);
                    
                    %Peak I_Na
                    [peak_INa(i),peak_INa_ind(i)] = min(Sweeps_I{i});
                    peak_INa(i) = peak_INa(i)-base_I(i);

%                     if plotting
%                         subplot(2,2,1), plot(Sweeps_T{i}, Sweeps_I{i}), hold on
%                         subplot(2,2,2), plot(Sweeps_T{i}, Sweeps_I_in{i}), hold on
%                         subplot(2,2,3), plot(Sweeps_T{i}, Sweeps_I_CaL{i},'b'), hold on
%                         subplot(2,2,3), plot(Sweeps_T{i}, Sweeps_I_CaPQ{i},'k')
%                         subplot(2,2,3), plot(Sweeps_T{i}, Sweeps_I_CaT{i},'r')
%                         if func2str(model)=="Riz2014_rhs_INa_low"
%                             subplot(2,2,3), plot(Sweeps_T{i}, Sweeps_I_Na_high{i},'g')
%                             subplot(2,2,3), plot(Sweeps_T{i}, Sweeps_I_Na_low{i},'g')
%                         else
%                             subplot(2,2,3), plot(Sweeps_T{i}, Sweeps_I_Na{i},'g')
%                         end
%                         subplot(2,2,4), plot(Sweeps_T{i}, Sweeps_V{i}), hold on
% %                         subplot(2,2,4), plot(Sweeps_T{i}(peak_ICa_ind(i)), peak_ICa(i), 'r*')
% %                         subplot(2,2,4), plot(Sweeps_T{i}(end), late_ICa(i), 'r*')
% %                         subplot(2,2,4), plot(Sweeps_T{i}(peak_INa_ind(i)), peak_INa(i), 'k*')
%                         subplot(2,2,1), ylabel('Total current (pA/pF)')
%                         subplot(2,2,2), ylabel('Total inward current (pA/pF)')
%                         subplot(2,2,3), ylabel('Individual inward currents (pA/pF)')
%                         subplot(2,2,4), ylabel('V_{m} (mV)')
%                         subplot(2,2,3), legend('I_{CaL}','I_{CaPQ}', 'I_{CaT}', 'I_{Na}')
%                         subplot(2,2,3), xlabel('time (ms)'), subplot(2,2,4), xlabel('time (ms)')
% %                         pause
%                     end
                end

                % Pseudo-activation curve
                peak_INa_subt = peak_INa-(max(peak_INa));
                peak_INa_norm = peak_INa_subt./min(peak_INa);
                boltzmann = '1/(1+exp((a-x)/b))';
                start = [-20,5];
                [f1,gof] = fit(V',peak_INa_norm',boltzmann,'Start',start);
                % if plotting
                %     figure
                %     plot(f1,V,peak_INa_norm,'r*')
                %     disp(strcat('v_half',num2str(f1.a)))
                % end
  %              pause

                % Find overall peaks and assign outputs
                peak_INa = min(peak_INa);
                peak_ICa = min(peak_ICa);
                late_ICa = min(late_ICa);
                outs(1) = {peak_INa};
                outs(2) = {f1.a};
                outs(3) = {gof.rsquare};
                outs(4) = {peak_ICa};
                outs(5) = {late_ICa};
            catch
                disp(strcat(num2str(i), ':', num2str(min(Sweeps_T{i,:})),':', num2str(max(Sweeps_T{i,:})),':',num2str((i-1)*CL),':',num2str((i)*CL), ':', num2str(max(Sweeps_T{1,:})), ':',num2str((1)*CL)));
                outs(1) = {NaN};
                outs(2) = {NaN};
                outs(3) = {NaN};
                outs(4) = {NaN};
                outs(5) = {NaN};
            end
        case 'ramp'
            try
               [I_Kmax, I_Kmax_ind] = max(I_tot); 
               outs(1) = {I_Kmax};
%                figure
%                subplot(2,1,1), plot(T,V);
%                subplot(2,1,2), plot(T, I_tot,T(I_Kmax_ind),I_Kmax,'r*');
            catch
                outs(1) = {NaN};
            end
        case 'K_square_1'
            try
                V = (-60:10:50);
%                 figure
                for i = 1:length(V)
                    start_ind = find(Sweeps_T{i}>=1000,1,'first');
                    [peak_K_square(i),peak_K_square_ind(i)]= max(Sweeps_I{i,:}(start_ind:end));
                    peak_K_square_ind(i) = peak_K_square_ind(i)+start_ind-1;
%                     subplot(2,1,1), plot(Sweeps_T{i},Sweeps_V{i}); hold on
%                     subplot(2,1,2), plot(Sweeps_T{i},Sweeps_I{i},Sweeps_T{i}(peak_K_square_ind(i)),Sweeps_I{i}(peak_K_square_ind(i)),'r*'); hold on
                end
%                 figure
%                 plot(V, peak_K_square)
                outs(1) = {peak_K_square(5)};
                outs(2) = {max(peak_K_square)};
            catch
                outs(1) = {NaN};
                outs(2) = {NaN};
            end

         case 'K_tail_1'
            try
                V = (-70:10:40);
%                 figure
                for i = 1:length(V)
                    start_ind = find(Sweeps_T{i}>=5100.1,1,'first');
                    [peak_K_tail_1(i),peak_K_tail_1_ind(i)]= max(Sweeps_I{i,:}(start_ind:end));
                    peak_K_tail_1_ind(i) = peak_K_tail_1_ind(i)+start_ind-1;
%                     subplot(2,1,1), plot(Sweeps_T{i},Sweeps_V{i}); hold on
%                     subplot(2,1,2), plot(Sweeps_T{i},Sweeps_I{i},Sweeps_T{i}(peak_K_tail_1_ind(i)),Sweeps_I{i}(peak_K_tail_1_ind(i)),'r*'); hold on
                end
                max_K_tail_1 = max(peak_K_tail_1);
                peak_K_tail_1_norm = peak_K_tail_1./max_K_tail_1;
                boltzmann = '1/(1+exp(-(x-a)/b))';
                start = [-5,9];
                f1 = fit(V',peak_K_tail_1_norm',boltzmann,'Start',start);
%                 if plotting
%                     figure
%                     plot(f1,V,peak_K_tail_1_norm,'r*')
%                 end
                v_half_act_K = f1.a;
                n_act_K = f1.b;
                outs(1) = {v_half_act_K};
                outs(2) = {n_act_K};
            catch
                outs(1) = {NaN};
                outs(2) = {NaN};
            end
         case 'K_inact_1'
             try
                V = (-120:10:-30);
%                 figure
                for i = 1:length(V)
                    start_ind = find(Sweeps_T{i}>=10010.1,1,'first');
                    [peak_K_inact(i),peak_K_inact_ind(i)]= max(Sweeps_I{i,:}(start_ind:end));
                    peak_K_inact_ind(i) = peak_K_inact_ind(i)+start_ind-1;
%                     if plotting
%                         subplot(2,1,1), plot(Sweeps_T{i}, Sweeps_I{i}), hold on 
%                         subplot(2,1,1), plot(Sweeps_T{i}(peak_K_inact_ind(i)), peak_K_inact(i), 'r*')
%                         subplot(2,1,2), plot(Sweeps_T{i}, Sweeps_V{i}), hold on
%                     end
                end
                max_IK_inact = max(peak_K_inact);
                peak_K_inact_norm = peak_K_inact./max_IK_inact;
                boltzmann = '1/(1+exp((x-a)/b))';
                start = [-45,7];
                f1 = fit(V',peak_K_inact_norm',boltzmann,'Start',start);
%                 if plotting
%                     figure
%                     plot(f1,V,peak_K_inact_norm,'r*')
%                 end
                v_half_K_inact = f1.a;
                outs(1) = {v_half_K_inact};
             catch
                outs(1) = {NaN};
             end
         case 'K_inact_2'
             try
                start_ind = find(T>=10000.1,1,'first');
                [peak_K_inact2,peak_K_inact2_ind]= max(I_tot(start_ind:end));
                fit_interval_T = T(peak_K_inact2_ind:end);            
                fit_interval_I = I_tot(peak_K_inact2_ind:end)./peak_K_inact2;
                double_expo = 'A*exp(a*x)+B*exp(b*x)';
                start = [0.8,2500,0.2,130];
                f1 = fit(fit_interval_T',fit_interval_I',double_expo,'Start',start);
%                 if plotting
%                     figure
%                     plot(f1,fit_interval_T',fit_interval_I','k')
%                 end
                Amp_1 = f1.A;
                tau_1 = f1.a;
                Amp_2 = f1.B;
                tau_2 = f1.b;
                outs(1) = {Amp_1};
                outs(2) = {tau_1};
                outs(3) = {Amp_2};
                outs(4) = {tau_2};
             catch
                outs(1) = {NaN};
                outs(2) = {NaN};
                outs(3) = {NaN};
                outs(4) = {NaN};
             end
    end
%     pause
%     close all
end



