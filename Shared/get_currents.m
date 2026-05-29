function [outs] = get_currents(time, values, parameters, model, varargin)
% Extraction of time course of ion currents

    if ~isempty(varargin)
        v_clamp_amps = varargin{1};
        v_clamp_times = varargin{2};
        R_clamp = varargin{3};
    end
    if func2str(model)=="Riz2014_rhs_INa_low"
        I_Na_high = zeros(size(time));
        I_Na_low = zeros(size(time));
    else
        I_Na = zeros(size(time));
    end
    I_BK = zeros(size(time));
    I_CaL = zeros(size(time));
    I_CaPQ = zeros(size(time));
    I_CaL = zeros(size(time));
    I_CaT = zeros(size(time));
    I_GABAR = zeros(size(time));
    I_HERG = zeros(size(time));
    I_KATP = zeros(size(time));
    I_Kv = zeros(size(time));
    I_SK = zeros(size(time));
    I_leak = zeros(size(time));
    I_stim = zeros(size(time));
    I_vclamp = zeros(size(time));

    for i= 1:size(values,1)
        if ~isempty(varargin)
            [~, data] = model(time(i), values(i,:), parameters, v_clamp_amps,v_clamp_times, R_clamp);
        else
            [~, data] = model(time(i), values(i,:), parameters);
        end
        if func2str(model)=="Riz2014_rhs_INa_low"
            I_Na_high(i) = data(1);
            I_Na_low(i) = data(2);
            I_BK(i) = data(3);
            I_CaL(i) = data(4);
            I_CaPQ(i) = data(5);
            I_CaT(i) = data(6);
            I_GABAR(i) = data(7);
            I_HERG(i)= data(8);
            I_KATP(i) = data(9);
            I_Kv(i) = data(10);
            I_SK(i) = data(11);
            I_leak(i) = data(12);
            I_stim(i) = data(13);
            I_vclamp(i) = data(14); 
        else
            I_Na(i) = data(1);
            I_BK(i) = data(2);
            I_CaL(i) = data(3);
            I_CaPQ(i) = data(4);
            I_CaT(i) = data(5);
            I_GABAR(i) = data(6);
            I_HERG(i)= data(7);
            I_KATP(i) = data(8);
            I_Kv(i) = data(9);
            I_SK(i) = data(10);
            I_leak(i) = data(11);
            I_stim(i) = data(12);
            I_vclamp(i) = data(13);
        end
    end
    if func2str(model)=="Riz2014_rhs_INa_low"
        outs = {I_Na_high, I_Na_low, I_BK, I_CaL, I_CaPQ, I_CaT, I_GABAR, I_HERG, I_KATP, I_Kv, I_SK, I_leak, I_stim, I_vclamp};
    else
        outs = {I_Na, I_BK, I_CaL, I_CaPQ, I_CaT, I_GABAR, I_HERG, I_KATP, I_Kv, I_SK, I_leak, I_stim, I_vclamp};
    end
end