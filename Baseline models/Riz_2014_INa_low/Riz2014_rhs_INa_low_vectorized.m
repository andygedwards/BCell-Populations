function [values] = Riz2014_rhs_INa_low_vectorized(t, states, parameters, varargin)
  % Compute the right hand side of the Riz2014 ODE

  h_Na=states(1,:); h_Na_low=states(2,:); h_CaL=states(3,:); h_CaT=states(4,:); m_BK=states(5,:);...
    m_Kv=states(6,:); m_HERG=states(7,:); h_HERG=states(8,:); G6PF6P=states(9,:);...
    FBP=states(10,:); DHAPG3P=states(11,:); a=states(12,:); Ca_m=states(13,:);...
    Ca_c=states(14,:); V=states(15,:);

  % voltage clamp parameters
  if ~isempty(varargin)
      if length(varargin) < 2
          error('Voltage clamp parameters must be passed as a 3-element cell array: varargin{1} = command step amplitude (mV); varargin{2} = step duration (ms); varargin{3} = Clamp resistance (GOhm).');
      else
        v_clamp_amps = varargin{1};
        v_clamp_times = varargin{2};
        R_clamp = varargin{3};
      end
  end


  V_hNa=parameters(1,:); V_mNa=parameters(2,:); g_Na=parameters(3,:);...
    n_hNa=parameters(4,:); n_mNa=parameters(5,:); tau_hNa=parameters(6,:);...
    V_hNa_low=parameters(7,:); V_mNa_low=parameters(8,:); g_Na_low=parameters(9,:);...
    n_hNa_low=parameters(10,:); n_mNa_low=parameters(11,:); tau_hNa_low=parameters(12,:);...
    V_mCaL=parameters(13,:); g_CaL=parameters(14,:); n_mCaL=parameters(15,:);...
    tau_hCaL=parameters(16,:); V_mCaPQ=parameters(17,:); g_CaPQ=parameters(18,:);...
    n_mCaPQ=parameters(19,:); V_hCaT=parameters(20,:); V_mCaT=parameters(21,:);...
    g_CaT=parameters(22,:); n_hCaT=parameters(23,:); n_mCaT=parameters(24,:);...
    tau_hCaT=parameters(25,:); K_SK=parameters(26,:); g_SK=parameters(27,:);...
    n=parameters(28,:); B_BK=parameters(29,:); V_mBK=parameters(30,:);...
    g_BK=parameters(31,:); n_mBK=parameters(32,:); tau_mBK=parameters(33,:);...
    V_mKv=parameters(34,:); g_Kv=parameters(35,:); n_mKv=parameters(36,:);...
    tau_mKv0=parameters(37,:); V_hHERG=parameters(38,:); V_mHERG=parameters(39,:);...
    g_HERG=parameters(40,:); n_hHERG=parameters(41,:); n_mHERG=parameters(42,:);...
    tau_hHERG=parameters(43,:); tau_mHERG=parameters(44,:);...
    g_KATP0=parameters(45,:); g_KATP_hat=parameters(46,:);...
    glycolysis=parameters(47,:); g_GABAR=parameters(48,:); V_leak=parameters(49,:);...
    g_leak=parameters(50,:); J_SERCA_max=parameters(51,:);...
    K_SERCA=parameters(52,:); J_PMCA_max=parameters(53,:); K_PMCA=parameters(54,:);...
    J_leak=parameters(55,:); J_NCX0=parameters(56,:); G=parameters(57,:);...
    K_FBA=parameters(58,:); K_GAPDH=parameters(59,:); K_GK=parameters(60,:);...
    K_GPI=parameters(61,:); K_PFK=parameters(62,:); K_TPI=parameters(63,:);...
    P_FBA=parameters(64,:); Q_FBA=parameters(65,:); V_FBA_max=parameters(66,:);...
    V_GAPDH_max=parameters(67,:); V_GK_max=parameters(68,:);...
    V_PFK_max=parameters(69,:); X_PFK=parameters(70,:); alpha_G=parameters(71,:);...
    h_GK=parameters(72,:); h_PFK=parameters(73,:); h_X=parameters(74,:);...
    h_act=parameters(75,:); k_A=parameters(76,:); B=parameters(77,:);...
    Cm=parameters(78,:); Vol_c=parameters(79,:); Vol_m=parameters(80,:);...
    alpha=parameters(81,:); f=parameters(82,:); V_Ca=parameters(83,:);...
    V_Cl=parameters(84,:); V_K=parameters(85,:); V_Na=parameters(86,:);

  % Init return args
  values = zeros(size(states));

  % Expressions for the I_Na component
  m_Na_inf = 1.0./(1 + exp((-V_mNa + V)./n_mNa));
  h_Na_inf = 1.0./(1 + exp((-V_hNa + V)./n_hNa));
  values(1,:) = (-h_Na + h_Na_inf)./tau_hNa;
  I_Na_high = g_Na.*(-V_Na + V).*h_Na.*m_Na_inf;

  % Expressions for the I_Na_low component
  m_Na_inf_low = 1.0./(1 + exp((-V_mNa_low + V)./n_mNa_low));
  h_Na_inf_low = 1.0./(1 + exp((-V_hNa_low + V)./n_hNa_low));
  values(2,:) = (-h_Na_low + h_Na_inf_low)./tau_hNa_low;
  I_Na_low = g_Na_low.*(-V_Na + V).*h_Na_low.*m_Na_inf_low;
  I_Na_tot = I_Na_high+I_Na_low;

  % Expressions for the I_CaL component
  m_CaL_inf = 1.0./(1 + exp((-V_mCaL + V)./n_mCaL));
  h_CaL_inf = ((1 + 0.0175438596491.*(-V_Ca + V).*m_CaL_inf > 0).*(1) + ~(1 +...
    0.0175438596491.*(-V_Ca + V).*m_CaL_inf > 0).*(0)).*((1 +...
    0.0175438596491.*(-V_Ca + V).*m_CaL_inf < 1).*(1 + 0.0175438596491.*(-V_Ca +...
    V).*m_CaL_inf) + ~(1 + 0.0175438596491.*(-V_Ca + V).*m_CaL_inf < 1).*(1));
  values(3,:) = (-h_CaL + h_CaL_inf)./tau_hCaL;
  I_CaL = g_CaL.*(-V_Ca + V).*h_CaL.*m_CaL_inf;

  % Expressions for the I_CaPQ component
  m_CaPQ_inf = 1.0./(1 + exp((-V_mCaPQ + V)./n_mCaPQ));
  I_CaPQ = g_CaPQ.*(-V_Ca + V).*m_CaPQ_inf;

  % Expressions for the I_CaT component
  m_CaT_inf = 1.0./(1 + exp((-V_mCaT + V)./n_mCaT));
  h_CaT_inf = 1.0./(1 + exp((-V_hCaT + V)./n_hCaT));
  values(4,:) = (-h_CaT + h_CaT_inf)./tau_hCaT;
  I_CaT = g_CaT.*(-V_Ca + V).*h_CaT.*m_CaT_inf;
  I_Ca = I_CaL + I_CaPQ + I_CaT;

  % Expressions for the I_SK component
  I_SK = g_SK.*Ca_m.^n.*(-V_K + V)./(K_SK.^n + Ca_m.^n);

  % Expressions for the I_BK component
  m_BK_inf = 1.0./(1 + exp((-V_mBK + V)./n_mBK));
  values(5,:) = (-m_BK + m_BK_inf)./tau_mBK;
  I_BK = g_BK.*(B_BK - I_Ca).*(-V_K + V).*m_BK;

  % Expressions for the I_Kv component
  tau_mKv = ((V < -26.6).*(30 + tau_mKv0) + ~(V < -26.6).*(tau_mKv0 +...
    10.*exp(-10./3 - V./6)));
  m_Kv_inf = 1.0./(1 + exp((-V_mKv + V)./n_mKv));
  values(6,:) = (-m_Kv + m_Kv_inf)./tau_mKv;
  I_Kv = g_Kv.*(-V_K + V).*m_Kv;

  % Expressions for the I_HERG component
  m_HERG_inf = 1.0./(1 + exp((-V_mHERG + V)./n_mHERG));
  h_HERG_inf = 1.0./(1 + exp((-V_hHERG + V)./n_hHERG));
  values(7,:) = (-m_HERG + m_HERG_inf)./tau_mHERG;
  values(8,:) = (-h_HERG + h_HERG_inf)./tau_hHERG;
  I_HERG = g_HERG.*(-V_K + V).*h_HERG.*m_HERG;

  % Expressions for the I_KATP component
  g_KATP = ((glycolysis).*(g_KATP_hat./(1 + a)) + ~(glycolysis).*(g_KATP0));
  I_KATP = (-V_K + V).*g_KATP;

  % Expressions for the I_leak component
  I_leak = g_leak.*(-V_leak + V);

  % Expressions for the I_GABAR component
  I_GABAR = g_GABAR.*(-V_Cl + V);

  % Expressions for the J_SERCA component
  J_SERCA = J_SERCA_max.*Ca_c.^2./(K_SERCA.^2 + Ca_c.^2);

  % Expressions for the J_PMCA component
  J_PMCA = J_PMCA_max.*Ca_m./(K_PMCA + Ca_m);

  % Expressions for the J_NCX component
  J_NCX = J_NCX0.*Ca_m;

  % Expressions for the Metabolism component
  F6P = K_GPI.*G6PF6P./(1 + K_GPI);
  G3P = K_TPI.*DHAPG3P./(1 + K_TPI);
  DHAP = -G3P + DHAPG3P;
  h_FBP = h_PFK - (h_PFK - h_act).*FBP./(K_FBA + FBP);
  V_GK = V_GK_max.*G.^h_GK./(G.^h_GK + K_GK.^h_GK);
  V_PFK = V_PFK_max.*(F6P./K_PFK).^h_FBP./((F6P./K_PFK).^h_FBP + (1 +...
    (FBP./X_PFK).^h_X)./(1 + alpha_G.^h_FBP.*(FBP./X_PFK).^h_X));
  V_FBA = V_FBA_max.*(FBP./K_FBA - DHAP.*G3P./(K_FBA.*P_FBA.*Q_FBA))./(1 + FBP./K_FBA...
    + DHAP./Q_FBA + DHAP.*G3P./(P_FBA.*Q_FBA));
  V_GAPDH = V_GAPDH_max.*G3P./(K_GAPDH + G3P);
  values(9,:) = -V_PFK + V_GK;
  values(10,:) = -V_FBA + V_PFK;
  values(11,:) = -V_GAPDH + 2.*V_FBA;
  values(12,:) = -k_A.*a + V_GAPDH;

  % Expressions for the Calcium concentrations component
  values(13,:) = -Vol_c.*f.*(B.*(-Ca_c + Ca_m) + J_NCX + J_PMCA)./Vol_m +...
    Cm.*alpha.*f.*(-I_CaL - I_CaPQ - I_CaT)./Vol_m;
  values(14,:) = f.*(J_leak - J_SERCA + B.*(-Ca_c + Ca_m));

  % Apply voltage clamp protocol
  if ~isempty(varargin)
      if t <= v_clamp_times(1)
          v_clamp_amp = v_clamp_amps(1);
          I_vclamp = (v_clamp_amp - V)/R_clamp;
          I_stim = 0;
      else 
          for i = 2:length(v_clamp_times)
               if and(t >= v_clamp_times(i-1),t <= v_clamp_times(i))
                    v_clamp_amp = v_clamp_amps(i);
                    I_vclamp = (v_clamp_amp - V)/R_clamp;
                    I_stim = 0;
                    break
               end
          end
      end
  else
      I_vclamp = 0;
      I_stim = 0;
  end

  % Expressions for the Membrane potential component
  I_tot = I_BK + I_CaL + I_CaPQ + I_CaT + I_GABAR + I_HERG + I_KATP + I_Kv +...
    I_Na_tot + I_SK + I_leak - I_stim - I_vclamp;
  values(15,:) = -I_tot;

  %currents = [I_Na_high, I_Na_low, I_BK, I_CaL, I_CaPQ, I_CaT, I_GABAR, I_HERG, I_KATP, I_Kv, I_SK, I_leak, I_stim, I_vclamp];
end