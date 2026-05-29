function [values, currents] = Riz2014_rhs(t, states, parameters, varargin)
  % Compute the right hand side of the Riz2014 ODE

  % Assign states
  if length(states)~=14
    error('Expected the states array to be of size 14.');
  end
  h_Na=states(1); h_CaL=states(2); h_CaT=states(3); m_BK=states(4);...
    m_Kv=states(5); m_HERG=states(6); h_HERG=states(7); G6PF6P=states(8);...
    FBP=states(9); DHAPG3P=states(10); a=states(11); Ca_m=states(12);...
    Ca_c=states(13); V=states(14);

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

  % Assign parameters
  if length(parameters)~=80
    error('Expected the parameters array to be of size 80.');
  end
  V_hNa=parameters(1); V_mNa=parameters(2); g_Na=parameters(3);...
    n_hNa=parameters(4); n_mNa=parameters(5); tau_hNa=parameters(6);...
    V_mCaL=parameters(7); g_CaL=parameters(8); n_mCaL=parameters(9);...
    tau_hCaL=parameters(10); V_mCaPQ=parameters(11); g_CaPQ=parameters(12);...
    n_mCaPQ=parameters(13); V_hCaT=parameters(14); V_mCaT=parameters(15);...
    g_CaT=parameters(16); n_hCaT=parameters(17); n_mCaT=parameters(18);...
    tau_hCaT=parameters(19); K_SK=parameters(20); g_SK=parameters(21);...
    n=parameters(22); B_BK=parameters(23); V_mBK=parameters(24);...
    g_BK=parameters(25); n_mBK=parameters(26); tau_mBK=parameters(27);...
    V_mKv=parameters(28); g_Kv=parameters(29); n_mKv=parameters(30);...
    tau_mKv0=parameters(31); V_hHERG=parameters(32); V_mHERG=parameters(33);...
    g_HERG=parameters(34); n_hHERG=parameters(35); n_mHERG=parameters(36);...
    tau_hHERG=parameters(37); tau_mHERG=parameters(38);...
    g_KATP0=parameters(39); g_KATP_hat=parameters(40);...
    glycolysis=parameters(41); g_GABAR=parameters(42); V_leak=parameters(43);...
    g_leak=parameters(44); J_SERCA_max=parameters(45);...
    K_SERCA=parameters(46); J_PMCA_max=parameters(47); K_PMCA=parameters(48);...
    J_leak=parameters(49); J_NCX0=parameters(50); G=parameters(51);...
    K_FBA=parameters(52); K_GAPDH=parameters(53); K_GK=parameters(54);...
    K_GPI=parameters(55); K_PFK=parameters(56); K_TPI=parameters(57);...
    P_FBA=parameters(58); Q_FBA=parameters(59); V_FBA_max=parameters(60);...
    V_GAPDH_max=parameters(61); V_GK_max=parameters(62);...
    V_PFK_max=parameters(63); X_PFK=parameters(64); alpha_G=parameters(65);...
    h_GK=parameters(66); h_PFK=parameters(67); h_X=parameters(68);...
    h_act=parameters(69); k_A=parameters(70); B=parameters(71);...
    Cm=parameters(72); Vol_c=parameters(73); Vol_m=parameters(74);...
    alpha=parameters(75); f=parameters(76); V_Ca=parameters(77);...
    V_Cl=parameters(78); V_K=parameters(79); V_Na=parameters(80);

  % Init return args
  values = zeros(14, 1);

  % Expressions for the I_Na component
  m_Na_inf = 1.0/(1 + exp((-V_mNa + V)/n_mNa));
  h_Na_inf = 1.0/(1 + exp((-V_hNa + V)/n_hNa));
  values(1) = (-h_Na + h_Na_inf)/tau_hNa;
  I_Na = g_Na*(-V_Na + V)*h_Na*m_Na_inf;

  % Expressions for the I_CaL component
  m_CaL_inf = 1.0/(1 + exp((-V_mCaL + V)/n_mCaL));
  h_CaL_inf = ((1 + 0.0175438596491*(-V_Ca + V)*m_CaL_inf > 0)*(1) + ~(1 +...
    0.0175438596491*(-V_Ca + V)*m_CaL_inf > 0)*(0))*((1 +...
    0.0175438596491*(-V_Ca + V)*m_CaL_inf < 1)*(1 + 0.0175438596491*(-V_Ca +...
    V)*m_CaL_inf) + ~(1 + 0.0175438596491*(-V_Ca + V)*m_CaL_inf < 1)*(1));
  values(2) = (-h_CaL + h_CaL_inf)/tau_hCaL;
  I_CaL = g_CaL*(-V_Ca + V)*h_CaL*m_CaL_inf;

  % Expressions for the I_CaPQ component
  m_CaPQ_inf = 1.0/(1 + exp((-V_mCaPQ + V)/n_mCaPQ));
  I_CaPQ = g_CaPQ*(-V_Ca + V)*m_CaPQ_inf;

  % Expressions for the I_CaT component
  m_CaT_inf = 1.0/(1 + exp((-V_mCaT + V)/n_mCaT));
  h_CaT_inf = 1.0/(1 + exp((-V_hCaT + V)/n_hCaT));
  values(3) = (-h_CaT + h_CaT_inf)/tau_hCaT;
  I_CaT = g_CaT*(-V_Ca + V)*h_CaT*m_CaT_inf;
  I_Ca = I_CaL + I_CaPQ + I_CaT;

  % Expressions for the I_SK component
  I_SK = g_SK*Ca_m^n*(-V_K + V)/(K_SK^n + Ca_m^n);

  % Expressions for the I_BK component
  m_BK_inf = 1.0/(1 + exp((-V_mBK + V)/n_mBK));
  values(4) = (-m_BK + m_BK_inf)/tau_mBK;
  I_BK = g_BK*(B_BK - I_Ca)*(-V_K + V)*m_BK;

  % Expressions for the I_Kv component
  tau_mKv = ((V < -26.6)*(30 + tau_mKv0) + ~(V < -26.6)*(tau_mKv0 +...
    10*exp(-10/3 - V/6)));
  m_Kv_inf = 1.0/(1 + exp((-V_mKv + V)/n_mKv));
  values(5) = (-m_Kv + m_Kv_inf)/tau_mKv;
  I_Kv = g_Kv*(-V_K + V)*m_Kv;

  % Expressions for the I_HERG component
  m_HERG_inf = 1.0/(1 + exp((-V_mHERG + V)/n_mHERG));
  h_HERG_inf = 1.0/(1 + exp((-V_hHERG + V)/n_hHERG));
  values(6) = (-m_HERG + m_HERG_inf)/tau_mHERG;
  values(7) = (-h_HERG + h_HERG_inf)/tau_hHERG;
  I_HERG = g_HERG*(-V_K + V)*h_HERG*m_HERG;

  % Expressions for the I_KATP component
  g_KATP = ((glycolysis)*(g_KATP_hat/(1 + a)) + ~(glycolysis)*(g_KATP0));
  I_KATP = (-V_K + V)*g_KATP;

  % Expressions for the I_leak component
  I_leak = g_leak*(-V_leak + V);

  % Expressions for the I_GABAR component
  I_GABAR = g_GABAR*(-V_Cl + V);

  % Expressions for the J_SERCA component
  J_SERCA = J_SERCA_max*Ca_c^2/(K_SERCA^2 + Ca_c^2);

  % Expressions for the J_PMCA component
  J_PMCA = J_PMCA_max*Ca_m/(K_PMCA + Ca_m);

  % Expressions for the J_NCX component
  J_NCX = J_NCX0*Ca_m;

  % Expressions for the Metabolism component
  F6P = K_GPI*G6PF6P/(1 + K_GPI);
  G3P = K_TPI*DHAPG3P/(1 + K_TPI);
  DHAP = -G3P + DHAPG3P;
  h_FBP = h_PFK - (h_PFK - h_act)*FBP/(K_FBA + FBP);
  V_GK = V_GK_max*G^h_GK/(G^h_GK + K_GK^h_GK);
  V_PFK = V_PFK_max*(F6P/K_PFK)^h_FBP/((F6P/K_PFK)^h_FBP + (1 +...
    (FBP/X_PFK)^h_X)/(1 + alpha_G^h_FBP*(FBP/X_PFK)^h_X));
  V_FBA = V_FBA_max*(FBP/K_FBA - DHAP*G3P/(K_FBA*P_FBA*Q_FBA))/(1 + FBP/K_FBA...
    + DHAP/Q_FBA + DHAP*G3P/(P_FBA*Q_FBA));
  V_GAPDH = V_GAPDH_max*G3P/(K_GAPDH + G3P);
  values(8) = -V_PFK + V_GK;
  values(9) = -V_FBA + V_PFK;
  values(10) = -V_GAPDH + 2*V_FBA;
  values(11) = -k_A*a + V_GAPDH;

  % Expressions for the Calcium concentrations component
  values(12) = -Vol_c*f*(B*(-Ca_c + Ca_m) + J_NCX + J_PMCA)/Vol_m +...
    Cm*alpha*f*(-I_CaL - I_CaPQ - I_CaT)/Vol_m;
  values(13) = f*(J_leak - J_SERCA + B*(-Ca_c + Ca_m));

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
    I_Na + I_SK + I_leak - I_stim - I_vclamp;
  values(14) = -I_tot;

  currents = [I_Na, I_BK, I_CaL, I_CaPQ, I_CaT, I_GABAR, I_HERG, I_KATP, I_Kv, I_SK, I_leak, I_stim, I_vclamp];
end