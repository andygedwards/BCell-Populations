function [parameters, varargout] = Riz2014_init_parameters_INa_low()
  % % Default parameter values for ODE model: Riz2014
  % % -----------------------------------------------
  % %
  % % parameters = Riz2014_init_parameters();
  % % [parameters, parameters_names] = Riz2014_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(86, 1);

  % --- I_Na ---
  parameters(1) = -42.0; % V_hNa;
  parameters(2) = -18; % V_mNa;
  parameters(3) = 0.4; % g_Na;
  parameters(4) = 6.0; % n_hNa;
  parameters(5) = -5.0; % n_mNa;
  parameters(6) = 2.0; % tau_hNa;

  % --- I_Na_low ---
  parameters(7) = -71.0; % V_hNa_low;
  parameters(8) = -43.0; % V_mNa_low;
  parameters(9) = 0.4;  % g_Na_low;
  parameters(10) = 2.5; % n_hNa_low;
  parameters(11) = -5.0; % n_mNa_low;
  parameters(12) = 2.0; % tau_hNa_low;

  % --- I_CaL ---
  parameters(13) = -25.0; % V_mCaL;
  parameters(14) = 0.14; % g_CaL;
  parameters(15) = -6.0; % n_mCaL;
  parameters(16) = 20.0; % tau_hCaL;

  % --- I_CaPQ ---
  parameters(17) = -10.0; % V_mCaPQ;
  parameters(18) = 0.17; % g_CaPQ;
  parameters(19) = -6.0; % n_mCaPQ;

  % --- I_CaT ---
  parameters(20) = -64.0; % V_hCaT;
  parameters(21) = -40.0; % V_mCaT;
  parameters(22) = 0.05; % g_CaT;
  parameters(23) = 8.0; % n_hCaT;
  parameters(24) = -4.0; % n_mCaT;
  parameters(25) = 7.0; % tau_hCaT;

  % --- I_SK ---
  parameters(26) = 0.57; % K_SK;
  parameters(27) = 0.1; % g_SK;
  parameters(28) = 5.2; % n;

  % --- I_BK ---
  parameters(29) = 20.0; % B_BK;
  parameters(30) = 0.0; % V_mBK;
  parameters(31) = 0.02; % g_BK;
  parameters(32) = -10.0; % n_mBK;
  parameters(33) = 2; % tau_mBK;

  % --- I_Kv ---
  parameters(34) = 0.0; % V_mKv;
  parameters(35) = 1.0; % g_Kv;
  parameters(36) = -10.0; % n_mKv;
  parameters(37) = 2.0; % tau_mKv0;

  % --- I_HERG ---
  parameters(38) = -42.0; % V_hHERG;
  parameters(39) = -30.0; % V_mHERG;
  parameters(40) = 0.0; % g_HERG;
  parameters(41) = 17.5; % n_hHERG;
  parameters(42) = -10.0; % n_mHERG;
  parameters(43) = 50.0; % tau_hHERG;
  parameters(44) = 100.0; % tau_mHERG;

  % --- I_KATP ---
  parameters(45) = 0.01; % g_KATP0;
  parameters(46) = 0.05; % g_KATP_hat;
  parameters(47) = 1.0; % glycolysis; %original value
  %parameters(47) = 1.0; % glycolysis;
  

  % --- I_GABAR ---
  parameters(48) = 0.0; % g_GABAR;

  % --- I_leak ---
  parameters(49) = -30.0; % V_leak;
  parameters(50) = 0.015; % g_leak;

  % --- J_SERCA ---
  parameters(51) = 0.06; % J_SERCA_max;
  parameters(52) = 0.27; % K_SERCA;

  % --- J_PMCA ---
  parameters(53) = 0.021; % J_PMCA_max;
  parameters(54) = 0.5; % K_PMCA;

  % --- J_leak ---
  parameters(55) = 0.00094; % J_leak;

  % --- J_NCX ---
  parameters(56) = 0.01867; % J_NCX0;

  % --- Metabolism ---
  parameters(57) = 10.0; % G;
  parameters(58) = 0.005; % K_FBA;
  parameters(59) = 0.005; % K_GAPDH;
  parameters(60) = 8.0; % K_GK;
  parameters(61) = 0.3; % K_GPI;
  parameters(62) = 4.0; % K_PFK;
  parameters(63) = 0.045455; % K_TPI;
  parameters(64) = 0.5; % P_FBA;
  parameters(65) = 0.275; % Q_FBA;
  parameters(66) = 0.000139; % V_FBA_max;
  parameters(67) = 0.00139; % V_GAPDH_max;
  parameters(68) = 5.56e-05; % V_GK_max;
  parameters(69) = 0.000556; % V_PFK_max;
  parameters(70) = 0.01; % X_PFK;
  parameters(71) = 5.0; % alpha_G;
  parameters(72) = 1.7; % h_GK;
  parameters(73) = 2.5; % h_PFK;
  parameters(74) = 2.5; % h_X;
  parameters(75) = 1.0; % h_act;
  parameters(76) = 0.0001; % k_A;

  % --- Calcium concentrations ---
  parameters(77) = 0.1; % B;
  parameters(78) = 10; % Cm;
  parameters(79) = 1.15e-12; % Vol_c;
  parameters(80) = 1e-13; % Vol_m;
  parameters(81) = 5.18e-15; % alpha;
  parameters(82) = 0.01; % f;

  % --- Membrane potential ---
  parameters(83) = 65.0; % V_Ca;
  parameters(84) = -40.0; % V_Cl;
  parameters(85) = -75.0; % V_K;
  parameters(86) = 70.0; % V_Na;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(86, 1);

    % --- I_Na ---
    parameter_names{1} = 'V_hNa';
    parameter_names{2} = 'V_mNa';
    parameter_names{3} = 'g_Na';
    parameter_names{4} = 'n_hNa';
    parameter_names{5} = 'n_mNa';
    parameter_names{6} = 'tau_hNa';

    % --- I_Na_low ---
    parameter_names{7} = 'V_hNa_low';
    parameter_names{8} = 'V_mNa_low';
    parameter_names{9} = 'g_Na_low';
    parameter_names{10} = 'n_hNa_low';
    parameter_names{11} = 'n_mNa_low';
    parameter_names{12} = 'tau_hNa_low';

    % --- I_CaL ---
    parameter_names{13} = 'V_mCaL';
    parameter_names{14} = 'g_CaL';
    parameter_names{15} = 'n_mCaL';
    parameter_names{16} = 'tau_hCaL';

    % --- I_CaPQ ---
    parameter_names{17} = 'V_mCaPQ';
    parameter_names{18} = 'g_CaPQ';
    parameter_names{19} = 'n_mCaPQ';

    % --- I_CaT ---
    parameter_names{20} = 'V_hCaT';
    parameter_names{21} = 'V_mCaT';
    parameter_names{22} = 'g_CaT';
    parameter_names{23} = 'n_hCaT';
    parameter_names{24} = 'n_mCaT';
    parameter_names{25} = 'tau_hCaT';

    % --- I_SK ---
    parameter_names{26} = 'K_SK';
    parameter_names{27} = 'g_SK';
    parameter_names{28} = 'n';

    % --- I_BK ---
    parameter_names{29} = 'B_BK';
    parameter_names{30} = 'V_mBK';
    parameter_names{31} = 'g_BK';
    parameter_names{32} = 'n_mBK';
    parameter_names{33} = 'tau_mBK';

    % --- I_Kv ---
    parameter_names{34} = 'V_mKv';
    parameter_names{35} = 'g_Kv';
    parameter_names{36} = 'n_mKv';
    parameter_names{37} = 'tau_mKv0';

    % --- I_HERG ---
    parameter_names{38} = 'V_hHERG';
    parameter_names{39} = 'V_mHERG';
    parameter_names{40} = 'g_HERG';
    parameter_names{41} = 'n_hHERG';
    parameter_names{42} = 'n_mHERG';
    parameter_names{43} = 'tau_hHERG';
    parameter_names{44} = 'tau_mHERG';

    % --- I_KATP ---
    parameter_names{45} = 'g_KATP0';
    parameter_names{46} = 'g_KATP_hat';
    parameter_names{47} = 'glycolysis';

    % --- I_GABAR ---
    parameter_names{48} = 'g_GABAR';

    % --- I_leak ---
    parameter_names{49} = 'V_leak';
    parameter_names{50} = 'g_leak';

    % --- J_SERCA ---
    parameter_names{51} = 'J_SERCA_max';
    parameter_names{52} = 'K_SERCA';

    % --- J_PMCA ---
    parameter_names{53} = 'J_PMCA_max';
    parameter_names{54} = 'K_PMCA';

    % --- J_leak ---
    parameter_names{55} = 'J_leak';

    % --- J_NCX ---
    parameter_names{56} = 'J_NCX0';

    % --- Metabolism ---
    parameter_names{57} = 'G';
    parameter_names{58} = 'K_FBA';
    parameter_names{59} = 'K_GAPDH';
    parameter_names{60} = 'K_GK';
    parameter_names{61} = 'K_GPI';
    parameter_names{62} = 'K_PFK';
    parameter_names{63} = 'K_TPI';
    parameter_names{64} = 'P_FBA';
    parameter_names{65} = 'Q_FBA';
    parameter_names{66} = 'V_FBA_max';
    parameter_names{67} = 'V_GAPDH_max';
    parameter_names{68} = 'V_GK_max';
    parameter_names{69} = 'V_PFK_max';
    parameter_names{70} = 'X_PFK';
    parameter_names{71} = 'alpha_G';
    parameter_names{72} = 'h_GK';
    parameter_names{73} = 'h_PFK';
    parameter_names{74} = 'h_X';
    parameter_names{75} = 'h_act';
    parameter_names{76} = 'k_A';

    % --- Calcium concentrations ---
    parameter_names{77} = 'B';
    parameter_names{78} = 'Cm';
    parameter_names{79} = 'Vol_c';
    parameter_names{80} = 'Vol_m';
    parameter_names{81} = 'alpha';
    parameter_names{82} = 'f';

    % --- Membrane potential ---
    parameter_names{83} = 'V_Ca';
    parameter_names{84} = 'V_Cl';
    parameter_names{85} = 'V_K';
    parameter_names{86} = 'V_Na';
    varargout(1) = {parameter_names};
  end
end