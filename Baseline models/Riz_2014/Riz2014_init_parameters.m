function [parameters, varargout] = Riz2014_init_parameters()
  % % Default parameter values for ODE model: Riz2014
  % % -----------------------------------------------
  % %
  % % parameters = Riz2014_init_parameters();
  % % [parameters, parameters_names] = Riz2014_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(80, 1);

  % --- I_Na ---
  parameters(1) = -42.0; % V_hNa;
  parameters(2) = -18; % V_mNa;
  parameters(3) = 0.4; % g_Na;
  parameters(4) = 6.0; % n_hNa;
  parameters(5) = -5.0; % n_mNa;
  parameters(6) = 2.0; % tau_hNa;

  % --- I_CaL ---
  parameters(7) = -25.0; % V_mCaL;
  parameters(8) = 0.14; % g_CaL;
  parameters(9) = -6.0; % n_mCaL;
  parameters(10) = 20.0; % tau_hCaL;

  % --- I_CaPQ ---
  parameters(11) = -10.0; % V_mCaPQ;
  parameters(12) = 0.17; % g_CaPQ;
  parameters(13) = -6.0; % n_mCaPQ;

  % --- I_CaT ---
  parameters(14) = -64.0; % V_hCaT;
  parameters(15) = -40.0; % V_mCaT;
  parameters(16) = 0.05; % g_CaT;
  parameters(17) = 8.0; % n_hCaT;
  parameters(18) = -4.0; % n_mCaT;
  parameters(19) = 7.0; % tau_hCaT;

  % --- I_SK ---
  parameters(20) = 0.57; % K_SK;
  parameters(21) = 0.1; % g_SK;
  parameters(22) = 5.2; % n;

  % --- I_BK ---
  parameters(23) = 20.0; % B_BK;
  parameters(24) = 0.0; % V_mBK;
  parameters(25) = 0.02; % g_BK;
  parameters(26) = -10.0; % n_mBK;
  parameters(27) = 2; % tau_mBK;

  % --- I_Kv ---
  parameters(28) = 0.0; % V_mKv;
  parameters(29) = 1.0; % g_Kv;
  parameters(30) = -10.0; % n_mKv;
  parameters(31) = 2.0; % tau_mKv0;

  % --- I_HERG ---
  parameters(32) = -42.0; % V_hHERG;
  parameters(33) = -30.0; % V_mHERG;
  parameters(34) = 0.0; % g_HERG;
  parameters(35) = 17.5; % n_hHERG;
  parameters(36) = -10.0; % n_mHERG;
  parameters(37) = 50.0; % tau_hHERG;
  parameters(38) = 100.0; % tau_mHERG;

  % --- I_KATP ---
  parameters(39) = 0.01; % g_KATP0;
  parameters(40) = 0.05; % g_KATP_hat;
  parameters(41) = 1; % glycolysis;

  % --- I_GABAR ---
  parameters(42) = 0.0; % g_GABAR;

  % --- I_leak ---
  parameters(43) = -30.0; % V_leak;
  parameters(44) = 0.015; % g_leak;

  % --- J_SERCA ---
  parameters(45) = 0.06; % J_SERCA_max;
  parameters(46) = 0.27; % K_SERCA;

  % --- J_PMCA ---
  parameters(47) = 0.021; % J_PMCA_max;
  parameters(48) = 0.5; % K_PMCA;

  % --- J_leak ---
  parameters(49) = 0.00094; % J_leak;

  % --- J_NCX ---
  parameters(50) = 0.01867; % J_NCX0;

  % --- Metabolism ---
  parameters(51) = 10.0; % G;
  parameters(52) = 0.005; % K_FBA;
  parameters(53) = 0.005; % K_GAPDH;
  parameters(54) = 8.0; % K_GK;
  parameters(55) = 0.3; % K_GPI;
  parameters(56) = 4.0; % K_PFK;
  parameters(57) = 0.045455; % K_TPI;
  parameters(58) = 0.5; % P_FBA;
  parameters(59) = 0.275; % Q_FBA;
  parameters(60) = 0.000139; % V_FBA_max;
  parameters(61) = 0.00139; % V_GAPDH_max;
  parameters(62) = 5.56e-05; % V_GK_max;
  parameters(63) = 0.000556; % V_PFK_max;
  parameters(64) = 0.01; % X_PFK;
  parameters(65) = 5.0; % alpha_G;
  parameters(66) = 1.7; % h_GK;
  parameters(67) = 2.5; % h_PFK;
  parameters(68) = 2.5; % h_X;
  parameters(69) = 1.0; % h_act;
  parameters(70) = 0.0001; % k_A;

  % --- Calcium concentrations ---
  parameters(71) = 0.1; % B;
  parameters(72) = 10; % Cm;
  parameters(73) = 1.15e-12; % Vol_c;
  parameters(74) = 1e-13; % Vol_m;
  parameters(75) = 5.18e-15; % alpha;
  parameters(76) = 0.01; % f;

  % --- Membrane potential ---
  parameters(77) = 65.0; % V_Ca;
  parameters(78) = -40.0; % V_Cl;
  parameters(79) = -75.0; % V_K;
  parameters(80) = 70.0; % V_Na;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(80, 1);

    % --- I_Na ---
    parameter_names{1} = 'V_hNa';
    parameter_names{2} = 'V_mNa';
    parameter_names{3} = 'g_Na';
    parameter_names{4} = 'n_hNa';
    parameter_names{5} = 'n_mNa';
    parameter_names{6} = 'tau_hNa';

    % --- I_CaL ---
    parameter_names{7} = 'V_mCaL';
    parameter_names{8} = 'g_CaL';
    parameter_names{9} = 'n_mCaL';
    parameter_names{10} = 'tau_hCaL';

    % --- I_CaPQ ---
    parameter_names{11} = 'V_mCaPQ';
    parameter_names{12} = 'g_CaPQ';
    parameter_names{13} = 'n_mCaPQ';

    % --- I_CaT ---
    parameter_names{14} = 'V_hCaT';
    parameter_names{15} = 'V_mCaT';
    parameter_names{16} = 'g_CaT';
    parameter_names{17} = 'n_hCaT';
    parameter_names{18} = 'n_mCaT';
    parameter_names{19} = 'tau_hCaT';

    % --- I_SK ---
    parameter_names{20} = 'K_SK';
    parameter_names{21} = 'g_SK';
    parameter_names{22} = 'n';

    % --- I_BK ---
    parameter_names{23} = 'B_BK';
    parameter_names{24} = 'V_mBK';
    parameter_names{25} = 'g_BK';
    parameter_names{26} = 'n_mBK';
    parameter_names{27} = 'tau_mBK';

    % --- I_Kv ---
    parameter_names{28} = 'V_mKv';
    parameter_names{29} = 'g_Kv';
    parameter_names{30} = 'n_mKv';
    parameter_names{31} = 'tau_mKv0';

    % --- I_HERG ---
    parameter_names{32} = 'V_hHERG';
    parameter_names{33} = 'V_mHERG';
    parameter_names{34} = 'g_HERG';
    parameter_names{35} = 'n_hHERG';
    parameter_names{36} = 'n_mHERG';
    parameter_names{37} = 'tau_hHERG';
    parameter_names{38} = 'tau_mHERG';

    % --- I_KATP ---
    parameter_names{39} = 'g_KATP0';
    parameter_names{40} = 'g_KATP_hat';
    parameter_names{41} = 'glycolysis';

    % --- I_GABAR ---
    parameter_names{42} = 'g_GABAR';

    % --- I_leak ---
    parameter_names{43} = 'V_leak';
    parameter_names{44} = 'g_leak';

    % --- J_SERCA ---
    parameter_names{45} = 'J_SERCA_max';
    parameter_names{46} = 'K_SERCA';

    % --- J_PMCA ---
    parameter_names{47} = 'J_PMCA_max';
    parameter_names{48} = 'K_PMCA';

    % --- J_leak ---
    parameter_names{49} = 'J_leak';

    % --- J_NCX ---
    parameter_names{50} = 'J_NCX0';

    % --- Metabolism ---
    parameter_names{51} = 'G';
    parameter_names{52} = 'K_FBA';
    parameter_names{53} = 'K_GAPDH';
    parameter_names{54} = 'K_GK';
    parameter_names{55} = 'K_GPI';
    parameter_names{56} = 'K_PFK';
    parameter_names{57} = 'K_TPI';
    parameter_names{58} = 'P_FBA';
    parameter_names{59} = 'Q_FBA';
    parameter_names{60} = 'V_FBA_max';
    parameter_names{61} = 'V_GAPDH_max';
    parameter_names{62} = 'V_GK_max';
    parameter_names{63} = 'V_PFK_max';
    parameter_names{64} = 'X_PFK';
    parameter_names{65} = 'alpha_G';
    parameter_names{66} = 'h_GK';
    parameter_names{67} = 'h_PFK';
    parameter_names{68} = 'h_X';
    parameter_names{69} = 'h_act';
    parameter_names{70} = 'k_A';

    % --- Calcium concentrations ---
    parameter_names{71} = 'B';
    parameter_names{72} = 'Cm';
    parameter_names{73} = 'Vol_c';
    parameter_names{74} = 'Vol_m';
    parameter_names{75} = 'alpha';
    parameter_names{76} = 'f';

    % --- Membrane potential ---
    parameter_names{77} = 'V_Ca';
    parameter_names{78} = 'V_Cl';
    parameter_names{79} = 'V_K';
    parameter_names{80} = 'V_Na';
    varargout(1) = {parameter_names};
  end
end
