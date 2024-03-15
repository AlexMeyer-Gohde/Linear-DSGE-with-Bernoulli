%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'US_HL16_NEP';
M_.dynare_version = '4.5.7';
oo_.dynare_version = '4.5.7';
options_.dynare_version = '4.5.7';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('US_HL16_NEP.log');
M_.exo_names = 'epsilon_p';
M_.exo_names_tex = 'epsilon\_p';
M_.exo_names_long = 'epsilon_p';
M_.exo_names = char(M_.exo_names, 'epsilon_z');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_z');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_z');
M_.exo_names = char(M_.exo_names, 'epsilon_i');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_i');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_i');
M_.exo_names = char(M_.exo_names, 'epsilon_d');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_d');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_d');
M_.exo_names = char(M_.exo_names, 'epsilon_h');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_h');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_h');
M_.exo_names = char(M_.exo_names, 'epsilon_e');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_e');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_e');
M_.exo_names = char(M_.exo_names, 'epsilon_nu_h');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_nu\_h');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_nu_h');
M_.exo_names = char(M_.exo_names, 'epsilon_nu_e');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_nu\_e');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_nu_e');
M_.endo_names = 'c_b';
M_.endo_names_tex = 'c\_b';
M_.endo_names_long = 'c_b';
M_.endo_names = char(M_.endo_names, 'i_d');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_d');
M_.endo_names_long = char(M_.endo_names_long, 'i_d');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'i_h');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_h');
M_.endo_names_long = char(M_.endo_names_long, 'i_h');
M_.endo_names = char(M_.endo_names, 'lambda');
M_.endo_names_tex = char(M_.endo_names_tex, 'lambda');
M_.endo_names_long = char(M_.endo_names_long, 'lambda');
M_.endo_names = char(M_.endo_names, 'l_h');
M_.endo_names_tex = char(M_.endo_names_tex, 'l\_h');
M_.endo_names_long = char(M_.endo_names_long, 'l_h');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'nu_h');
M_.endo_names_tex = char(M_.endo_names_tex, 'nu\_h');
M_.endo_names_long = char(M_.endo_names_long, 'nu_h');
M_.endo_names = char(M_.endo_names, 'h_s');
M_.endo_names_tex = char(M_.endo_names_tex, 'h\_s');
M_.endo_names_long = char(M_.endo_names_long, 'h_s');
M_.endo_names = char(M_.endo_names, 'h_b');
M_.endo_names_tex = char(M_.endo_names_tex, 'h\_b');
M_.endo_names_long = char(M_.endo_names_long, 'h_b');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.endo_names = char(M_.endo_names, 'c_s');
M_.endo_names_tex = char(M_.endo_names_tex, 'c\_s');
M_.endo_names_long = char(M_.endo_names_long, 'c_s');
M_.endo_names = char(M_.endo_names, 'x');
M_.endo_names_tex = char(M_.endo_names_tex, 'x');
M_.endo_names_long = char(M_.endo_names_long, 'x');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'k_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'k\_e');
M_.endo_names_long = char(M_.endo_names_long, 'k_e');
M_.endo_names = char(M_.endo_names, 'h');
M_.endo_names_tex = char(M_.endo_names_tex, 'h');
M_.endo_names_long = char(M_.endo_names_long, 'h');
M_.endo_names = char(M_.endo_names, 'l_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'l\_e');
M_.endo_names_long = char(M_.endo_names_long, 'l_e');
M_.endo_names = char(M_.endo_names, 'q_k');
M_.endo_names_tex = char(M_.endo_names_tex, 'q\_k');
M_.endo_names_long = char(M_.endo_names_long, 'q_k');
M_.endo_names = char(M_.endo_names, 'i_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_e');
M_.endo_names_long = char(M_.endo_names_long, 'i_e');
M_.endo_names = char(M_.endo_names, 'nu_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'nu\_e');
M_.endo_names_long = char(M_.endo_names_long, 'nu_e');
M_.endo_names = char(M_.endo_names, 'lambda_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'lambda\_e');
M_.endo_names_long = char(M_.endo_names_long, 'lambda_e');
M_.endo_names = char(M_.endo_names, 'v');
M_.endo_names_tex = char(M_.endo_names_tex, 'v');
M_.endo_names_long = char(M_.endo_names_long, 'v');
M_.endo_names = char(M_.endo_names, 'k_B');
M_.endo_names_tex = char(M_.endo_names_tex, 'k\_B');
M_.endo_names_long = char(M_.endo_names_long, 'k_B');
M_.endo_names = char(M_.endo_names, 'pi_B');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_B');
M_.endo_names_long = char(M_.endo_names_long, 'pi_B');
M_.endo_names = char(M_.endo_names, 'i_l');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_l');
M_.endo_names_long = char(M_.endo_names_long, 'i_l');
M_.endo_names = char(M_.endo_names, 'l');
M_.endo_names_tex = char(M_.endo_names_tex, 'l');
M_.endo_names_long = char(M_.endo_names_long, 'l');
M_.endo_names = char(M_.endo_names, 'k_BL');
M_.endo_names_tex = char(M_.endo_names_tex, 'k\_BL');
M_.endo_names_long = char(M_.endo_names_long, 'k_BL');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'xi_z');
M_.endo_names_tex = char(M_.endo_names_tex, 'xi\_z');
M_.endo_names_long = char(M_.endo_names_long, 'xi_z');
M_.endo_names = char(M_.endo_names, 'xi_d');
M_.endo_names_tex = char(M_.endo_names_tex, 'xi\_d');
M_.endo_names_long = char(M_.endo_names_long, 'xi_d');
M_.endo_names = char(M_.endo_names, 'xi_i');
M_.endo_names_tex = char(M_.endo_names_tex, 'xi\_i');
M_.endo_names_long = char(M_.endo_names_long, 'xi_i');
M_.endo_names = char(M_.endo_names, 'varepsilon_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'varepsilon\_e');
M_.endo_names_long = char(M_.endo_names_long, 'varepsilon_e');
M_.endo_names = char(M_.endo_names, 'varepsilon_h');
M_.endo_names_tex = char(M_.endo_names_tex, 'varepsilon\_h');
M_.endo_names_long = char(M_.endo_names_long, 'varepsilon_h');
M_.endo_names = char(M_.endo_names, 'xi_p');
M_.endo_names_tex = char(M_.endo_names_tex, 'xi\_p');
M_.endo_names_long = char(M_.endo_names_long, 'xi_p');
M_.endo_names = char(M_.endo_names, 'pi_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'pi_obs');
M_.endo_names = char(M_.endo_names, 'y_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'y\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'y_obs');
M_.endo_names = char(M_.endo_names, 'l_h_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'l\_h\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'l_h_obs');
M_.endo_names = char(M_.endo_names, 'l_e_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'l\_e\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'l_e_obs');
M_.endo_names = char(M_.endo_names, 'd_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'd\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'd_obs');
M_.endo_names = char(M_.endo_names, 'i_d_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_d\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'i_d_obs');
M_.endo_names = char(M_.endo_names, 'i_h_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_h\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'i_h_obs');
M_.endo_names = char(M_.endo_names, 'i_e_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_e\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'i_e_obs');
M_.endo_partitions = struct();
M_.param_names = 'beta_b';
M_.param_names_tex = 'beta\_b';
M_.param_names_long = 'beta_b';
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'R_d');
M_.param_names_tex = char(M_.param_names_tex, 'R\_d');
M_.param_names_long = char(M_.param_names_long, 'R_d');
M_.param_names = char(M_.param_names, 'R_h');
M_.param_names_tex = char(M_.param_names_tex, 'R\_h');
M_.param_names_long = char(M_.param_names_long, 'R_h');
M_.param_names = char(M_.param_names, 'Nu_h');
M_.param_names_tex = char(M_.param_names_tex, 'Nu\_h');
M_.param_names_long = char(M_.param_names_long, 'Nu_h');
M_.param_names = char(M_.param_names, 'phi_w');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_w');
M_.param_names_long = char(M_.param_names_long, 'phi_w');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'gamma_b');
M_.param_names_tex = char(M_.param_names_tex, 'gamma\_b');
M_.param_names_long = char(M_.param_names_long, 'gamma_b');
M_.param_names = char(M_.param_names, 'beta_s');
M_.param_names_tex = char(M_.param_names_tex, 'beta\_s');
M_.param_names_long = char(M_.param_names_long, 'beta_s');
M_.param_names = char(M_.param_names, 'theta_R');
M_.param_names_tex = char(M_.param_names_tex, 'theta\_R');
M_.param_names_long = char(M_.param_names_long, 'theta_R');
M_.param_names = char(M_.param_names, 'varepsilon_p');
M_.param_names_tex = char(M_.param_names_tex, 'varepsilon\_p');
M_.param_names_long = char(M_.param_names_long, 'varepsilon_p');
M_.param_names = char(M_.param_names, 'gamma_p');
M_.param_names_tex = char(M_.param_names_tex, 'gamma\_p');
M_.param_names_long = char(M_.param_names_long, 'gamma_p');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'Nu_e');
M_.param_names_tex = char(M_.param_names_tex, 'Nu\_e');
M_.param_names_long = char(M_.param_names_long, 'Nu_e');
M_.param_names = char(M_.param_names, 'phi_k');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_k');
M_.param_names_long = char(M_.param_names_long, 'phi_k');
M_.param_names = char(M_.param_names, 'R_e');
M_.param_names_tex = char(M_.param_names_tex, 'R\_e');
M_.param_names_long = char(M_.param_names_long, 'R_e');
M_.param_names = char(M_.param_names, 'beta_e');
M_.param_names_tex = char(M_.param_names_tex, 'beta\_e');
M_.param_names_long = char(M_.param_names_long, 'beta_e');
M_.param_names = char(M_.param_names, 'delta_e');
M_.param_names_tex = char(M_.param_names_tex, 'delta\_e');
M_.param_names_long = char(M_.param_names_long, 'delta_e');
M_.param_names = char(M_.param_names, 'kappa_v');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_v');
M_.param_names_long = char(M_.param_names_long, 'kappa_v');
M_.param_names = char(M_.param_names, 'phi_s');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_s');
M_.param_names_long = char(M_.param_names_long, 'phi_s');
M_.param_names = char(M_.param_names, 'delta_B');
M_.param_names_tex = char(M_.param_names_tex, 'delta\_B');
M_.param_names_long = char(M_.param_names_long, 'delta_B');
M_.param_names = char(M_.param_names, 'kappa_k');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_k');
M_.param_names_long = char(M_.param_names_long, 'kappa_k');
M_.param_names = char(M_.param_names, 'tau');
M_.param_names_tex = char(M_.param_names_tex, 'tau');
M_.param_names_long = char(M_.param_names_long, 'tau');
M_.param_names = char(M_.param_names, 'kappa_e');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_e');
M_.param_names_long = char(M_.param_names_long, 'kappa_e');
M_.param_names = char(M_.param_names, 'varepsilon_ess');
M_.param_names_tex = char(M_.param_names_tex, 'varepsilon\_ess');
M_.param_names_long = char(M_.param_names_long, 'varepsilon_ess');
M_.param_names = char(M_.param_names, 'beta_B');
M_.param_names_tex = char(M_.param_names_tex, 'beta\_B');
M_.param_names_long = char(M_.param_names_long, 'beta_B');
M_.param_names = char(M_.param_names, 'kappa_h');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_h');
M_.param_names_long = char(M_.param_names_long, 'kappa_h');
M_.param_names = char(M_.param_names, 'varepsilon_hss');
M_.param_names_tex = char(M_.param_names_tex, 'varepsilon\_hss');
M_.param_names_long = char(M_.param_names_long, 'varepsilon_hss');
M_.param_names = char(M_.param_names, 'L_hL');
M_.param_names_tex = char(M_.param_names_tex, 'L\_hL');
M_.param_names_long = char(M_.param_names_long, 'L_hL');
M_.param_names = char(M_.param_names, 'L_eL');
M_.param_names_tex = char(M_.param_names_tex, 'L\_eL');
M_.param_names_long = char(M_.param_names_long, 'L_eL');
M_.param_names = char(M_.param_names, 'kappa_i');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_i');
M_.param_names_long = char(M_.param_names_long, 'kappa_i');
M_.param_names = char(M_.param_names, 'kappa_pi');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_pi');
M_.param_names_long = char(M_.param_names_long, 'kappa_pi');
M_.param_names = char(M_.param_names, 'kappa_y');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_y');
M_.param_names_long = char(M_.param_names_long, 'kappa_y');
M_.param_names = char(M_.param_names, 'CY');
M_.param_names_tex = char(M_.param_names_tex, 'CY');
M_.param_names_long = char(M_.param_names_long, 'CY');
M_.param_names = char(M_.param_names, 'K_BY');
M_.param_names_tex = char(M_.param_names_tex, 'K\_BY');
M_.param_names_long = char(M_.param_names_long, 'K_BY');
M_.param_names = char(M_.param_names, 'LY');
M_.param_names_tex = char(M_.param_names_tex, 'LY');
M_.param_names_long = char(M_.param_names_long, 'LY');
M_.param_names = char(M_.param_names, 'rho_z');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_z');
M_.param_names_long = char(M_.param_names_long, 'rho_z');
M_.param_names = char(M_.param_names, 'rho_d');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_d');
M_.param_names_long = char(M_.param_names_long, 'rho_d');
M_.param_names = char(M_.param_names, 'rho_i');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_i');
M_.param_names_long = char(M_.param_names_long, 'rho_i');
M_.param_names = char(M_.param_names, 'rho_e');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_e');
M_.param_names_long = char(M_.param_names_long, 'rho_e');
M_.param_names = char(M_.param_names, 'rho_h');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_h');
M_.param_names_long = char(M_.param_names_long, 'rho_h');
M_.param_names = char(M_.param_names, 'rho_nuh');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_nuh');
M_.param_names_long = char(M_.param_names_long, 'rho_nuh');
M_.param_names = char(M_.param_names, 'rho_nue');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_nue');
M_.param_names_long = char(M_.param_names_long, 'rho_nue');
M_.param_names = char(M_.param_names, 'rho_p');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_p');
M_.param_names_long = char(M_.param_names_long, 'rho_p');
M_.param_names = char(M_.param_names, 'i_h_ss');
M_.param_names_tex = char(M_.param_names_tex, 'i\_h\_ss');
M_.param_names_long = char(M_.param_names_long, 'i_h_ss');
M_.param_names = char(M_.param_names, 'i_e_ss');
M_.param_names_tex = char(M_.param_names_tex, 'i\_e\_ss');
M_.param_names_long = char(M_.param_names_long, 'i_e_ss');
M_.param_names = char(M_.param_names, 'i_d_ss');
M_.param_names_tex = char(M_.param_names_tex, 'i\_d\_ss');
M_.param_names_long = char(M_.param_names_long, 'i_d_ss');
M_.param_names = char(M_.param_names, 'pi_ss');
M_.param_names_tex = char(M_.param_names_tex, 'pi\_ss');
M_.param_names_long = char(M_.param_names_long, 'pi_ss');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 8;
M_.endo_nbr = 42;
M_.param_nbr = 49;
M_.orig_endo_nbr = 42;
M_.aux_vars = [];
M_.Sigma_e = zeros(8, 8);
M_.Correlation_matrix = eye(8, 8);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('US_HL16_NEP_static');
erase_compiled_function('US_HL16_NEP_dynamic');
M_.orig_eq_nbr = 42;
M_.eq_nbr = 42;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 22 64;
 2 23 0;
 3 24 65;
 4 25 66;
 0 26 0;
 5 27 0;
 0 28 67;
 6 29 0;
 0 30 0;
 0 31 0;
 7 32 0;
 8 33 68;
 0 34 69;
 9 35 70;
 10 36 0;
 0 37 0;
 11 38 0;
 0 39 71;
 12 40 72;
 13 41 0;
 0 42 0;
 0 43 73;
 14 44 0;
 15 45 0;
 0 46 0;
 0 47 0;
 0 48 0;
 0 49 0;
 16 50 0;
 17 51 0;
 18 52 0;
 19 53 0;
 20 54 0;
 21 55 0;
 0 56 0;
 0 57 0;
 0 58 0;
 0 59 0;
 0 60 0;
 0 61 0;
 0 62 0;
 0 63 0;]';
M_.nstatic = 17;
M_.nfwrd   = 4;
M_.npred   = 15;
M_.nboth   = 6;
M_.nsfwrd   = 10;
M_.nspred   = 21;
M_.ndynamic   = 25;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:8];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(42, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(8, 1);
M_.params = NaN(49, 1);
M_.NNZDerivatives = [170; -1; -1];
M_.params( 1 ) = 0.96;
beta_b = M_.params( 1 );
M_.params( 2 ) = 0.5;
phi = M_.params( 2 );
M_.params( 3 ) = 1.01;
R_d = M_.params( 3 );
M_.params( 4 ) = 1.033;
R_h = M_.params( 4 );
M_.params( 5 ) = 0.75;
Nu_h = M_.params( 5 );
M_.params( 6 ) = 1;
phi_w = M_.params( 6 );
M_.params( 7 ) = 1;
eta = M_.params( 7 );
M_.params( 9 ) = 2;
gamma_b = M_.params( 9 );
M_.params( 10 ) = 0.99;
beta_s = M_.params( 10 );
M_.params( 8 ) = 2;
gamma = M_.params( 8 );
M_.params( 18 ) = 0.95;
beta_e = M_.params( 18 );
M_.params( 17 ) = 1.039;
R_e = M_.params( 17 );
M_.params( 15 ) = 0.65;
Nu_e = M_.params( 15 );
M_.params( 16 ) = 1;
phi_k = M_.params( 16 );
M_.params( 14 ) = 0.33;
alpha = M_.params( 14 );
M_.params( 19 ) = 0.025;
delta_e = M_.params( 19 );
M_.params( 12 ) = 11;
varepsilon_p = M_.params( 12 );
M_.params( 20 ) = 2;
kappa_v = M_.params( 20 );
M_.params( 11 ) = 0.65;
theta_R = M_.params( 11 );
M_.params( 13 ) = 0.25;
gamma_p = M_.params( 13 );
M_.params( 21 ) = 0.53;
phi_s = M_.params( 21 );
M_.params( 27 ) = 0.99;
beta_B = M_.params( 27 );
M_.params( 22 ) = 0.4;
delta_B = M_.params( 22 );
M_.params( 23 ) = 4;
kappa_k = M_.params( 23 );
M_.params( 25 ) = 8;
kappa_e = M_.params( 25 );
M_.params( 28 ) = 8;
kappa_h = M_.params( 28 );
M_.params( 24 ) = 0.11;
tau = M_.params( 24 );
M_.params( 26 ) = 1.341;
varepsilon_ess = M_.params( 26 );
M_.params( 29 ) = 1.427;
varepsilon_hss = M_.params( 29 );
M_.params( 30 ) = 0.45;
L_hL = M_.params( 30 );
M_.params( 31 ) = 0.55;
L_eL = M_.params( 31 );
M_.params( 32 ) = 0.65;
kappa_i = M_.params( 32 );
M_.params( 33 ) = 1.5;
kappa_pi = M_.params( 33 );
M_.params( 34 ) = 0.25;
kappa_y = M_.params( 34 );
M_.params( 35 ) = 0.653;
CY = M_.params( 35 );
M_.params( 36 ) = 0.165;
K_BY = M_.params( 36 );
M_.params( 37 ) = 1.5;
LY = M_.params( 37 );
M_.params( 38 ) = 0.75;
rho_z = M_.params( 38 );
M_.params( 39 ) = 0.5;
rho_d = M_.params( 39 );
M_.params( 40 ) = 0.5;
rho_i = M_.params( 40 );
M_.params( 41 ) = 0.5;
rho_e = M_.params( 41 );
M_.params( 42 ) = 0.5;
rho_h = M_.params( 42 );
M_.params( 43 ) = 0.75;
rho_nuh = M_.params( 43 );
M_.params( 44 ) = 0.75;
rho_nue = M_.params( 44 );
M_.params( 45 ) = 0.5;
rho_p = M_.params( 45 );
M_.params( 49 ) = 0.005887;
pi_ss = M_.params( 49 );
M_.params( 46 ) = 0.078998;
i_h_ss = M_.params( 46 );
M_.params( 47 ) = 0.08487;
i_e_ss = M_.params( 47 );
M_.params( 48 ) = 0.04558;
i_d_ss = M_.params( 48 );
load median_param_NEP.mat;
coeffs = median_param;
M_.params( 8 ) = coeffs(1);
gamma = M_.params( 8 );
M_.params( 9 ) = coeffs(2);
gamma_b = M_.params( 9 );
M_.params( 2 ) = coeffs(3);
phi = M_.params( 2 );
M_.params( 11 ) = coeffs(4);
theta_R = M_.params( 11 );
M_.params( 13 ) = coeffs(5);
gamma_p = M_.params( 13 );
M_.params( 5 ) = coeffs(7);
Nu_h = M_.params( 5 );
M_.params( 15 ) = coeffs(9);
Nu_e = M_.params( 15 );
M_.params( 23 ) = coeffs(11);
kappa_k = M_.params( 23 );
M_.params( 25 ) = coeffs(12);
kappa_e = M_.params( 25 );
M_.params( 28 ) = coeffs(13);
kappa_h = M_.params( 28 );
M_.params( 32 ) = coeffs(14);
kappa_i = M_.params( 32 );
M_.params( 33 ) = coeffs(15);
kappa_pi = M_.params( 33 );
M_.params( 34 ) = coeffs(16);
kappa_y = M_.params( 34 );
M_.params( 38 ) = coeffs(17);
rho_z = M_.params( 38 );
M_.params( 39 ) = coeffs(18);
rho_d = M_.params( 39 );
M_.params( 40 ) = coeffs(19);
rho_i = M_.params( 40 );
M_.params( 41 ) = coeffs(20);
rho_e = M_.params( 41 );
M_.params( 42 ) = coeffs(21);
rho_h = M_.params( 42 );
M_.params( 43 ) = coeffs(22);
rho_nuh = M_.params( 43 );
M_.params( 44 ) = coeffs(23);
rho_nue = M_.params( 44 );
M_.params( 45 ) = coeffs(25);
rho_p = M_.params( 45 );
load median_std_NEP.mat;
coeffs_std = median_std;
std_z	=	coeffs_std(1);
std_i	=	coeffs_std(2);
std_d	=	coeffs_std(3);
std_e	=	coeffs_std(4);
std_h	=	coeffs_std(5);
std_nu_h	=	coeffs_std(6);
std_nu_e	=	coeffs_std(7);
std_p	=	coeffs_std(9);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = std_p^2;
M_.Sigma_e(2, 2) = std_z^2;
M_.Sigma_e(3, 3) = std_i^2;
M_.Sigma_e(4, 4) = std_d^2;
M_.Sigma_e(5, 5) = std_h^2;
M_.Sigma_e(6, 6) = std_e^2;
M_.Sigma_e(7, 7) = std_nu_h^2;
M_.Sigma_e(8, 8) = std_nu_e^2;
options_.irf = 15;
options_.nograph = 1;
options_.order = 1;
var_list_ = char();
info = stoch_simul(var_list_);
save('US_HL16_NEP_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('US_HL16_NEP_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('US_HL16_NEP_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('US_HL16_NEP_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('US_HL16_NEP_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('US_HL16_NEP_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('US_HL16_NEP_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off