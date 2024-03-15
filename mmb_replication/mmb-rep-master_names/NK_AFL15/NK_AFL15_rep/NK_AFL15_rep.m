%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'NK_AFL15_rep';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('NK_AFL15_rep.log');
M_.exo_names = 'ua';
M_.exo_names_tex = 'ua';
M_.exo_names_long = 'ua';
M_.exo_names = char(M_.exo_names, 'ug');
M_.exo_names_tex = char(M_.exo_names_tex, 'ug');
M_.exo_names_long = char(M_.exo_names_long, 'ug');
M_.exo_names = char(M_.exo_names, 'ur');
M_.exo_names_tex = char(M_.exo_names_tex, 'ur');
M_.exo_names_long = char(M_.exo_names_long, 'ur');
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names_long = 'c';
M_.endo_names = char(M_.endo_names, 'pai');
M_.endo_names_tex = char(M_.endo_names_tex, 'pai');
M_.endo_names_long = char(M_.endo_names_long, 'pai');
M_.endo_names = char(M_.endo_names, 'rn');
M_.endo_names_tex = char(M_.endo_names_tex, 'rn');
M_.endo_names_long = char(M_.endo_names_long, 'rn');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'mc');
M_.endo_names_tex = char(M_.endo_names_tex, 'mc');
M_.endo_names_long = char(M_.endo_names_long, 'mc');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'inv');
M_.endo_names_tex = char(M_.endo_names_tex, 'inv');
M_.endo_names_long = char(M_.endo_names_long, 'inv');
M_.endo_names = char(M_.endo_names, 'rok');
M_.endo_names_tex = char(M_.endo_names_tex, 'rok');
M_.endo_names_long = char(M_.endo_names_long, 'rok');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'g');
M_.endo_names_tex = char(M_.endo_names_tex, 'g');
M_.endo_names_long = char(M_.endo_names_long, 'g');
M_.endo_names = char(M_.endo_names, 'uc');
M_.endo_names_tex = char(M_.endo_names_tex, 'uc');
M_.endo_names_long = char(M_.endo_names_long, 'uc');
M_.endo_names = char(M_.endo_names, 'un');
M_.endo_names_tex = char(M_.endo_names_tex, 'un');
M_.endo_names_long = char(M_.endo_names_long, 'un');
M_.endo_names = char(M_.endo_names, 'Fk');
M_.endo_names_tex = char(M_.endo_names_tex, 'Fk');
M_.endo_names_long = char(M_.endo_names_long, 'Fk');
M_.endo_names = char(M_.endo_names, 'Fn');
M_.endo_names_tex = char(M_.endo_names_tex, 'Fn');
M_.endo_names_long = char(M_.endo_names_long, 'Fn');
M_.endo_names = char(M_.endo_names, 'bk');
M_.endo_names_tex = char(M_.endo_names_tex, 'bk');
M_.endo_names_long = char(M_.endo_names_long, 'bk');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.endo_names = char(M_.endo_names, 'deprat');
M_.endo_names_tex = char(M_.endo_names_tex, 'deprat');
M_.endo_names_long = char(M_.endo_names_long, 'deprat');
M_.endo_names = char(M_.endo_names, 'ra');
M_.endo_names_tex = char(M_.endo_names_tex, 'ra');
M_.endo_names_long = char(M_.endo_names_long, 'ra');
M_.endo_names = char(M_.endo_names, 'br');
M_.endo_names_tex = char(M_.endo_names_tex, 'br');
M_.endo_names_long = char(M_.endo_names_long, 'br');
M_.endo_names = char(M_.endo_names, 'fai');
M_.endo_names_tex = char(M_.endo_names_tex, 'fai');
M_.endo_names_long = char(M_.endo_names_long, 'fai');
M_.endo_names = char(M_.endo_names, 'rd');
M_.endo_names_tex = char(M_.endo_names_tex, 'rd');
M_.endo_names_long = char(M_.endo_names_long, 'rd');
M_.endo_names = char(M_.endo_names, 'crun');
M_.endo_names_tex = char(M_.endo_names_tex, 'crun');
M_.endo_names_long = char(M_.endo_names_long, 'crun');
M_.endo_names = char(M_.endo_names, 'cpai');
M_.endo_names_tex = char(M_.endo_names_tex, 'cpai');
M_.endo_names_long = char(M_.endo_names_long, 'cpai');
M_.endo_names = char(M_.endo_names, 'rsh');
M_.endo_names_tex = char(M_.endo_names_tex, 'rsh');
M_.endo_names_long = char(M_.endo_names_long, 'rsh');
M_.param_names = 'PSI';
M_.param_names_tex = 'PSI';
M_.param_names_long = 'PSI';
M_.param_names = char(M_.param_names, 'calvo');
M_.param_names_tex = char(M_.param_names_tex, 'calvo');
M_.param_names_long = char(M_.param_names_long, 'calvo');
M_.param_names = char(M_.param_names, 'ritcapss');
M_.param_names_tex = char(M_.param_names_tex, 'ritcapss');
M_.param_names_long = char(M_.param_names_long, 'ritcapss');
M_.param_names = char(M_.param_names, 'YoK');
M_.param_names_tex = char(M_.param_names_tex, 'YoK');
M_.param_names_long = char(M_.param_names_long, 'YoK');
M_.param_names = char(M_.param_names, 'IoK');
M_.param_names_tex = char(M_.param_names_tex, 'IoK');
M_.param_names_long = char(M_.param_names_long, 'IoK');
M_.param_names = char(M_.param_names, 'CoY');
M_.param_names_tex = char(M_.param_names_tex, 'CoY');
M_.param_names_long = char(M_.param_names_long, 'CoY');
M_.param_names = char(M_.param_names, 'PAIss');
M_.param_names_tex = char(M_.param_names_tex, 'PAIss');
M_.param_names_long = char(M_.param_names_long, 'PAIss');
M_.param_names = char(M_.param_names, 'css');
M_.param_names_tex = char(M_.param_names_tex, 'css');
M_.param_names_long = char(M_.param_names_long, 'css');
M_.param_names = char(M_.param_names, 'rnss');
M_.param_names_tex = char(M_.param_names_tex, 'rnss');
M_.param_names_long = char(M_.param_names_long, 'rnss');
M_.param_names = char(M_.param_names, 'kss');
M_.param_names_tex = char(M_.param_names_tex, 'kss');
M_.param_names_long = char(M_.param_names_long, 'kss');
M_.param_names = char(M_.param_names, 'zss');
M_.param_names_tex = char(M_.param_names_tex, 'zss');
M_.param_names_long = char(M_.param_names_long, 'zss');
M_.param_names = char(M_.param_names, 'yss');
M_.param_names_tex = char(M_.param_names_tex, 'yss');
M_.param_names_long = char(M_.param_names_long, 'yss');
M_.param_names = char(M_.param_names, 'mcss');
M_.param_names_tex = char(M_.param_names_tex, 'mcss');
M_.param_names_long = char(M_.param_names_long, 'mcss');
M_.param_names = char(M_.param_names, 'nss');
M_.param_names_tex = char(M_.param_names_tex, 'nss');
M_.param_names_long = char(M_.param_names_long, 'nss');
M_.param_names = char(M_.param_names, 'qss');
M_.param_names_tex = char(M_.param_names_tex, 'qss');
M_.param_names_long = char(M_.param_names_long, 'qss');
M_.param_names = char(M_.param_names, 'invss');
M_.param_names_tex = char(M_.param_names_tex, 'invss');
M_.param_names_long = char(M_.param_names_long, 'invss');
M_.param_names = char(M_.param_names, 'rokss');
M_.param_names_tex = char(M_.param_names_tex, 'rokss');
M_.param_names_long = char(M_.param_names_long, 'rokss');
M_.param_names = char(M_.param_names, 'gss');
M_.param_names_tex = char(M_.param_names_tex, 'gss');
M_.param_names_long = char(M_.param_names_long, 'gss');
M_.param_names = char(M_.param_names, 'ucss');
M_.param_names_tex = char(M_.param_names_tex, 'ucss');
M_.param_names_long = char(M_.param_names_long, 'ucss');
M_.param_names = char(M_.param_names, 'unss');
M_.param_names_tex = char(M_.param_names_tex, 'unss');
M_.param_names_long = char(M_.param_names_long, 'unss');
M_.param_names = char(M_.param_names, 'Fkss');
M_.param_names_tex = char(M_.param_names_tex, 'Fkss');
M_.param_names_long = char(M_.param_names_long, 'Fkss');
M_.param_names = char(M_.param_names, 'Fnss');
M_.param_names_tex = char(M_.param_names_tex, 'Fnss');
M_.param_names_long = char(M_.param_names_long, 'Fnss');
M_.param_names = char(M_.param_names, 'bkss');
M_.param_names_tex = char(M_.param_names_tex, 'bkss');
M_.param_names_long = char(M_.param_names_long, 'bkss');
M_.param_names = char(M_.param_names, 'dss');
M_.param_names_tex = char(M_.param_names_tex, 'dss');
M_.param_names_long = char(M_.param_names_long, 'dss');
M_.param_names = char(M_.param_names, 'depratss');
M_.param_names_tex = char(M_.param_names_tex, 'depratss');
M_.param_names_long = char(M_.param_names_long, 'depratss');
M_.param_names = char(M_.param_names, 'rass');
M_.param_names_tex = char(M_.param_names_tex, 'rass');
M_.param_names_long = char(M_.param_names_long, 'rass');
M_.param_names = char(M_.param_names, 'brss');
M_.param_names_tex = char(M_.param_names_tex, 'brss');
M_.param_names_long = char(M_.param_names_long, 'brss');
M_.param_names = char(M_.param_names, 'faiss');
M_.param_names_tex = char(M_.param_names_tex, 'faiss');
M_.param_names_long = char(M_.param_names_long, 'faiss');
M_.param_names = char(M_.param_names, 'rdss');
M_.param_names_tex = char(M_.param_names_tex, 'rdss');
M_.param_names_long = char(M_.param_names_long, 'rdss');
M_.param_names = char(M_.param_names, 'crunss');
M_.param_names_tex = char(M_.param_names_tex, 'crunss');
M_.param_names_long = char(M_.param_names_long, 'crunss');
M_.param_names = char(M_.param_names, 'cpaiss');
M_.param_names_tex = char(M_.param_names_tex, 'cpaiss');
M_.param_names_long = char(M_.param_names_long, 'cpaiss');
M_.param_names = char(M_.param_names, 'XXss');
M_.param_names_tex = char(M_.param_names_tex, 'XXss');
M_.param_names_long = char(M_.param_names_long, 'XXss');
M_.param_names = char(M_.param_names, 'OMP');
M_.param_names_tex = char(M_.param_names_tex, 'OMP');
M_.param_names_long = char(M_.param_names_long, 'OMP');
M_.param_names = char(M_.param_names, 'SIG');
M_.param_names_tex = char(M_.param_names_tex, 'SIG');
M_.param_names_long = char(M_.param_names_long, 'SIG');
M_.param_names = char(M_.param_names, 'PHI');
M_.param_names_tex = char(M_.param_names_tex, 'PHI');
M_.param_names_long = char(M_.param_names_long, 'PHI');
M_.param_names = char(M_.param_names, 'BETTA');
M_.param_names_tex = char(M_.param_names_tex, 'BETTA');
M_.param_names_long = char(M_.param_names_long, 'BETTA');
M_.param_names = char(M_.param_names, 'RHOa');
M_.param_names_tex = char(M_.param_names_tex, 'RHOa');
M_.param_names_long = char(M_.param_names_long, 'RHOa');
M_.param_names = char(M_.param_names, 'ALFA');
M_.param_names_tex = char(M_.param_names_tex, 'ALFA');
M_.param_names_long = char(M_.param_names_long, 'ALFA');
M_.param_names = char(M_.param_names, 'RHOg');
M_.param_names_tex = char(M_.param_names_tex, 'RHOg');
M_.param_names_long = char(M_.param_names_long, 'RHOg');
M_.param_names = char(M_.param_names, 'GY');
M_.param_names_tex = char(M_.param_names_tex, 'GY');
M_.param_names_long = char(M_.param_names_long, 'GY');
M_.param_names = char(M_.param_names, 'OMK');
M_.param_names_tex = char(M_.param_names_tex, 'OMK');
M_.param_names_long = char(M_.param_names_long, 'OMK');
M_.param_names = char(M_.param_names, 'THETA');
M_.param_names_tex = char(M_.param_names_tex, 'THETA');
M_.param_names_long = char(M_.param_names_long, 'THETA');
M_.param_names = char(M_.param_names, 'DELTA');
M_.param_names_tex = char(M_.param_names_tex, 'DELTA');
M_.param_names_long = char(M_.param_names_long, 'DELTA');
M_.param_names = char(M_.param_names, 'EPSI');
M_.param_names_tex = char(M_.param_names_tex, 'EPSI');
M_.param_names_long = char(M_.param_names_long, 'EPSI');
M_.param_names = char(M_.param_names, 'BET');
M_.param_names_tex = char(M_.param_names_tex, 'BET');
M_.param_names_long = char(M_.param_names_long, 'BET');
M_.param_names = char(M_.param_names, 'HH');
M_.param_names_tex = char(M_.param_names_tex, 'HH');
M_.param_names_long = char(M_.param_names_long, 'HH');
M_.param_names = char(M_.param_names, 'vP');
M_.param_names_tex = char(M_.param_names_tex, 'vP');
M_.param_names_long = char(M_.param_names_long, 'vP');
M_.param_names = char(M_.param_names, 'vY');
M_.param_names_tex = char(M_.param_names_tex, 'vY');
M_.param_names_long = char(M_.param_names_long, 'vY');
M_.param_names = char(M_.param_names, 'vQ');
M_.param_names_tex = char(M_.param_names_tex, 'vQ');
M_.param_names_long = char(M_.param_names_long, 'vQ');
M_.param_names = char(M_.param_names, 'vR');
M_.param_names_tex = char(M_.param_names_tex, 'vR');
M_.param_names_long = char(M_.param_names_long, 'vR');
M_.param_names = char(M_.param_names, 'VV');
M_.param_names_tex = char(M_.param_names_tex, 'VV');
M_.param_names_long = char(M_.param_names_long, 'VV');
M_.param_names = char(M_.param_names, 'CR');
M_.param_names_tex = char(M_.param_names_tex, 'CR');
M_.param_names_long = char(M_.param_names_long, 'CR');
M_.param_names = char(M_.param_names, 'a_ss');
M_.param_names_tex = char(M_.param_names_tex, 'a\_ss');
M_.param_names_long = char(M_.param_names_long, 'a_ss');
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 27;
M_.param_nbr = 53;
M_.orig_endo_nbr = 27;
M_.aux_vars = [];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('NK_AFL15_rep_static');
erase_compiled_function('NK_AFL15_rep_dynamic');
M_.lead_lag_incidence = [
 0 6 0;
 0 7 33;
 1 8 0;
 2 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 0;
 0 14 34;
 0 15 0;
 0 16 35;
 3 17 0;
 4 18 0;
 0 19 36;
 0 20 0;
 0 21 0;
 0 22 0;
 0 23 37;
 0 24 0;
 0 25 0;
 0 26 38;
 0 27 0;
 0 28 0;
 0 29 0;
 0 30 0;
 0 31 0;
 5 32 0;]';
M_.nstatic = 16;
M_.nfwrd   = 6;
M_.npred   = 5;
M_.nboth   = 0;
M_.nsfwrd   = 6;
M_.nspred   = 5;
M_.ndynamic   = 11;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(27, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(53, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 107;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
global mcss SIG ALFA PHI CoY zss;  
M_.params( 53 ) = 0;
a_ss = M_.params( 53 );
M_.params( 34 ) = 1;
SIG = M_.params( 34 );
M_.params( 35 ) = 3;
PHI = M_.params( 35 );
M_.params( 41 ) = 2;
OMK = M_.params( 41 );
M_.params( 44 ) = 6;
EPSI = M_.params( 44 );
M_.params( 37 ) = 0.95;
RHOa = M_.params( 37 );
M_.params( 39 ) = 0.9;
RHOg = M_.params( 39 );
M_.params( 36 ) = 0.995;
BETTA = M_.params( 36 );
M_.params( 38 ) = 0.3333333333333333;
ALFA = M_.params( 38 );
M_.params( 7 ) = 1;
PAIss = M_.params( 7 );
M_.params( 43 ) = 0.025;
DELTA = M_.params( 43 );
M_.params( 40 ) = 0.2;
GY = M_.params( 40 );
M_.params( 2 ) = 0.75;
calvo = M_.params( 2 );
M_.params( 52 ) = 0.10;
CR = M_.params( 52 );
M_.params( 45 ) = 0.45;
BET = M_.params( 45 );
M_.params( 46 ) = 0.40;
HH = M_.params( 46 );
M_.params( 51 ) = M_.params(35);
VV = M_.params( 51 );
M_.params( 42 ) = 0.97;
THETA = M_.params( 42 );
M_.params( 47 ) = 1.5;
vP = M_.params( 47 );
M_.params( 48 ) = 0.125;
vY = M_.params( 48 );
M_.params( 49 ) = 0;
vQ = M_.params( 49 );
M_.params( 50 ) = 0;
vR = M_.params( 50 );
stst2=fsolve(@(stst2) (1-BETTA*(stst2+HH)/(2-BET+CR*(1+BET)))-THETA*((1-BETTA*(stst2+HH)/(2-BET+CR*(1+BET)))+(stst2+HH-(stst2+HH)/(2-BET        ))^2/(8*HH)-0.135)-0, [1]);
M_.params( 26 ) = stst2;
rass = M_.params( 26 );
M_.params( 15 ) = 1;
qss = M_.params( 15 );
M_.params( 11 ) = M_.params(26)-(1-M_.params(43));
zss = M_.params( 11 );
M_.params( 17 ) = 1-M_.params(43)+M_.params(11);
rokss = M_.params( 17 );
M_.params( 9 ) = M_.params(7)/M_.params(36);
rnss = M_.params( 9 );
M_.params( 13 ) = (M_.params(44)-1)/M_.params(44);
mcss = M_.params( 13 );
M_.params( 4 ) = M_.params(11)/(M_.params(13)*M_.params(38));
YoK = M_.params( 4 );
M_.params( 5 ) = M_.params(43);
IoK = M_.params( 5 );
IoY=IoK/YoK;
M_.params( 6 ) = 1-M_.params(40)-IoY;
CoY = M_.params( 6 );
CoK=CoY*YoK;
M_.params( 32 ) = (1-M_.params(38))*M_.params(13)/M_.params(51)*M_.params(4)/CoK;
XXss = M_.params( 32 );
M_.params( 14 ) = M_.params(32)/(1+M_.params(32));
nss = M_.params( 14 );
M_.params( 10 ) = M_.params(4)^(1/(M_.params(38)-1))*M_.params(14);
kss = M_.params( 10 );
M_.params( 23 ) = (1-M_.params(36)*(M_.params(26)+M_.params(46))/(2-M_.params(45)+M_.params(52)*(1+M_.params(45))))*M_.params(15)*M_.params(10);
bkss = M_.params( 23 );
M_.params( 25 ) = M_.params(36)*(M_.params(26)+M_.params(46))/(2-M_.params(45)+M_.params(52)*(1+M_.params(45)));
depratss = M_.params( 25 );
M_.params( 24 ) = M_.params(10)*M_.params(15)*M_.params(25);
dss = M_.params( 24 );
M_.params( 27 ) = 0.5*(1-(M_.params(26)-M_.params(25)*M_.params(9))/M_.params(46));
brss = M_.params( 27 );
M_.params( 12 ) = M_.params(10)^M_.params(38)*M_.params(14)^(1-M_.params(38))-M_.params(10)*M_.params(26)*M_.params(52)*M_.params(27);
yss = M_.params( 12 );
M_.params( 21 ) = M_.params(38)*M_.params(4);
Fkss = M_.params( 21 );
M_.params( 22 ) = (1-M_.params(38))*M_.params(12)/M_.params(14);
Fnss = M_.params( 22 );
M_.params( 16 ) = M_.params(10)*M_.params(5);
invss = M_.params( 16 );
M_.params( 18 ) = M_.params(40)*M_.params(12);
gss = M_.params( 18 );
M_.params( 8 ) = M_.params(12)*M_.params(6);
css = M_.params( 8 );
M_.params( 19 ) = M_.params(8)^(-M_.params(34));
ucss = M_.params( 19 );
M_.params( 20 ) = M_.params(35)/(M_.params(14)-1);
unss = M_.params( 20 );
M_.params( 3 ) = (M_.params(26)+M_.params(46)-M_.params(25)*M_.params(9))^2/(M_.params(46)*8);
ritcapss = M_.params( 3 );
M_.params( 28 ) = M_.params(27)*(1+M_.params(45))*.25*(1-M_.params(52))*(M_.params(9)+(M_.params(26)-M_.params(46))/M_.params(25));
faiss = M_.params( 28 );
M_.params( 29 ) = M_.params(9)*(1-M_.params(28));
rdss = M_.params( 29 );
M_.params( 33 ) = (M_.params(44)-1)*M_.params(12)*M_.params(2)/((1-M_.params(36)*M_.params(2))*(1-M_.params(2)));
OMP = M_.params( 33 );
M_.params( 1 ) = M_.params(10)*M_.params(15)*M_.params(42)*(-0.135);
PSI = M_.params( 1 );
M_.params( 30 ) = M_.params(10)*M_.params(15)*M_.params(26)*M_.params(52)*M_.params(27);
crunss = M_.params( 30 );
M_.params( 31 ) = 0;
cpaiss = M_.params( 31 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 0;
M_.Sigma_e(3, 3) = 1;
options_.steadystate.nocheck = 1;
steady;
options_.irf = 41;
options_.nograph = 1;
options_.order = 1;
var_list_=[];
var_list_ = 'y';
var_list_ = char(var_list_, 'pai');
var_list_ = char(var_list_, 'c');
var_list_ = char(var_list_, 'inv');
var_list_ = char(var_list_, 'deprat');
var_list_ = char(var_list_, 'br');
info = stoch_simul(var_list_);
save('NK_AFL15_rep_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('NK_AFL15_rep_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('NK_AFL15_rep_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('NK_AFL15_rep_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('NK_AFL15_rep_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
