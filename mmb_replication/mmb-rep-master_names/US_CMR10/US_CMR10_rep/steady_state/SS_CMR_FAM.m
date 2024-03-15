%---------------------------------------------------------------------------------------
%      THE STEADY STATE FOR THE CMR Financial Accelerator Model
%      Code based on SS_CMR.m program for the baseline model by Fabio Verona
%      July 2012 
%---------------------------------------------------------------------------------------

clear all
clc

%Standard Parameters (Table 1 of the paper)
%% Household sector
betaUU        =   0.9966;
sigmaLXUU     =   1; 
bUU           =   0.63; 
lambdawUU     =   1.05;

%% Goods producing sector 
muzstarUU     =   1.0036; %as in Fabio Verona code
deltaUU       =   0.025; 
alphaUU       =   0.4; 
lambdafUU     =   1.2; 
phiUU         =   0.07;

%% Entrepreneurs
gammaUU       =   0.9762; 
muUU          =   0.94;
sigmaUU       =   sqrt (0.24); 
weUU          =   0.009;       % as in Fabio Verona code, UNKNOWN SOURCE  
bigthetaUU    =   0.1;         % as in Fabio Verona code, UNKNOWN SOURCE  

%% Policy 
taucUU        =   0.05;
taukUU        =   0.32;
taulUU        =   0.24;

%% PARAMETERS ESTIMATED - see table A.2
xipUU         =  0.702; % Calvo SP  
xiwUU         =  0.771; % Calvo SW  
iota1UU       =  0.159; % Weight on steady state inflation (CSP)0.362
iotaw1UU      =  0.285; % Weight on steady state inflation (CSW)0.641
iotamuUU      =  0.917;  % Weight on technology growth0.93
SdouprXUU     =  29.31; % Investment adjustament cost 26.64;
sigmaaUU      =  18.85;% Capacity utilization 19.718
aptilUU       =  1.817; % Weight on inflation in Taylor rule 1.849
aytilUU       =  0.310; % Weight on output growth in Taylor rule0.311;
apdeltilUU       =  0.181; % Weight on inflation change in Taylor rule 
rhotilUU      =  0.877;  % Coefficient on lagged interest rate 0.88


%% OTHER PARAMETERS 
etagUU        =   0.2; 
muupUU        =   1;
epsilUU       =   1;
zetacUU       =   1;
zetaiUU       =   1;
tauoUU        =   1;
epsilUU       =   1;
upsilUU       =   1.0035;       %   UNKNOWN SOURCE

%% ASSUMPTION to compute the steady state (see CMR JMCB 2003 paper);
rkUU          =   0.059; 
vlUU          =   1;



%% -------------- COMPUTE THE SS for some endogenous variables -------------------------------
qUU           =  1/muupUU;
piUU          =  1; % assumption
pitargetUU    =  piUU;
pibarUU       =  piUU;
ReXUU         =  (piUU*muzstarUU)/betaUU - 1;
uUU           =  1;
pstarUU       =  1;
wstarUU       =  1;
sUU           =  1/lambdafUU;
RkXUU         = ((1-taukUU)*rkUU+1-deltaUU)*piUU/upsilUU + taukUU*deltaUU - 1; 


Annual_inflation_rate=(piUU-1)*400;
Annualized_Re_SS=((ReXUU+1)^4-1)*100;
Annualized_Rk_SS=((RkXUU+1)^4-1)*100;
%% -------------------------- COMPUTE THE SS for the "contract" -------------------------------
x0=0.12;

p = fsolve(@(x) funcontract2CMR(x,RkXUU,sigmaUU,ReXUU,muUU), x0, ...
                            optimset('TolFun',1e-10,'MaxIter',100000000,'MaxFunEvals',10000000));
omegabarUU=p(1);


%omegabarUU = 0.225;
x0=[5;4];

x = fsolve(@(x) funcontractCMR(x,qUU,RkXUU,sigmaUU,ReXUU,gammaUU,piUU,muUU,weUU,muzstarUU,omegabarUU), x0, ...
                            optimset('TolFun',1e-10,'MaxIter',100000000,'MaxFunEvals',10000000));
kbarUU=x(1);
nUU=x(2);

if nUU>kbarUU
    disp('warning: nUU>kbarUU')
end

if min(x)<0
    disp('**************************************************');
    disp('this is not the true steady state for the contract');
    disp('**************************************************');
end

iUU = (1-(1-deltaUU)/upsilUU/muzstarUU)*kbarUU;

k1 = ((1/(1-alphaUU))^(1-alphaUU))*((rkUU/alphaUU)^alphaUU);
k2 = rkUU/alphaUU/((upsilUU*muzstarUU)^(1-alphaUU));
wUU = (sUU/k1)^(1/(1-alphaUU));
hUU = kbarUU*((k2/sUU)^(1/(1-alphaUU)));
cUU = (1-etagUU)*(vlUU*((kbarUU/(muzstarUU*upsilUU))^alphaUU)*(hUU^(1-alphaUU)) - phiUU) - bigthetaUU*(1-gammaUU)*(nUU-weUU)/gammaUU -...
              ( ( ( normcdf(( (log(omegabarUU)+sigmaUU^2/2)/sigmaUU )-sigmaUU)+omegabarUU*(1-normcdf(( (log(omegabarUU)+sigmaUU^2/2)/sigmaUU ))) )-...
              ( (1-muUU)*normcdf(( (log(omegabarUU)+sigmaUU^2/2)/sigmaUU )-sigmaUU)+...
              omegabarUU*(1-normcdf(( (log(omegabarUU)+sigmaUU^2/2)/sigmaUU ))) ) )*(1+RkXUU)*qUU*kbarUU/(muzstarUU*piUU) ) -iUU ;
     
uzcUU = 1/cUU*((muzstarUU-betaUU*bUU)/(muzstarUU-bUU));
lambdazUU = uzcUU/(1+taucUU); 
YUU=( (pstarUU^(lambdafUU/(lambdafUU-1))) * ( epsilUU * vlUU * ( (uUU*kbarUU/(muzstarUU*upsilUU))^alphaUU ) * ( ( hUU*(wstarUU ^ (lambdawUU/(lambdawUU-1)) ) )^(1-alphaUU) ) - phiUU ) );
gUU = etagUU*YUU;
BUU = (kbarUU*qUU-nUU);
PrUU=( ( normcdf(( (log(omegabarUU)+sigmaUU^2/2)/sigmaUU )-sigmaUU)+omegabarUU*(1-normcdf(( (log(omegabarUU)+sigmaUU^2/2)/sigmaUU ))) )-( (1-muUU)*normcdf(( (log(omegabarUU)+sigmaUU^2/2)/sigmaUU )-sigmaUU)+omegabarUU*(1-normcdf(( (log(omegabarUU)+sigmaUU^2/2)/sigmaUU ))) ) )*(1+RkXUU)*kbarUU*qUU/( kbarUU*qUU-nUU);
if lambdazUU<0
    disp('warning: lambdazUU<0')
end

%% other variables 

FpXUU = lambdazUU*YUU/(1-betaUU*xipUU);
FwXUU = hUU*lambdazUU*(1-taulUU)/(lambdawUU*(1-betaUU*xiwUU));
psiLXUU = FwXUU*wUU*(1-betaUU*xiwUU)/(hUU^(1+sigmaLXUU));

%% parameters for shocks - see Table A2 in CMR ECB WP
rho1_lambdafU = 0.627;
rho1_muupU    = 0.983;
rho1_gU       = 0.930;
rho1_muzstarU = 0.166;
rho1_gammaU   = 0.598;
rho1_epsilU   = 0.816;
rho1_sigmaU   = 0.722;
rho1_zetaiU   = 0.434;
rho1_tauoU    = 0.929;

std1_lambdafU = 0.022;
std1_muupU    = 0.003;
std1_gU       = 0.021;
std1_muzstarU = 0.007;
std1_gammaU   = 0.014;
std1_epsilU   = 0.008;
std1_sigmaU   = 0.076;
std1_zetaiU   = 0.019;
std1_tauoU    = 0.135;
std1_xpU      = 0.519; %monetary policy shock

% additional parameters (as in Fabio Verona)
rho1_pitargetU = 0.975; 
std1_pitargetU = 0.0007785;
zetaUU       =   1;
taudUU = 0;
actilUU = 0;



%% save the steady state values

save ss_parameters_US_FAM lambdawUU sigmaLXUU betaUU xiwUU bUU psiLXUU iotaw1UU iotamuUU muUU weUU bigthetaUU sigmaaUU SdouprXUU...
    zetaUU upsilUU alphaUU deltaUU xipUU iota1UU etagUU taukUU taulUU taucUU taudUU aptilUU aytilUU actilUU rhotilUU phiUU pibarUU vlUU...
    piUU sUU rkUU iUU uUU omegabarUU RkXUU nUU qUU  uzcUU lambdazUU cUU wUU hUU kbarUU ReXUU ...
    pstarUU wstarUU FpXUU FwXUU YUU PrUU BUU lambdafUU pitargetUU muupUU gUU muzstarUU gammaUU epsilUU sigmaUU zetacUU zetaiUU tauoUU...
    std1_xpU std1_lambdafU std1_muupU std1_gU std1_muzstarU std1_gammaU std1_epsilU std1_sigmaU std1_zetaiU std1_tauoU ...
    rho1_lambdafU rho1_muupU rho1_gU rho1_muzstarU rho1_gammaU rho1_epsilU rho1_sigmaU rho1_zetaiU rho1_tauoU rho1_pitargetU std1_pitargetU;

