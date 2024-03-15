%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program simulates Smets and Wouters (AER) from steady state for 12100 periods 	%
% Parameter values make the model bimonthly and reflect estimates based on the          %
% 	modern (89-09) sample and sector-specific (rather than GDP) deflators.		%
%	Also, re-estimated to match BLS price change frequency and no price indexation.	%
%	Removes kimball kink for final goods and for labor.				%
%											%
% BAM 08/03/11, this version 8/11/11							%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var   labobs robs pinfobs dy dc dinve dw  
	ewma epinfma  zcapf rkf kf pkf    
	cf invef yf labf wf rrf mc zcap rk k pk    
	c inve y lab pinf w r a  b g qs  ms  spinf sw kpf kp ;

varexo ea eb eg  eqs  em  epinf ew  ;

parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa
 czcap cbeta csadjcost ctou csigma chabb cfc
 cindw cprobw cindp cprobp csigl clandaw
 crpi crdy cry crr
 crhoa crhob crhog crhoqs crhoms crhopinf crhow
 ctrend
 conster cg cgamma clandap cbetabar cr cpie crk cw cikbar cik clk cky ciy ccy crkky cwhlc cwly ;

// fixed parameters
ctou=.017;
clandaw= 1.5;
cg=0.18;
curvw=0; %- Note this is set to 0 to get the Kimball Kink = 0 for labor.;

// Estimated parameters (from usmodel_bimonth_est_err.mod, for estimated curvp)

crhoa=    0.928;
crhob=    0.958;
crhog=    0.952;
crhoqs=   0.982;
crhoms= 0.692;
crhopinf= 0.491;
crhow= 0.210;
cmap = 0.394;
cmaw  = 0.315;
cgy = 1.228;
csadjcost= 5.799;
csigma=   1.269;
chabb=    0.637;
cprobw=   0.856;
csigl=    1.060;
cprobp=   0.688; //Set to match BLS frequency;
cindw=    0.466;
cindp=    0.000; //No price indexation;
czcap=    0.647;
cfc=	1.723;
crpi=     1.398;
crr=      0.884;
cry=      0.109;
crdy=     0.004;
constepinf = 0.522;
constebeta = 0.154;
constelab = -1.747; 
ctrend = 0.240;    
calfa=.234;
curvp=0; %- Note this is set to 0 to get the Kimball Kink = 0 for final goods.;

% derived parameters - shares - etc
cpie=1+constepinf/100;
cgamma=1+ctrend/100 ;
cbeta=1/(1+constebeta/100);

clandap=cfc;
cbetabar=cbeta*cgamma^(-csigma);
cr=cpie/(cbeta*cgamma^(-csigma));
crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
cikbar=(1-(1-ctou)/cgamma);
cik=(1-(1-ctou)/cgamma)*cgamma;
clk=((1-calfa)/calfa)*(crk/cw);
cky=cfc*(clk)^(calfa-1);
ciy=cik*cky;
ccy=1-cg-cik*cky;
crkky=crk*cky;
cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
cwly=1-crk*cky;

conster=(cr-1)*100;

model(linear);

// flexible economy

	      0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
	      rkf =  (wf)+labf-kf ;
	      kf =  kpf(-1)+zcapf ;
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
	      cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
	      yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
	      yf = cfc*( calfa*kf+(1-calfa)*labf +a );
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

// sticky price - wage economy

	      mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
	      zcap =  (1/(czcap/(1-czcap)))* rk ;
	      rk =  w+lab-k ;
	      k =  kp(-1)+zcap ;
	      inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          pk = -r+pinf(1)-0*b +(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
	      y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
	      y = cfc*( calfa*k+(1-calfa)*lab +a );

//	      pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1)
//               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ;

	      pinf = cbetabar*cgamma* (cprobp + (1-cprobp*cpie^(1/(clandap-1)))^(clandap)/(1-cprobp)^(clandap-1)/cpie^(1+1/(clandap-1))) * pinf(1)
		+((1-cbetabar*cgamma*cprobp)*(1-cprobp*cpie^(1/(clandap-1)))^(clandap)/cprobp/(1-cprobp)^(clandap-1)/cpie^(1/(clandap-1)))/((clandap-1)*curvp+1)*(mc)
		+ spinf ;

	      w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w)
               + 1*sw ;
	      r =  crpi*(1-crr)*pinf
               +cry*(1-crr)*(y-yf)
               +crdy*(y-yf-y(-1)+yf(-1))
               +crr*r(-1)
               +ms  ;
	      a = crhoa*a(-1)  + ea;
	      b = crhob*b(-1) + eb;
	      g = crhog*(g(-1)) + eg + cgy*ea;
	      qs = crhoqs*qs(-1) + eqs;
	      ms = crhoms*ms(-1) + em;
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=epinf;
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=ew;
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

// measurment equations

dy=y-y(-1)+ctrend;
dc=c-c(-1)+ctrend;
dinve=inve-inve(-1)+ctrend;
dw=w-w(-1)+ctrend;
pinfobs = 1*(pinf) + constepinf;
robs =    1*(r) + conster;
labobs = lab + constelab;

end;

%
% We need to initialize the variables
%

initval;

labobs = constelab;
robs = conster;
pinfobs = constepinf;
dy = ctrend;
dc = ctrend;
dinve = ctrend;
dw = ctrend;
ewma = 0;
epinfma = 0;
zcapf = 0;
rkf = 0;
kf = 0;
pkf = 0;
cf = 0;
invef = 0;
yf = 0;
labf = 0;
wf = 0;
rrf = 0;
mc = 0;
zcap = 0;
rk = 0;
k = 0;
pk = 0;
c = 0;
inve = 0;
y = 0;
lab = 0;
pinf = 0;
w = 0;
r = 0;
a = 0;
b = 0;
g = 0;
qs = 0;
ms = 0;
spinf = 0;
sw = 0;
kpf = 0;
kp = 0;

end;

%
% The standard deviations of the shocks comes from the same file as the parameters
%

%  (usmodel_bimonth_est_err.mod, estimated curvp)
shocks;
var ea;
stderr 1.172;
var eb;
stderr .031;
var eg;
stderr .586;
var eqs;
stderr .269;
var em;
stderr 0.040;
var epinf;
stderr 0;    % -- No price markup shock;
var ew;
stderr 0.483;
end;

steady;

// The execution command we run is:
stoch_simul(order = 1, noprint, IRF=0, periods = 12102, drop = 200);
%stoch_simul(order = 1, noprint, IRF=100, periods = 12102, drop = 200) y;

