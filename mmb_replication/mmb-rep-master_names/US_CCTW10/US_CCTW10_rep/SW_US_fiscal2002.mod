// Title: Shocks and Frictions in US Business Cycles: A Bayesian DSGE-Approach
// Authors: Smets, Frank and Raf Wouters
// Publication: The American Economic Review, June 2007, 97(3), 586-606.


var labobs robs pinfobs dy dc dinve dw ewma epinfma zcapf rkf kf pkf cf invef 
    yf labf wf rrf mc zcap rk k pk c inve y lab pinf w r a  b g qs  ms  spinf 
    sw kpf kp pinf4

//**************************************************************************
// Modelbase Variables                                                   
        interest inflation inflationq outputgap output pinflag1 pinflag2 pinfobs4 robs4 realinterest;            
//**************************************************************************

varexo ea eb eqs em epinf ew  

//**************************************************************************
// Modelbase Shocks                                                    
  fiscal_ dummy_MP;                                        
//**************************************************************************
  
 
parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
           czcap cbeta csadjcost ctou csigma chabb ccs cinvs cfc cindw cprobw 
           cindp cprobp csigl clandaw crdpi crhoa crhob 
           crhog crhoqs crhoms crhopinf crhow ctrend cg cgamma clandap cbetabar 
           cr cpie crk cw cikbar cik clk cky ciy ccy crkky cwhlc cwly conster

//************************************************************************** 
// Modelbase Parameters                                                  
                                                                    
        cofintintb1 cofintintb2 cofintintb3 cofintintb4              
        cofintinf0 cofintinfb1 cofintinfb2 cofintinfb3 cofintinfb4      
        cofintinff1 cofintinff2 cofintinff3 cofintinff4                  
        cofintout cofintoutb1 cofintoutb2 cofintoutb3 cofintoutb4        
        cofintoutf1 cofintoutf2 cofintoutf3 cofintoutf4;                                            
                                                                                                             
//**************************************************************************


// fixed parameters
ctou=.025; //depreciation rate
clandaw=1.5; // SS markup labor market
cg=0.18; //exogenous spending GDP-ratio
curvp=10; //curvature Kimball aggregator goods market
curvw=10; //curvature Kimball aggregator labor market

// estimated parameters initialisation
ctrend=0.4312; //quarterly trend growth rate to GDP
cgamma=ctrend/100+1;
constebeta=0.1657;
cbeta=100/(constebeta+100); //discount factor
constepinf=0.7869; //quarterly SS inflation rate
cpie=constepinf/100+1;
constelab=0.5509;

calfa=0.1901; //labor share in production

csigma=1.3808;//coefficient of relative risk aversion
cfc=1.6064; 
cgy=0.5187;

csadjcost= 5.7606; //investment adjustment cost
chabb=    0.7133;  // habit persistence 
cprobw=   0.7061;  //calvo parameter labor market
csigl=    1.8383; 
cprobp=   0.6523; //calvo parameter goods market
cindw=    0.5845; //indexation labor market
cindp=    0.2432; //indexation goods market
czcap=    0.5462;//capital utilization
crhoa=    0.9577;
crhob=    0.2194;
crhog=    0.9767;
crhoqs=   0.7113;
crhoms=0.1479; 
crhopinf=0.8895;
crhow=0.9688;
cmap = 0.7010;
cmaw  = 0.8503;

// derived from steady state
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

//Policy parameters (Taylor Rule)
cofintintb1 =  0; cofintintb2 = 0; cofintintb3 = 0; cofintintb4 = 0;
cofintinf0 = 1.5; cofintinfb1 = 0; cofintinfb2 =0; cofintinfb3 = 0; cofintinfb4 = 0; cofintinff1 = 0; cofintinff2 = 0; cofintinff3 = 0; cofintinff4 = 0;
cofintout = 0.5; cofintoutb1 = 0; cofintoutb2 = 0; cofintoutb3 = 0; cofintoutb4 = 0; cofintoutf1 = 0; cofintoutf2 = 0; cofintoutf3 = 0; cofintoutf4 = 0;


model; 

//**************************************************************************
// Definition of Modelbase Variables in Terms of Original Model Variables //*

interest   = r*4;                                                        //*
inflation  = pinf4;                                                 	 //*
inflationq = 4*pinf;                                                 	 //*
outputgap  = y-yf;                                                       //*
output     = y;                                                          //*
//**************************************************************************


//**************************************************************************                                                                    
// Policy Rule                                                           //*

interest =   dummy_MP*0+(1-dummy_MP)*(cofintintb1*interest(-1)                                    //* 
           + cofintinf0*inflation                                       //* 
           + cofintinfb1*inflation(-1)                                  //* 
           + cofintout*outputgap 	                                     //* 
           + cofintoutb1*outputgap(-1));                                   //*          
//**************************************************************************

// Original Model Code:

// flexible economy

	      0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
	      rkf =  (wf)+labf-kf ;
	      kf =  kpf(-1)+zcapf ;
          //investment Euler equation
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
          //consumption Euler equation
	      cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
          //aggregate resource constraint
	      yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
          // aggregate production function
	      yf = cfc*( calfa*kf+(1-calfa)*labf +a );
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
        //accumulation of installed capital
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

// sticky price - wage economy
          //marginal cost
	      mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
          // capital utilization 
	      zcap =  (1/(czcap/(1-czcap)))* rk ;
          //rental rate of capital
	      rk =  w+lab-k ;
          // Capital installed used one period later in production 
	      k =  kp(-1)+zcap ;
          //investment Euler equation
	      inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          //arbitrage equation for value of capital
          pk = -r+pinf(1)-0*b +(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
          //consumption Euler equation
	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
          //aggregate resource constraint
	      y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
          // aggregate production function
	      y = cfc*( calfa*k+(1-calfa)*lab +a );
          //Phillips Curve
	      pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 
	      w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w) 
               + 1*sw ;
	      a = crhoa*a(-1)  + ea;
	      b = crhob*b(-1) + eb;
          // exogenous spending (also including net exports)
	      g = fiscal_;
	      qs = crhoqs*qs(-1) + eqs;
	      ms = crhoms*ms(-1) + em;
          //cost push shock
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=epinf;
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=ew; 
          //accumulation of installed capital
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;
        dy=y-y(-1)+ctrend;
        dc=c-c(-1)+ctrend;
        dinve=inve-inve(-1)+ctrend;
        dw=w-w(-1)+ctrend;
        pinfobs = 1*(pinf) + constepinf;
        pinfobs4=pinfobs*4;
        pinf4= pinf + pinf(-1) + pinflag2 + pinflag2(-1);
        pinflag1=pinf(-1);
        pinflag2=pinflag1(-1);
        robs =    1*(r) + conster;
        robs4=robs*4;
        labobs = lab + constelab;
        //Real interest rate
        realinterest=interest-inflationq;

end; 

shocks;
var fiscal_;
// fiscal package, numbers 20/02/09  
// federal and state 60% one quarter slip
periods 1:1   2:3       4:4  5:7  8:8  9:11  12:12  13:15  16:16  17:20 21:23  24:27  28:31    32:35   36:39    40:43;
values  0.1964  0.4657  0.59 0.77 0.73  0.49  0.40   0.25   0.17   0.16  0.15   0.06  ,-0.01  ,-0.01  ,-0.01,   0.01;

var dummy_MP;
periods 1:4;
values 1;
end;
