%Fiscal policy parameters - these are taken from a estimated SVAR model -
%this is a quick solution - these are going to be imported from another
%.m-file asap! 
%The order of these in SVAR: tauy tauc g_gap
%the order in SOE          : tauk tauy tauc tauw g_gap
A0_fis = [ 1   0   0;  ...
           0   1   0;  ... 
           0   0   1];
B0_fis = [ 0.132418        0           0;  ...
           0               0.142992    0;  ... 
           0               0           0.414848];     
     
%NOTE: the lag polynomial parameters are reported in reduced form in Eviews, so no need to multiply them...        
A1_fis = [ 1.574508    0.063984   -0.039498;...            
           0.036577    0.762311    0.028616;...            
          -0.048547    0.053786    0.562102];
        
A2_fis = [-0.646174    0.037071   -0.017812;...            
          -0.029465    0.211986   -0.052633;...            
           0.105095   -0.139486    0.199290];            
         
%Forming reduced form covariance matrix    
B0_fis = inv(A0_fis)*B0_fis;

A1_fis = [rho_tauk 0               0         0       0            ;...
          0       A1_fis(1,1)  A1_fis(1,2)  0       A1_fis(1,3)  ;...
          0       A1_fis(2,1)  A1_fis(2,2)  0       A1_fis(2,3)  ;...
          0       0               0         rho_tauw 0            ;...
          0       A1_fis(3,1)  A1_fis(3,2)  0       A1_fis(3,3)  ];
        
A2_fis = [0       0            0            0       0            ;...
          0       A2_fis(1,1)  A2_fis(1,2)  0       A2_fis(1,3)  ;...
          0       A2_fis(2,1)  A2_fis(2,2)  0       A2_fis(2,3)  ;...
          0       0            0            0       0            ;...
          0       A2_fis(3,1)  A2_fis(3,2)  0       A2_fis(3,3)  ];
 
B0_fis = [epstauk 0            0            0       0            ;...
          0       B0_fis(1,1)  B0_fis(1,2)  0       B0_fis(1,3)  ;...
          0       B0_fis(2,1)  B0_fis(2,2)  0       B0_fis(2,3)  ;...
          0       0            0            epstauw 0            ;...
          0       B0_fis(3,1)  B0_fis(3,2)  0       B0_fis(3,3)  ];
     
%Foreign economy parameters
A0_for = [ 1.000000    0.000000    0.000000;...            
           0.000000    1.000000    0.000000;...            
           0.099868   -0.168144    1.000000];            


B0_for = [ 0.258501    0.000000    0.000000;...             
           0.000000    0.337480    0.000000;...             
           0.000000    0.000000    0.225671];             
                
A1_for = [ 0.176229    0.279843   -0.208028;...             
          -0.062804    1.048208    0.284264;...             
           0.321600    0.132657    0.919272];             
                   
A2_for = [ 0.264370   -0.156559    0.081247;...   
           0.153104   -0.130402   -0.578732;...             
           0.143343   -0.200941   -0.195221];             
                   
A3_for = [ 0.358584    0.041918    0.239725;...             
          -0.165318    0.100930    0.475599;...             
          -0.096743    0.168364    0.215159];             
                    
A4_for = [ 0.027420   -0.056139   -0.087750;...             
           0.135853   -0.214142   -0.281115;...             
          -0.079107   -0.160620   -0.115767];
          
%Forming the reduced form covariance matrix    
B0_for   = inv(A0_for)*B0_for;