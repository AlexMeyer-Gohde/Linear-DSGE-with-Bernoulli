%Description: ....
%....
% set_param_value('crpi',10)
% [dr_new, info, M, options, oo_new] = resol(0, M_, options_, oo_);
% ALPHA_ZS_dynare=[zeros(n,nstatic) oo_new.dr.ghx zeros(n,nfwrd)];
% ALPHA_ZS_dynare=ALPHA_ZS_dynare(oo_.dr.inv_order_var,oo_.dr.inv_order_var);
% ALPHA_ZS_dynare(2,29)

clear 

STEPS=-1:0.25:8;
STEPS=[STEPS 9:0.25:13];

for KK=1:length(STEPS)

run_time_reps=1;
step_size=STEPS(KK);
%newton_options.dynare_reduced_sylvester=1;
newton_options.sylvester_method="dlyap_stripped";
newton_options.maximum_iterations=1000;
addpath('C:\dynare\5.1\matlab')
addpath('..\algorithm\')
%YourPath=''
%YourPath='C:\Users\saecker.ITS\PowerFolders\Newton_Saecker (Alexander Meyer-Gohde)\code\workspace\2021_10_15_js';



%run dynare
%dynare ([mmb_vec{k} '_rep']) 
dynare US_SW07_rep_SMG noclearall nograph nostrict

CRPI_vector=linspace(1.5,1.5*(1+10^(-step_size)),10);
CRY_vector=linspace(0.125,0.125*(1+10^(-step_size)),10);
AMG_Results=zeros(16,7,length(CRPI_vector),length(CRY_vector));

times=zeros(1,run_time_reps);

for JJ=1:length(CRPI_vector)
    for II=1:length(CRY_vector)
        set_param_value('crpi',CRPI_vector(JJ))
        set_param_value('cry',CRY_vector(II))
if JJ==1 && II==1
    X_0_1=oo_.dr.ghx(M_.nstatic+1:end,:);%zeros( M_.ndynamic,M_.nspred);
    X_0_2=X_0_1;
    X_0_3=X_0_1;
    X_0_4=X_0_1;
    X_0_5=X_0_1;
    X_0_6=X_0_1;
    X_0_7=X_0_1;
    X_0_8=X_0_1;
    X_0_9=X_0_1;
    X_0_10=X_0_1;
    X_0_11=X_0_1;
    X_0_12=X_0_1;
    X_0_13=X_0_1;
    X_0_14=X_0_1;
elseif JJ~=1 && II==1
    X_0_1=X_0_1_begin;
    X_0_2=X_0_2_begin;
    X_0_3=X_0_3_begin;
    X_0_4=X_0_4_begin;
    X_0_5=X_0_5_begin;
    X_0_6=X_0_6_begin;
    X_0_7=X_0_7_begin;
    X_0_8=X_0_8_begin;
    X_0_9=X_0_9_begin;
    X_0_10=X_0_10_begin;
    X_0_11=X_0_11_begin;
    X_0_12=X_0_12_begin;
    X_0_13=X_0_13_begin;
    X_0_14=X_0_14_begin;
end

[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, M_.endo_names);%steady;[dr_new, info, M, options, oo_new] = resol(0, M_, options_, oo_);

%%%%% current problem: dynare_to_matrix_quadratic needs to be located in
%%%%% mmb-rep-folders (e.g. BRA_SAMBA08_rep)
AMG_Results(1,:,JJ,II)=[M_.nstatic, M_.nfwrd, M_.npred, M_.nboth, M_.nsfwrd, M_.nspred, M_.ndynamic];

[matrix_quadratic, jacobia_]=create_reduced_matrix_quadratic_from_dynare(M_,oo_);

%tic; [info, oo_, options_]  = stoch_simul(M_, options_, oo_, var_list_); toc    
 for jj=1:run_time_reps;tic;[dr,info] = dyn_first_order_solver(jacobia_,M_,oo_.dr,options_,0); total_time(jj)=toc;end;   AMG_Results(2,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));   

 ALPHA_ZS_dynare=[zeros(M_.endo_nbr,M_.nstatic) oo_.dr.ghx zeros(M_.endo_nbr,M_.nfwrd)];
X_dynare=ALPHA_ZS_dynare;
matrix_quadratic.X=ALPHA_ZS_dynare;
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(2,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(2,3:5,JJ,II)=NaN(1,3);
end
 

try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(2,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(2,3:5,JJ,II)=NaN(1,3);
end

try
 newton_options.algorithm='baseline';
 newton_options.M_=M_;
 newton_options.convergence_metric="reldiff";
 newton_options.initial=X_0_1;
 total_time=[]; for jj=1:run_time_reps;tic;[X,X_additional] = newton_matrix_quadratic(matrix_quadratic,newton_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);total_time(jj)=toc;end;   
AMG_Results(3,1,JJ,II) =mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(3,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(3,end,JJ,II)=X_additional;
X_0_1=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(3,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(3,3:5,JJ,II)=NaN(1,3);
end
if X_additional==newton_options.maximum_iterations; 
AMG_Results(3,:,JJ,II)=NaN(1,7); 
end
catch
AMG_Results(3,:,JJ,II)=NaN(1,7);
end

clear bernoulli_options
bernoulli_options=[1 0 0 0 0 0 0 1 0 4 100 matrix_quadratic.ndynamic*eps 0];
X_0_bernoulli=X_0_1;
[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_1,bernoulli_options); 

matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);max(max(abs(X_dynare-matrix_quadratic.X)))
maximum_iterations_total=bernoulli_options(11)+1;
while X_additional==bernoulli_options(11)+1%&&max(max(abs(X_dynare-matrix_quadratic.X)))>1e-5
    if maximum_iterations_total<50000
    bernoulli_options(11)=2*bernoulli_options(11);
        [X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options);
        matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);max(max(abs(X_dynare-matrix_quadratic.X)))
        %X_additional=X_additional+X_additional_plus;
        maximum_iterations_total=maximum_iterations_total+bernoulli_options(11)
        X_0_bernoulli=X(:,1:M_.nspred);
    else
        matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);
        if max(max(abs(X_dynare-matrix_quadratic.X)))<1e-5
            break
        else
            bernoulli_options(11)=1;
            break
        end
    end
end
if bernoulli_options(11)~=1
   bernoulli_options(11)=maximum_iterations_total;
end
run_time_reps_old=run_time_reps;
if bernoulli_options(11)>=500
    run_time_reps=run_time_reps/5;
end
if bernoulli_options(11)>=1000
        run_time_reps=run_time_reps/2;
end
bernoulli_options(11)=10*bernoulli_options(11);

run_time_reps=ceil(run_time_reps);



try
 %   bernoulli_options.baseline=1;
    X_0_bernoulli=X_0_2;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(4,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(4,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(4,end,JJ,II)=X_additional;
X_0_2=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(4,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(4,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(4,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(4,:,JJ,II)=NaN(1,7);
end


%bernoulli_options.maximum_iterations=1;


try
 
%bernoulli_options.baseline=0;bernoulli_options.mbi=1;
bernoulli_options(1)=0;bernoulli_options(2)=1;
    X_0_bernoulli=X_0_3;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(5,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8))); 
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(5,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(5,end,JJ,II)=X_additional;
X_0_3=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(5,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(5,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(5,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(5,:,JJ,II)=NaN(1,7);
end

try
    %bernoulli_options.maximum_iterations=maximum_iterations_total;
%bernoulli_options.mbi=0; bernoulli_options.line_search=1;

bernoulli_options(2)=0;bernoulli_options(9)=1;
    X_0_bernoulli=X_0_4;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(6,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(6,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(6,end,JJ,II)=X_additional;
X_0_4=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(6,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(6,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iter12,ations; 
% AMG_Results(6,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(6,:,JJ,II)=NaN(1,7);
end

try
%bernoulli_options.line_search=0; bernoulli_options.newton='block';
bernoulli_options(9)=0;bernoulli_options(3)=1;bernoulli_options(4)=1;
    X_0_bernoulli=X_0_5;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(7,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(7,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(7,end,JJ,II)=X_additional;
X_0_5=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(7,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(7,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(7,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(7,:,JJ,II)=NaN(1,7);
end


try
 
%bernoulli_options.line_search=0; bernoulli_options.newton='column';
bernoulli_options(4)=0;bernoulli_options(3)=1;bernoulli_options(5)=1;
    X_0_bernoulli=X_0_6;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(8,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(8,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(8,end,JJ,II)=X_additional;
X_0_6=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(8,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(8,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(8,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(8,:,JJ,II)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=0; bernoulli_options.newton='row';
bernoulli_options(5)=0;bernoulli_options(3)=1;bernoulli_options(6)=1;
    X_0_bernoulli=X_0_7;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(9,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(9,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(9,end,JJ,II)=X_additional;
X_0_7=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(9,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(9,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(9,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(9,:,JJ,II)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='block';
bernoulli_options(6)=0;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(4)=1;
    X_0_bernoulli=X_0_8;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(10,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(10,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(10,end,JJ,II)=X_additional;
X_0_8=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(10,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(10,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(10,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(10,:,JJ,II)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='column';
bernoulli_options(4)=0;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(5)=1;
    X_0_bernoulli=X_0_9;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(11,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(11,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(11,end,JJ,II)=X_additional;
X_0_9=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(11,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(11,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(11,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(11,:,JJ,II)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='row';
bernoulli_options(5)=0;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(6)=1;
    X_0_bernoulli=X_0_10;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(12,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(12,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(12,end,JJ,II)=X_additional;
X_0_10=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(12,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(12,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(12,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(12,:,JJ,II)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=0; bernoulli_options.newton='block'; bernoulli_options.newton_power=1/3;
bernoulli_options(6)=0;bernoulli_options(3)=1;bernoulli_options(9)=0;bernoulli_options(8)=1/3;bernoulli_options(4)=1;
    X_0_bernoulli=X_0_11;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(13,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(13,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(13,end,JJ,II)=X_additional;
X_0_11=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(13,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(13,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(13,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(13,:,JJ,II)=NaN(1,7);
end
%bernoulli_options = rmfield(bernoulli_options,'newton_power');

try
 
%bernoulli_options.line_search=0; bernoulli_options.newton='opt';
bernoulli_options(7)=1;bernoulli_options(3)=1;bernoulli_options(9)=0;bernoulli_options(8)=1;bernoulli_options(4)=0;
    X_0_bernoulli=X_0_12;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(14,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(14,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(14,end,JJ,II)=X_additional;
X_0_12=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(14,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(14,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(14,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(14,:,JJ,II)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='block'; bernoulli_options.newton_power=1/3;
bernoulli_options(7)=0;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(8)=1/3;bernoulli_options(4)=1;
    X_0_bernoulli=X_0_13;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(15,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(15,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(15,end,JJ,II)=X_additional;
X_0_13=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(15,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(15,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(15,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(15,:,JJ,II)=NaN(1,7);
end
%bernoulli_options = rmfield(bernoulli_options,'newton_power');
try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='opt';
bernoulli_options(7)=1;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(8)=1;bernoulli_options(4)=0;
    X_0_bernoulli=X_0_14;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(16,1,JJ,II) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(16,2,JJ,II)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(16,end,JJ,II)=X_additional;
X_0_14=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(16,3:5,JJ,II)=errors([4,7,8],1)';
catch
AMG_Results(16,3:5,JJ,II)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(16,:,JJ,II)=NaN(1,7); 
% end
catch
AMG_Results(16,:,JJ,II)=NaN(1,7);
end

if  II==1
    X_0_1_begin=X_0_1;
     X_0_2_begin=X_0_2;
     X_0_3_begin=X_0_3;
    X_0_4_begin=X_0_4;
    X_0_5_begin=X_0_5;
    X_0_6_begin=X_0_6;
     X_0_7_begin=X_0_7;
    X_0_8_begin=X_0_8;
    X_0_9_begin=X_0_9;
    X_0_10_begin=X_0_10;
    X_0_11_begin=X_0_11;
    X_0_12_begin=X_0_12;
    X_0_13_begin=X_0_13;
    X_0_14_begin=X_0_14;
end





%run different Newton methods
%newton_solvent_2(matrix_quadratic);
%newton_solvent_3(matrix_quadratic);


sprintf('Iteration %d and %d of %d and %d',JJ,II,length(CRPI_vector),length(CRY_vector))
%save certain results somewhere
%clearvars -except JJ II loop_n AMG_Results YourPath mmb_vec run_time_reps max_newton_it 


    end
end
results=median(AMG_Results,4);results=median(results,3)

total_AMG_Results{KK}=AMG_Results;
total_results{KK}=results;
end
%cd([YourPath])
%save Policy_Run_AMG_JS_new
combined_results=[];
for j=1:length(STEPS)
combined_results(:,:,j)=total_results{j};
end
save('Combined_Policy_Run_fill_2_errors_2.mat')

