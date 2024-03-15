%Description: ....
%....
%Alexander Meyer-Gohde


%clear all
 %load First_Run_AMG_test
 %loop_k_start=loop_k;
loop_k_start=1;%7;%93;%92;%loop_k_start=93


run_time_reps=10;
%newton_options.dynare_reduced_sylvester=1;
newton_options.sylvester_method="dlyap_stripped";
%newton_options.sylvester_method="slicot";
newton_options.maximum_iterations=1000;
addpath('C:\dynare\5.1\matlab')
addpath('..\algorithm\')
addpath('..\mmb_replication\')
%YourPath=''
%YourPath='C:\Users\saecker.ITS\PowerFolders\Newton_Saecker (Alexander Meyer-Gohde)\code\workspace\2021_10_15_js';
YourPath=pwd;
cd (YourPath)
addpath(pwd)

%fileID = fopen('mmb_test.txt','r');
fileID = fopen('mmb_names.txt','r');
mmbline = fgetl(fileID);        
%geht hier vielleicht auch besser (ff. 5 Zeilen von mathworks-website)
mmb_vec = cell(0,1);            
while ischar(mmbline)           
    mmb_vec{end+1,1} = mmbline; 
    mmbline = fgetl(fileID);    
end    
loop_n=size(mmb_vec,1);

if ~exist('AMG_Results')
AMG_Results=NaN(12,7,loop_n);
end

for loop_k=loop_k_start:loop_n
    %k=1;    %for testing
    [loop_k loop_k/loop_n]
%change directory to folder path in MMB

cd([YourPath '..\..\mmb_replication\mmb-rep-master_names\' mmb_vec{loop_k} '\' mmb_vec{loop_k} '_rep'])

%run dynare
%dynare ([mmb_vec{k} '_rep']) 
eval(['dynare ', mmb_vec{loop_k}, '_rep noclearall nograph nostrict'])
%%%%% current problem: dynare_to_matrix_quadratic needs to be located in
%%%%% mmb-rep-folders (e.g. BRA_SAMBA08_rep)
AMG_Results(1,:,loop_k)=[M_.nstatic, M_.nfwrd, M_.npred, M_.nboth, M_.nsfwrd, M_.nspred, M_.ndynamic];


[matrix_quadratic, jacobia_]=create_reduced_matrix_quadratic_from_dynare(M_,oo_);

%tic; [info, oo_, options_]  = stoch_simul(M_, options_, oo_, var_list_); toc    
total_time=[]; for jj=1:run_time_reps;tic;[dr,info] = dyn_first_order_solver(jacobia_,M_,oo_.dr,options_,0); total_time(jj)=toc;end;   AMG_Results(2,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));   

 
 
ALPHA_ZS_dynare=[zeros(M_.endo_nbr,M_.nstatic) oo_.dr.ghx zeros(M_.endo_nbr,M_.nfwrd)];
X_dynare=ALPHA_ZS_dynare;
matrix_quadratic.X=ALPHA_ZS_dynare;
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(2,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(2,3:5,loop_k)=NaN(1,3);
end



try
 newton_options.algorithm='baseline';
 newton_options.M_=M_;
 newton_options.convergence_metric="reldiff";
 total_time=[]; for jj=1:run_time_reps;tic;[X,X_additional] = newton_matrix_quadratic(matrix_quadratic,newton_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);total_time(jj)=toc;end;   
AMG_Results(3,1,loop_k) =mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(3,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(3,end,loop_k)=X_additional;
X_0_1=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(3,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(3,3:5,loop_k)=NaN(1,3);
end
if X_additional==newton_options.maximum_iterations; 
AMG_Results(3,:,loop_k)=NaN(1,7); 
end
catch
AMG_Results(3,:,loop_k)=NaN(1,7);
end

clear bernoulli_options eps
bernoulli_options=[1 0 0 0 0 0 0 1 0 4 100 matrix_quadratic.ndynamic*eps 0];
X_0_bernoulli=zeros(M_.ndynamic,M_.nspred);%X_0_1;
[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); 

matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);max(max(abs(X_dynare-matrix_quadratic.X)))
maximum_iterations_total=bernoulli_options(11)+1;
while X_additional==bernoulli_options(11)+1%&&max(max(abs(X_dynare-matrix_quadratic.X)))>1e-5
    if maximum_iterations_total<500000
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
bernoulli_options(11)=2*bernoulli_options(11);




try
 %   bernoulli_options.baseline=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_2;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(4,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(4,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(4,end,loop_k)=X_additional;
X_0_2=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(4,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(4,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(4,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(4,:,loop_k)=NaN(1,7);
end


%bernoulli_options.maximum_iterations=1;


try
 
%bernoulli_options.baseline=0;bernoulli_options.mbi=1;
bernoulli_options(1)=0;bernoulli_options(2)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_3;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(5,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8))); 
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(5,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(5,end,loop_k)=X_additional;
X_0_3=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(5,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(5,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(5,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(5,:,loop_k)=NaN(1,7);
end

try
    %bernoulli_options.maximum_iterations=maximum_iterations_total;
%bernoulli_options.mbi=0; bernoulli_options.line_search=1;

bernoulli_options(2)=0;bernoulli_options(9)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_4;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(6,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(6,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(6,end,loop_k)=X_additional;
X_0_4=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(6,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(6,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iter12,ations; 
% AMG_Results(6,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(6,:,loop_k)=NaN(1,7);
end

try
%bernoulli_options.line_search=0; bernoulli_options.newton='block';
bernoulli_options(9)=0;bernoulli_options(3)=1;bernoulli_options(4)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_5;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(7,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(7,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(7,end,loop_k)=X_additional;
X_0_5=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(7,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(7,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(7,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(7,:,loop_k)=NaN(1,7);
end


try
 
%bernoulli_options.line_search=0; bernoulli_options.newton='column';
bernoulli_options(4)=0;bernoulli_options(3)=1;bernoulli_options(5)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_6;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(8,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(8,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(8,end,loop_k)=X_additional;
X_0_6=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(8,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(8,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(8,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(8,:,loop_k)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=0; bernoulli_options.newton='row';
bernoulli_options(5)=0;bernoulli_options(3)=1;bernoulli_options(6)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_7;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(9,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(9,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(9,end,loop_k)=X_additional;
X_0_7=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(9,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(9,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(9,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(9,:,loop_k)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='block';
bernoulli_options(6)=0;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(4)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_8;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(10,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(10,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(10,end,loop_k)=X_additional;
X_0_8=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(10,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(10,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(10,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(10,:,loop_k)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='column';
bernoulli_options(4)=0;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(5)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_9;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(11,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(11,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(11,end,loop_k)=X_additional;
X_0_9=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(11,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(11,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(11,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(11,:,loop_k)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='row';
bernoulli_options(5)=0;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(6)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_10;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(12,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(12,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(12,end,loop_k)=X_additional;
X_0_10=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(12,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(12,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(12,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(12,:,loop_k)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=0; bernoulli_options.newton='block'; bernoulli_options.newton_power=1/3;
bernoulli_options(6)=0;bernoulli_options(3)=1;bernoulli_options(9)=0;bernoulli_options(8)=1/3;bernoulli_options(4)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_11;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(13,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(13,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(13,end,loop_k)=X_additional;
X_0_11=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(13,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(13,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(13,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(13,:,loop_k)=NaN(1,7);
end
%bernoulli_options = rmfield(bernoulli_options,'newton_power');

try
 
%bernoulli_options.line_search=0; bernoulli_options.newton='opt';
bernoulli_options(7)=1;bernoulli_options(3)=1;bernoulli_options(9)=0;bernoulli_options(8)=1;bernoulli_options(4)=0;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_12;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(14,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(14,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(14,end,loop_k)=X_additional;
X_0_12=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(14,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(14,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(14,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(14,:,loop_k)=NaN(1,7);
end

try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='block'; bernoulli_options.newton_power=1/3;
bernoulli_options(7)=0;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(8)=1/3;bernoulli_options(4)=1;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_13;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(15,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(15,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(15,end,loop_k)=X_additional;
X_0_13=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(15,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(15,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(15,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(15,:,loop_k)=NaN(1,7);
end
%bernoulli_options = rmfield(bernoulli_options,'newton_power');
try
 
%bernoulli_options.line_search=1; bernoulli_options.newton='opt';
bernoulli_options(7)=1;bernoulli_options(3)=1;bernoulli_options(9)=1;bernoulli_options(8)=1;bernoulli_options(4)=0;
    X_0_bernoulli=zeros( M_.ndynamic,M_.nspred);%X_0_14;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_fast(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
AMG_Results(16,1,loop_k) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; AMG_Results(16,2,loop_k)=max(max(abs(X_dynare-matrix_quadratic.X)));end
AMG_Results(16,end,loop_k)=X_additional;
X_0_14=X(:,1:M_.nspred);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
AMG_Results(16,3:5,loop_k)=errors([4,7,8],1)';
catch
AMG_Results(16,3:5,loop_k)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% AMG_Results(16,:,loop_k)=NaN(1,7); 
% end
catch
AMG_Results(16,:,loop_k)=NaN(1,7);
end
run_time_reps=run_time_reps_old;

cd([YourPath])

%run different Newton methods
%newton_solvent_2(matrix_quadratic);
%newton_solvent_3(matrix_quadratic);


sprintf('Iteration %d of %d',loop_k,loop_n)
%save certain results somewhere
clearvars -except loop_k loop_n AMG_Results YourPath mmb_vec run_time_reps newton_options 
%save Second_Run_AMG_10_68
end

