%clear all
%load('First_Run_AMG_JS_test_ls.mat')
results=AMG_Final_Results(:,:,97);
table=[results(2,[1 2 4 5]) 1];
table=[table; [ results(3:end,1)./results(2,1) results(3:end,[ 2 4 5 7])]];
disp(compose('%1.2g',table))