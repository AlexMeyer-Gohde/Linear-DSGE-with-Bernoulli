load('First_run.mat')
AMG_Final_Results=AMG_Results;
AMG_Results=NaN(14,7);
AMG_Results(1:8,:)=AMG_Final_Results(1:8,:);
AMG_Results(9,:)=AMG_Final_Results(13,:);
AMG_Results(10:11,:)=AMG_Final_Results(10:11,:);
AMG_Results(12,:)=AMG_Final_Results(15,:);
AMG_Results(13:14,:)=AMG_Final_Results([14 16],:);
results=AMG_Results;
table=[results(2,[1 6 4 5]) 1];
table=[table; [ results(3:end,1)./results(2,1) results(3:end,[ 6 4 5 7])]];
disp(compose('%1.2g',table))