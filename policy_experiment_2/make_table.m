clear all;


load Combined_Policy_Run.mat

results_table=zeros(10,9,length(STEPS));

%grid size
for k=1:length(STEPS)

%algorithm type
for j=1:10
%run times
results_table(j,1,k)=median(median(total_AMG_Results{1,k}(j+1,1,:,:)./total_AMG_Results{1,k}(2,1,:,:)))
results_table(j,2,k)=min(min(total_AMG_Results{1,k}(j+1,1,:,:)./total_AMG_Results{1,k}(2,1,:,:)));
results_table(j,3,k)=max(max(total_AMG_Results{1,k}(j+1,1,:,:)./total_AMG_Results{1,k}(2,1,:,:)));
%fe1
results_table(j,4,k)=median(median(total_AMG_Results{1,k}(j+1,4,:,:)./total_AMG_Results{1,k}(2,4,:,:)));
results_table(j,5,k)=min(min(total_AMG_Results{1,k}(j+1,4,:,:)./total_AMG_Results{1,k}(2,4,:,:)));
results_table(j,6,k)=max(max(total_AMG_Results{1,k}(j+1,4,:,:)./total_AMG_Results{1,k}(2,4,:,:)));
%fe2
results_table(j,7,k)=median(median(total_AMG_Results{1,k}(j+1,5,:,:)./total_AMG_Results{1,k}(2,5,:,:)));
results_table(j,8,k)=min(min(total_AMG_Results{1,k}(j+1,5,:,:)./total_AMG_Results{1,k}(2,5,:,:)));
results_table(j,9,k)=max(max(total_AMG_Results{1,k}(j+1,5,:,:)./total_AMG_Results{1,k}(2,5,:,:)));
end

end

save('results_table.mat','results_table');