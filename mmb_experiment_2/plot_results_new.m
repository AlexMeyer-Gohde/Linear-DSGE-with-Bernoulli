
load('Second_Run_AMG_improvement_new_errors.mat')
AMG_Final_Results=AMG_Results;
% load('29-99_Second_Run_AMG_improvement.mat')
% AMG_Final_Results=AMG_Results;
% load('28_Second_Run_AMG_improvement.mat')
% AMG_Final_Results(:,:,28)=AMG_Results(:,:,28);
% load('19-27_Second_Run_AMG_improvement.mat')
% AMG_Final_Results(:,:,19:27)=AMG_Results(:,:,19:27);
% load('18_Second_Run_AMG_improvement.mat')
% AMG_Final_Results(:,:,18)=AMG_Results(:,:,18);
% load('10-17_Second_Run_AMG_improvement_fill.mat')
% AMG_Final_Results(:,:,10:17)=AMG_Results(:,:,10:17);
% load('9_Second_Run_AMG_improvement.mat')
% AMG_Final_Results(:,:,9)=AMG_Results(:,:,9);
% load('1-8_Second_Run_AMG_improvement.mat')
% AMG_Final_Results(:,:,1:8)=AMG_Results(:,:,1:8);



% %index(~isnan(AMG_Results(1,1,:)));
% load('Run_AMG_old.mat')
% 
% AMG_Final_Results_old=AMG_Final_Results;
% 
% load('1_First_Run_AMG_fill.mat')
% AMG_Final_Results=NaN(16,7,loop_n);
% AMG_Final_Results(:,:,1:91)=AMG_Results(:,:,1:91);
% load('2_First_Run_AMG_fill.mat')
% AMG_Final_Results(:,:,92)=AMG_Results(:,:,92);
% load('3_First_Run_AMG_fill.mat')
% AMG_Final_Results(:,:,93:end)=AMG_Results(:,:,93:end);
% AMG_Final_Results(3:12,:,:)=AMG_Final_Results_old(3:12,:,:);
% AMG_Final_Results(2:end,6,:)=NaN;

AMG_Results=NaN(14,7,loop_n);
AMG_Results(1:8,:,:)=AMG_Final_Results(1:8,:,:);
AMG_Results(9,:,:)=AMG_Final_Results(13,:,:);
AMG_Results(10:11,:,:)=AMG_Final_Results(10:11,:,:);
AMG_Results(12,:,:)=AMG_Final_Results(15,:,:);
AMG_Results(13:14,:,:)=AMG_Final_Results([14 16],:,:);
AMG_Final_Results=NaN(14,7,loop_n);
AMG_Final_Results=AMG_Results;

model_variable_number=sum(squeeze(AMG_Final_Results(1,1:4,:)));

times=squeeze(AMG_Final_Results(2:end,1,:));
iterations=squeeze(AMG_Final_Results(2:end,7,:)); iterations(1,:)=1;
rel_times=times./times(1,:);
diffs=squeeze(AMG_Final_Results(2:end,2,:));
diffs_select=0<diffs&diffs<1e-5;
diffs_select(1,:)=1;

fe_1=squeeze(AMG_Final_Results(2:end,4,:));
rel_fe_1=fe_1./fe_1(1,:);
fe_2=squeeze(AMG_Final_Results(2:end,5,:));
rel_fe_2=fe_2./fe_2(1,:);



results(1:13,1)=sum(diffs_select')';
for j=1:13
results(j,2)=median(rel_times(j,diffs_select(j,:)),'omitnan');
results(j,3)=min(rel_times(j,diffs_select(j,:)));
results(j,4)=max(rel_times(j,diffs_select(j,:)));
results(j,5)=median(rel_fe_1(j,diffs_select(j,:)),'omitnan');
results(j,6)=min(rel_fe_1(j,diffs_select(j,:)));
results(j,7)=max(rel_fe_1(j,diffs_select(j,:)));
results(j,8)=median(rel_fe_2(j,diffs_select(j,:)),'omitnan');
results(j,9)=min(rel_fe_2(j,diffs_select(j,:)));
results(j,10)=max(rel_fe_2(j,diffs_select(j,:)));
results(j,11)=median(iterations(j,diffs_select(j,:)),'omitnan');
end

disp(compose('%1.2g',results))
compose('%1.2g',results)

save('Run_AMG_2.mat')

figure
hold on
for j=2:13
scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled')
end
legend('Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
%xlim([-2 2])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off

legends={'Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS'};

for j=2:13
    hold on
    figure
scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled')
legend(legends{j-1},'AutoUpdate','off','location','northeast')
%legend('Baseline','Modified','Samanskii','Line Search', 'Occ. Line','Occ. Line Search Samanskii','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
xlim([-4 8])
ylim([-1 3])
h=lsline;
h.Color = 'b';
h.LineStyle = '--';
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off
end

figure
hold on
for j=2:13
scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled')
end
legend('Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
%xlim([-2 2])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off


for j=2:13
    hold on
    figure
scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled');
legend(legends{j-1},'AutoUpdate','off','location','northeast')

%legend('Baseline','Modified','Samanskii','Line Search', 'Occ. Line Search','Occ. Line Sam','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
xlim([-6 6])
ylim([-1 3])
h=lsline;
h.Color = 'b';
h.LineStyle = '--';
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off
end


figure

highlight=[2 3 5 6 9];
newcolors=[0 0 0;
    0 0 1
1 0 0
0.929000000000000	0.694000000000000	0.125000000000000
0.494000000000000	0.184000000000000	0.556000000000000
0.466000000000000	0.674000000000000	0.188000000000000
0.301000000000000	0.745000000000000	0.933000000000000
0.933000000000000	0.000000000000000	0.933000000000000];
newcolors=newcolors([1 2 3 4 5 6 6 6 7 7 7 8 8],:);
lineorder={'-','-','-','-','-','-','--',':','-','--',':','-','--',};
colororder(newcolors);

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=1:13
    [y,x]=ksdensity(log10(fe_1(j,diffs_select(j,:)))); 
    if sum(highlight==j)
        plot(x,y,'LineWidth',2)
    else
        plot(x,y,'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Dynare','Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 1,  Log10')
xlim([-18 -7])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=1:13
    [y,x]=ksdensity(log10(fe_2(j,diffs_select(j,:)))); 
    if sum(highlight==j)
        plot(x,y,'LineWidth',2)
    else
        plot(x,y,'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Dynare','Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 2,  Log10')
xlim([-18 -7])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off


newcolors=newcolors(2:end,:);
lineorder=lineorder(2:end);

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:13
    [y,x]=ksdensity(log10(rel_fe_1(j,diffs_select(j,:)))); 
    if sum(highlight==j)
        plot(x,y,'LineWidth',2)
    else
        plot(x,y,'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
xlim([-5 1])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(1),'-.','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
%xl.LabelHorizontalAlignment = 'left';
%legend('show');
hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:13
    [y,x]=ksdensity(log10(rel_fe_2(j,diffs_select(j,:)))); 
    if sum(highlight==j)
        plot(x,y,'LineWidth',2)
    else
        plot(x,y,'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
xlim([-5 1])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(1),'-.','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
%xl.LabelHorizontalAlignment = 'left';
%legend('show');
hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:13
scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled')
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
%xlim([-2 2])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:13
scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled')
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Baseline','MBI','Line Search','Newton-Bernoulli', 'Newton-Bernoulli Column','Newton-Bernoulli 1/3','Newton-Bernoulli LS', 'Newton-Bernoulli Column LS','Newton-Bernoulli LS 1/3','Newton-Bernoulli Opt','Newton-Bernoulli Opt LS','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
%xlim([-2 2])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off



