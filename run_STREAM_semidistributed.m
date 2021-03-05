%% crea file input
clear, clc, close all
name_bas= 'Mississippi'
TOPOL = load(['topology_',name_bas,'.csv']);
TOPOL(TOPOL(:,4)==0,4)=2;
TOPOL(TOPOL(:,4)==-1,4)=1;

%% load input data
load (['staz_check',name_bas,'.mat'])
load input.mat

%% File preparation for basin configuration 
EBRR_BASPAR = [TOPOL(:,1),load(['distance_',name_bas,'.txt']),TOPOL(:,6),TOPOL(:,4)];
BAS_PAR =[size(TOPOL,1); size(ID_bas_app,1); 0];

%% identification closure section according to basin configuration in "Mississippi_basin.png"
sez_outlet = 10; % Vicksburg 
bas_check = ID_bas_app(sez_outlet);
 
%% STREAM model calibration   
[X_OPT]=cal_MISDC_SEMIDISTR_basflow_pardistrGRACE(input,BAS_PAR,EBRR_BASPAR,sez_outlet,bas_check,ID_bas_app)    
save (['X_opt_',name_bas,'.txt'],'X_OPT','-ascii','-tabs')

%% STREAM model run
X_OPT= load(['X_opt_',name_bas,'.txt']);

[NS,KGE_sez,KGE_out,Qsim_out,QB_out,rr_tot]=STREAM_semidistributed(input,BAS_PAR,EBRR_BASPAR,X_OPT,sez_outlet,bas_check,ID_bas_app,1);

%% River discharge comparison for calibrated sections 
basin_data = input.basin_data;
temperature= input.temperature;
sez_check= [9,6,4,11];

close all
ID_bas_checkinn=ID_bas_app(sez_check); 

set(gcf,'position',[576    49   777   936])
subplot(5,1,1)
Qobs =  basin_data{bas_check}{:,4};

QQ=Qsim_out(:,sez_outlet);
QQB = QB_out(:,sez_outlet);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])

datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outlet),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])

subplot(5,1,2)
sez_outletinn = find(ID_bas_app==ID_bas_checkinn(1));
Qobs =  basin_data{ID_bas_checkinn(1)}{:,4};
QQ=Qsim_out(:,sez_outletinn);
QQB = QB_out(:,sez_outletinn);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])
datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outletinn),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])



subplot(5,1,3)
sez_outletinn = find(ID_bas_app==ID_bas_checkinn(2));
Qobs =  basin_data{ID_bas_checkinn(2)}{:,4};
QQ=Qsim_out(:,sez_outletinn);
QQB = QB_out(:,sez_outletinn);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])
datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outletinn),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])

subplot(5,1,4)
sez_outletinn = find(ID_bas_app==ID_bas_checkinn(3));
Qobs =  basin_data{ID_bas_checkinn(3)}{:,4};
QQ=Qsim_out(:,sez_outletinn);
QQB = QB_out(:,sez_outletinn);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])
datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outletinn),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])

subplot(5,1,5)
sez_outletinn = find(ID_bas_app==ID_bas_checkinn(4));
Qobs =  basin_data{ID_bas_checkinn(4)}{:,4};
QQ=Qsim_out(:,sez_outletinn);
QQB = QB_out(:,sez_outletinn);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])
datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outletinn),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])
saveas(gcf,'Results\Calibrated_sections.png')
%% River discharge comparison for non calibrated sections 
basin_data = input.basin_data;
temperature= input.temperature;
sez_check= [1,2,8];


close all
ID_bas_checkinn=ID_bas_app(sez_check); 

set(gcf,'position',[576    49   777   936])
subplot(4,1,1)
Qobs =  basin_data{bas_check}{:,4};
QQ=Qsim_out(:,sez_outlet);
QQB = QB_out(:,sez_outlet);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])

datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outlet),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])


subplot(4,1,2)
sez_outletinn = find(ID_bas_app==ID_bas_checkinn(1));
Qobs =  basin_data{ID_bas_checkinn(1)}{:,4};
QQ=Qsim_out(:,sez_outletinn);
QQB = QB_out(:,sez_outletinn);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])
datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outletinn),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])


subplot(4,1,3)
sez_outletinn = find(ID_bas_app==ID_bas_checkinn(2));
Qobs =  basin_data{ID_bas_checkinn(2)}{:,4};
QQ=Qsim_out(:,sez_outletinn);
QQB = QB_out(:,sez_outletinn);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])
datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
% legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outletinn),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])

subplot(4,1,4)
sez_outletinn = find(ID_bas_app==ID_bas_checkinn(3));
Qobs =  basin_data{ID_bas_checkinn(3)}{:,4};
QQ=Qsim_out(:,sez_outletinn);
QQB = QB_out(:,sez_outletinn);

plot(basin_data{1}{1},Qobs,'g-','linew',2)
hold on
plot(basin_data{1}{1},QQ,'r--')
plot(basin_data{1}{1},QQB,'b-')
axis([(basin_data{1}{1}(1)) (basin_data{1}{1}(end)) min(min(QQ,Qobs)) max(max(QQ,Qobs))])
datetick('x','keeplimits')
box on, grid on
ylabel('Discharge (m^3/s)')
% legend('Observed data','STREAM total discharge','STREAM slow flow')
[~,RMSE,~,RRQ]=perf(QQ,Qobs);
title(['Outlet section: ',num2str(sez_outletinn),' - KGE: ',num2str(klinggupta(QQ,Qobs),'%3.2f'),...
    ' | R: ',num2str(RRQ,'%3.2f'), ' | RRMSE: ',num2str(RMSE*100,'%3.2f'),'%'])
saveas(gcf,'Results\Notcalibrated_sections.png')

