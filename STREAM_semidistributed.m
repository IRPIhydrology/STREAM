% -------------------------------------------------------------------------------
%   STREAM: SATELLITE- BASED RUNOFF EVALUATION AND MAPPING
%   Author: Stefania Camici, Luca Brocca, Christian Massari
%   Research Institute for Geo-Hydrological Protection, National Research Council
% -------------------------------------------------------------------------------

function [NS,KGE_sez,KGE_out,Qsim_out,QBF_out,rr_tot]=STREAM_semidistributed(input,BAS_PAR,EBRR_BASPAR,PAR,sez_outlet,bas_check,ID_bas_app,FIG)

% Loading basin parameters
Nbas      = BAS_PAR(1); % number of subcatchments
Nsez      = BAS_PAR(2); % number of sections
Ninf      = BAS_PAR(3); % number of upstream inflows

basin_data = input.basin_data;
temperature= input.temperature;

delta_T=round(nanmean(diff(basin_data{1}{1}))*24*10000)/10000;
dt      = 0.1;          % computation time step in hour

% Morphological data
DIST=EBRR_BASPAR(:,2:Nsez+1);               % catchment distance to outlet sections
Ab=EBRR_BASPAR(:,2+Nsez);Ab(1:Ninf)=[];     % catchments area
A_DD=EBRR_BASPAR(:,3+Nsez);A_DD(1:Ninf)=[]; % catchments type (1: concentrater, 2: distributed)

% Initialization
M = size(basin_data{1}{1},1);
QQsim=zeros(M,Nbas);            % catchments runoff (mm/delta_T)
QQQsim=zeros(M,Nbas+Ninf,Nsez); % sections runoff (mm/delta_T)

rr_tot=cell(1,Nbas);


% Main ROUTINE
for i=1:Nbas
    
    alpha = PAR(1,i); % exponent of the infiltration 
    T1    = PAR(2,i); % characteristic time length
    gamma = PAR(3,i); % parameter of GIUH
    C     = PAR(4,i); % Celerity
    Diff  = PAR(5,i); % Diffusivity
    Ks    = PAR(6,i); % soil hydraulic conductivity
    m     = PAR(7,i); % exponent of drainage
    Cm    = PAR(8,i); % degree-day coefficient for snow module
    
    
    D=basin_data{i}{:,1};
    PIO_ =  basin_data{i}{:,2};
    try, for s=1:size(PIO_,2),
            PIO_(:,s)=interp1(D(~isnan(PIO_(:,s))),PIO_(~isnan(PIO_(:,s)),s),D,'linear','extrap');
        end
    end
    SM  =  basin_data{i}{:,3};
    try, for s=1:size(SM,2),
            SM(:,s)=interp1(D(~isnan(SM(:,s))),SM(~isnan(SM(:,s)),s),D,'linear','extrap');
        end
    end
    TEMPER = temperature{i};
    try, for s=1:size(TEMPER,2),
            TEMPER(:,s)=interp1(D(~isnan(TEMPER(:,s))),TEMPER(~isnan(TEMPER(:,s)),s),D,'linear','extrap');
        end
    end
    TWSA = basin_data{i}{:,5};
    try, for s=1:size(TWSA,2),
            TWSA(:,s)=interp1(D(~isnan(TWSA(:,s))),TWSA(~isnan(TWSA(:,s)),s),D,'linear','extrap');
        end
    end
    
    
    % Snow Module
    [PIO,SWE]=snow_model(PIO_, TEMPER, -0.5, 0.5, Cm);
    
    %  Soil Water Index comnputation
    SWI = SWIcomp_NAN([D,SM],T1);
    W1 =(SWI-nanmin(SWI))./((nanmax(SWI)-nanmin(SWI)));
    
    W2=(TWSA-nanmin(TWSA))./((nanmax(TWSA)-nanmin(TWSA)));
    
    peff=zeros(size(PIO));             % effective rainfall (mm/delta_T)
    BF=zeros(size(PIO));
    rr_tot{i}=zeros(size(PIO));
    
    for jj=1:size(PIO,2)
        peff(:,jj)=PIO(:,jj).*(W1(:,jj).^alpha); % Georgakakos and Baumer (1996)
        BF  (:,jj)=Ks.*(W2(:,jj)).^(m);
    end
    
    rr_tot {i}= (peff+BF);
    peffmean = nanmean(peff,2);
    BF_mean  = nanmean(BF,2);
    
    % Convolution (GIUH and NASH)
    if Ab(i)>0.0
        if A_DD(i)==1
            IUH=GIUH(gamma,Ab(i),dt,delta_T)*dt;
        else
            IUH=NASH(2*gamma,Ab(i),dt,delta_T,1)*dt;
        end
        % interpolation for dt time step
        peffint=interp1(1:M,peffmean,1:dt:M)';
        BFint=interp1(1:M,BF_mean,1:dt:M)';
        % convolution
        temp1=conv(IUH,peffint);
        temp2=conv(NASH(2*gamma,Ab(i),dt,delta_T,1)*dt,BFint);
        % saving outputs
        QQBF(:,i)= temp2(1:round(1/dt):M*round(1/dt));
        QQsim(:,i)= temp1(1:round(1/dt):M*round(1/dt))...
            + temp2(1:round(1/dt):M*round(1/dt)); % in delta_T
    else % if very low area (< 0 km^2) no convolution
        QQBF(:,i) = BF_mean;
        QQsim(:,i)= peffmean+BF_mean;
    end
    
    % Convolution (Hayami)
    for j=1:Nsez
        if DIST(i+Ninf,j)>0
            g=hayami(dt,DIST(i+Ninf,j),C,Diff,delta_T)*dt;
            % interpolation for dt time step
            QQsimint=interp1(1:M,QQsim(:,i),1:dt:M)';
            QQBFint=interp1(1:M,QQBF(:,i),1:dt:M)';
            % convolution
            temp1=(conv(g',QQsimint));
            temp2=(conv(g',QQBFint));
            % saving outputs
            QQQBF(:,i+Ninf,j)=temp2(1:round(1/dt):M*round(1/dt)); % in delta_T
            QQQsim(:,i+Ninf,j)=temp1(1:round(1/dt):M*round(1/dt));
        elseif DIST(i+Ninf,j)==0
            QQQBF(:,i+Ninf,j)=QQBF(:,i);
            QQQsim(:,i+Ninf,j)=QQsim(:,i);
        else % if negative distance no contribution
            QQQBF(:,i+Ninf,j)=zeros(M,1,1);
            QQQsim(:,i+Ninf,j)=zeros(M,1,1);
        end
    end
end

% Calculation of catchments outlet discharge in m^3/s
for i=1:Nbas
    QQBF(:,i)= QQBF(:,i).*(Ab(i)./delta_T./3.6);
    QQsim(:,i)=QQsim(:,i).*(Ab(i)./delta_T./3.6);
end

% Calculation of sections discharge in m^3/s
for i=1:Nbas
    for j=1:Nsez
        QQQBF(:,i+Ninf,j)=QQQBF(:,i+Ninf,j).*(Ab(i)./delta_T./3.6);
        QQQsim(:,i+Ninf,j)=QQQsim(:,i+Ninf,j).*(Ab(i)./delta_T./3.6);
        
    end
end


% Convolution (Hayami) of upstream discharge
for j=1:Nsez
    for i=1:Ninf
        if DIST(i,j)>0.0
            g=hayami(dt,DIST(i,j),C,Diff,delta_T)*dt;
            QQsimint=interp1(1:M,QM(:,i),1:dt:M)';
            temp1=conv(g',QQsimint);
            QQQsim(:,i,j)=temp1(1:round(1/dt):M*round(1/dt));
        end
    end
end

% Calculation of outlet sections discharge
for j=1:Nsez
    QBF1(:,j)=(nansum(QQQBF(:,:,j)'))';
    Qsim1(:,j)=(nansum(QQQsim(:,:,j)'))';
    KGE_sez(j,1) =klinggupta(Qsim1(:,j),basin_data{ID_bas_app(j)}{:,4});
end
Qsim=Qsim1(:,sez_outlet);
QBF_out=[QBF1,QQBF];
Qsim_out=[Qsim1,QQsim];


% identification Qobs
Qobs =  basin_data{bas_check}{:,4};

% Calculation of model performance
RMSE=nanmean((Qsim-Qobs).^2).^0.5;
NS=1-nansum((Qsim-Qobs).^2)./nansum((Qobs-nanmean(Qobs)).^2);
ANSE=1-nansum((Qobs+nanmean(Qobs)).*(Qsim-Qobs).^2)./...
    nansum((Qobs+nanmean(Qobs)).*(Qobs-nanmean(Qobs)).^2);
NS_radQ=1-nansum((sqrt(Qsim)-sqrt(Qobs)).^2)./nansum((sqrt(Qobs)-nanmean(sqrt(Qobs))).^2);
NS_lnQ=1-nansum((log(Qsim+0.00001)-log(Qobs+0.00001)).^2)...
    ./nansum((log(Qobs+0.00001)-nanmean(log(Qobs+0.00001))).^2);
X=[Qsim,Qobs]; X(any(isnan(X)'),:) = [];
RRQ=corrcoef(X).^2; RQ=RRQ(2);
KGE_out=klinggupta(Qsim,Qobs);

% -------------------------------------------------------------------------------
% Calculation of Geomorphological Instantaneous Unit Hydrograph
% -------------------------------------------------------------------------------
function IUH=GIUH(gamma,Ab,dt,deltaT)

Lag=(gamma*1.19*Ab^0.33)/deltaT;
hp=0.8/Lag;
data=load('GIUH');
t=data(:,1)*Lag;IUH_0=data(:,2)*hp;
ti=0:dt:max(t);
IUH=interp1(t,IUH_0,ti)';
IUH=IUH./(sum(IUH).*dt);

% -------------------------------------------------------------------------------
% Calculation of Nash Instantaneous Unit Hydrograph
% -------------------------------------------------------------------------------
function IUH=NASH(gamma,Ab,dt,deltaT,n)

K=(gamma*1.19*Ab^0.33)/deltaT;
Tmax=100; % (h)
time=0.00001:dt:Tmax;
IUH=((time/K).^(n-1).*exp(-time/K)/factorial(n-1)/K)';
IUH=IUH./(sum(IUH).*dt);

% -------------------------------------------------------------------------------
% Calculation of Hayami function (diffusive routing)
% -------------------------------------------------------------------------------
function g=hayami(dt,L,C,D,deltaT)

C=C*deltaT;D=D*deltaT;
Tmax=100; % (h)
tt=0.00001:dt:Tmax;
g=(L./sqrt(4*pi*D.*tt.^3)).*exp(-((L-C.*tt).^2)./(4*D.*tt));
g=g./(sum(g).*dt);

%--------------------------------------------------------------------------
% Snow accumulation-melting MODEL
%--------------------------------------------------------------------------

function [rainfall,SWE_melting,SWE_snowpack]=snow_model(precipitation, temperature, temp_min, temp_max, Cm)

rainfall = zeros(size(precipitation));
snowfall = zeros(size(precipitation));
SWE_snowpack = zeros(size(precipitation));
SWE_melting = zeros(size(precipitation));

% The precipitation is divided into rainfall and snowfall
% REFERENCES:
% U.S. Army Corps of Engineers (1956)
% Boscarello, L., Ravazzani, G., Pellegrini, M., Dedieu, J. P., & Mancini, M. (2014). Calibration of hydrological model FEST from MODIS images in Alpine Catchments. Politecnico di Milano, Dipartimento di Ingegneria Idraulica, Ambientale, Infrastrutture viarie, Rilevamento.
% Degree Day Method (Mockus, 1964)

% INITIALIZATION
for jj=1:size(precipitation,2)
    if precipitation(1,jj) == NaN || temperature(1,jj) == NaN
        rainfall(1,jj) = NaN;
        snowfall(1,jj) = NaN;
    elseif temperature(1,jj) <= temp_min
        snowfall(1,jj) = precipitation(1,jj);
        rainfall(1,jj) = 0;
        SWE_snowpack(1,jj) = snowfall(1,jj); % [mm]
        SWE_melting(1,jj) = 0; % [mm]
    elseif temperature(1,jj) >= temp_max
        snowfall(1,jj) = 0;
        rainfall(1,jj) = precipitation(1,jj);
        SWE_snowpack(1,jj) = 0; % [mm]
        SWE_melting(1,jj) = 0; % [mm]
    else
        rainfall(1,jj) = precipitation(1,jj) * ((temperature(1,jj)-temp_min)/(temp_max-temp_min));
        snowfall(1,jj) = precipitation(1,jj) - rainfall(1,jj);
        SWE_snowpack(1,jj) = snowfall(1,jj);
        SWE_melting(1,jj) = 0;
    end
    
    
    for i=2:length(precipitation)
        if precipitation(i,jj) == NaN || temperature(i,jj) == NaN
            % if the data is missing, it is equal to NaN
            rainfall(i,jj) = NaN;
            snowfall(i,jj) = NaN;
        elseif temperature(i,jj) <= temp_min
            % if the temperature is less than the low threshold,
            % the precipitation is entirely snowfall
            rainfall(i,jj) = 0;
            snowfall(i,jj) = precipitation(i,jj);
            SWE_snowpack(i,jj) = SWE_snowpack(i-1,jj) + snowfall(i,jj);
            SWE_melting(i,jj) = 0;
        elseif temperature(i,jj) > temp_max
            % if the temperature is more than the high threshold,
            % the precipitation is entirely rainfall
            rainfall(i,jj) = precipitation(i,jj);
            snowfall(i,jj) = 0;
            SWE_melting(i,jj) = Cm * (temperature(i,jj) - temp_max);
            % h_melting(i,1) = rho_water * SWE_melting(i,1) / rho_snow;
            % Check the snowpack SWE
            if SWE_snowpack(i-1,jj) >= SWE_melting(i,jj)
                SWE_snowpack(i,jj) = SWE_snowpack(i-1,jj) - SWE_melting(i,jj);
            else
                SWE_melting(i,jj) = SWE_snowpack(i-1,jj);
                SWE_snowpack(i,jj) = 0;
            end
        else
            rainfall(i,jj) = precipitation(i,jj) * ((temperature(i,jj)-temp_min)/(temp_max-temp_min));
            snowfall(i,jj) = precipitation(i,jj) - rainfall(i,jj);
            SWE_snowpack(i,jj) = SWE_snowpack(i-1,jj) + snowfall(i,jj);
            SWE_melting(i,jj) = 0;
        end
    end
end


%--------------------------------------------------------------------------
%  SWI computation, see Albergel et al., 2008 (HESS)
%--------------------------------------------------------------------------
function SSWI=SWIcomp_NAN(data,T)
dt=zeros(length(data),1);
Kn=ones(length(data),1);
SWI=data(:,2:end);
D=data(:,1);
SSM=data(:,2:end);
SSWI=nan(length(data(:,1)),size(data(:,2:end),2));

for j= 1:size(SWI,2)
    D_temp = D;
    D_temp(isnan(SWI(:,j)))=[];
    SSM_temp = SSM(:,j);
    SSM_temp(isnan(SWI(:,j)))=[];
    SWI_temp = SWI(:,j);
    SWI_temp(isnan(SWI(:,j)))=[];
    
    for i=2:length(SWI_temp)
        dt=(D_temp(i)-D_temp(i-1));
        Kn(i,j)=Kn(i-1)./(Kn(i-1)+exp(-dt/T));
        SWI_temp(i,1)=SWI_temp(i-1)+Kn(i,j).*(SSM_temp(i,1)-SWI_temp(i-1));
    end
    
    ID_notNaN=find(~isnan(SWI(:,j)));
    SSWI(ID_notNaN,j)=SWI_temp;
    clear SWI_temp SSM_temp D_temp
    SWIout{j}=[data(:,1),SSWI(:,j),data(:,2)];
end

