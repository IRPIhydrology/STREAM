function [PAR]=cal_STREAM_semidistributed(input,BAS_PAR,EBRR_BASPAR,sez_outlet,bas_check,ID_bas_app,X_ini)
NPAR=8;
Nbas      = BAS_PAR(1); % number of subcatchments

if nargin<7, X_ini=ones(NPAR,1)*.5; end

[RES,FVAL]=fmincon(@calibOK,X_ini,[],[],[],[],...
    zeros(NPAR,1),ones(NPAR,1),[],optimset('Display','iter','MaxIter',300,...
     'MaxFunEvals',500,'TolFun',1E-6,'TolCon',6,...
     'Largescale','off','Algorithm','active-set'),input,BAS_PAR,EBRR_BASPAR,sez_outlet,bas_check,ID_bas_app);
X_OPT=convert_adim(RES,1);
PAR=repmat(X_OPT',Nbas,1)';

%---------------------------------------------------------------------------------
function [errore]=calibOK(X_0,input,BAS_PAR,EBRR_BASPAR,sez_outlet,bas_check,ID_bas_app)
Nbas      = BAS_PAR(1); % number of subcatchments

X=convert_adim(X_0,1);
PAR=repmat(X,1,Nbas);

[NS,KGE_sez,KGE_out]=STREAM_semidistributed(input,BAS_PAR,EBRR_BASPAR,PAR,sez_outlet,bas_check,ID_bas_app,1);
errore=1-KGE_out;
%---------------------------------------------------------------------------------
function X=convert_adim(X_0,Nbas)
%     alpha    T   gamma  C   D  beta    m     Cm
LOW=[   1.0  0.01   0.5   1   1   0.1     1   0.1/24]';
UP =[  30.0    80   5.5  60  30    20    15     3   ]';


LOW =repmat(LOW,1,Nbas);
UP =repmat(UP,1,Nbas);

X=LOW+(UP-LOW).*X_0;


