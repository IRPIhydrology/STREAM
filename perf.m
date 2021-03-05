function [NS,RMSE,ANSE,RRQ,NS_lnQ,NS_radQ,KGE]=perf(Qsim,Qobs)

Qsim(isnan(Qobs))=NaN;
% Calculation of model performance
RMSE=nanmean((Qsim-Qobs).^2).^0.5/nanmean(Qobs);
NS=1-nansum((Qsim-Qobs).^2)./nansum((Qobs-nanmean(Qobs)).^2);
ANSE=1-nansum((Qobs+nanmean(Qobs)).*(Qsim-Qobs).^2)./...
    nansum((Qobs+nanmean(Qobs)).*(Qobs-nanmean(Qobs)).^2);
NS_radQ=1-nansum((sqrt(Qsim)-sqrt(Qobs)).^2)./nansum((sqrt(Qobs)-nanmean(sqrt(Qobs))).^2);
NS_lnQ=1-nansum((log(Qsim+0.00001)-log(Qobs+0.00001)).^2)...
    ./nansum((log(Qobs+0.00001)-nanmean(log(Qobs+0.00001))).^2);
X=[Qsim,Qobs]; X(any(isnan(X)'),:) = [];
RRQ=corr(Qsim,Qobs,'rows','pairwise');
KGE = klinggupta(Qsim,Qobs);
