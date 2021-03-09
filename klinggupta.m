function [ kge,r, relvar,bias] = klinggupta(modelled,observed)
%Nash Sutcliffe Efficiency measure
modelled(isnan(observed))=NaN;

cflow=[modelled,observed];

sdmodelled=nanstd(modelled);
sdobserved=nanstd(observed);
 
mmodelled=nanmean(modelled);
mobserved=nanmean(observed);

r=corrcoef(cflow,'rows','pairwise'); r=r(1,2);
relvar=sdmodelled/sdobserved;
bias=mmodelled/mobserved;

%KGE timeseries 
kge=1-  sqrt( ((r-1)^2) + ((relvar-1)^2)  + ((bias-1)^2) );
 
end