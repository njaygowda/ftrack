function [fi,ak]=tvlptoformants_akitofi(aki,Nx,npeaks,fs)

[q,p]=size(aki);
q=q-1;
tn=0:Nx-1;

akn=zeros(p,length(tn));
for k=1:p
    for i=0:q
        akn(k,:)=akn(k,:)+aki(i+1,k)*tn.^i;
    end
end

ak = [ones(1,length(tn));-akn];

fi=zeros(length(tn),npeaks);
for i=1:length(tn)        
    fitmp=angle(roots(ak(:,i)))*fs/2/pi;
    fitmp=sort(fitmp(fitmp>0));
    if(length(fitmp)>0)
        fi(i,1:min(npeaks,length(fitmp)))=fitmp(1:min(npeaks,length(fitmp)));
    end
end

return;
