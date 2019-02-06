clear;clc;
w = csvread('../madeleine_data/master_qlearning_K_T1_power_RUNS.csv',1);

clc;
i=1;
A.nI=50;
A.nT=17;
A.nJ=7;

A.choice=nan(A.nI,A.nT,A.nJ);
A.bust=nan(A.nI,A.nT,A.nJ);
A.RT=nan(A.nI,A.nT,A.nJ);
A.HE=nan(A.nI,A.nT,A.nJ);
A.PCPP=nan(A.nI,A.nT,A.nJ);
A.risky=nan(A.nI,A.nT,A.nJ);
A.safe=nan(A.nI,A.nT,A.nJ);

for x=1:length(w)
if x>1
if w(x,1)~=w(x-1,1)
    i=i+1;
    A.UID(i) = w(x,1);
end
end
A.risky(i,w(x,7),w(x,3))=w(x,4);
A.bust(i,w(x,7),w(x,3))=w(x,2);
A.HE(i,w(x,7),w(x,3))=w(x,8);
A.RT(i,w(x,7),w(x,3))=w(x,6)./1000;
A.PCPP(i,w(x,7),w(x,3))=w(x,9);
end
for i=1:A.nI
A.minT(i)=min(min(A.RT(i,:,:)));
end
A.RTACC = A.RT.*(2.*A.risky-1);
A.MX=7;
A.safe=1-A.risky;

for ii=1:A.nI
A.cummBust(ii,:)=cumsum(nansum(A.bust(ii,:,:),3));
end
A.CI = nansum(A.bust,3);

ccD=A;
save('ccD','ccD');
