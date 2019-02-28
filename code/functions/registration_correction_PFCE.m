function K1 = registration_correction_PFCE(K1,K3)
% outputs corrected K1 (translated) versus K3
L1=abs(ifftn(K1));
L3=abs(ifftn(K3));

[optimizer, metric] = imregconfig('monomodal');

[R_reg]=imregtform(L1,L3,'translation',optimizer, metric,'DisplayOptimization',0);
[dummy]=imregister(L1,L3,'translation',optimizer, metric,'DisplayOptimization',0);

shifts=R_reg.T(4,1:3);

% do the correction in k-space (on L1?)
corrx=exp(-linspace(-pi,pi,size(K1,1)).*1i.*shifts(2));
corrx=permute(corrx,[2 1]);
corry=exp(-linspace(-pi,pi,size(K1,2)).*1i.*shifts(1));
corrz=exp(-linspace(-pi,pi,size(K1,3)).*1i.*shifts(3));
corrz=permute(corrz,[1 3 2]);
K1_corr=K1;
K1_corr=bsxfun(@times,K1_corr,corrx);
K1_corr=bsxfun(@times,K1_corr,corry);
K1_corr=bsxfun(@times,K1_corr,corrz);

% check correction
L1_corr=abs(ifftn(K1_corr));
%apply correction
K1=K1_corr;