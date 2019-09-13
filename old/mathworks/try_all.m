% Test code for empca with missing data and weights
% 
% 2017 Vicente Parot
% Cohen Lab
% Harvard University

%% generate a test matrix
%rng(0,'twister') % for reproducibility
obj = poissrnd((phantom+1)*50);
obj = obj(:,76:180);
nan_fraction = .5; % 50% of values will be set to NaN
w = rand(size(obj)); % weights are uniform ~U(0,1)
obj_0 = obj;
obj_0(w<nan_fraction) = 0;
obj_missing = obj;
% obj_missing(ceil(rand(ceil(nan_fraction*end),1)*end)) = NaN;
obj_missing(w<nan_fraction) = nan;

  figure
    imshow(obj,[])
    hold on
    [yloc,xloc] = find(isnan(obj_missing));
    plot(xloc,yloc,'.r')
    title 'object with missing data'
    drawnow

%% calculate empca and svd
ncs = 15;
[u1, s1, v1] = svd(obj,'econ');
[u2, s2, v2] = empca(obj,ncs);
[u3, s3, v3] = empca_w(obj,w,ncs);
[u4, s4, v4] = empca_w(obj,~isnan(obj_missing),ncs);
[u5, s5, v5] = empca_w(obj_0,~isnan(obj_missing),ncs);
[u6, s6, v6] = empca_nan(obj_missing,ncs);

disp([
    diag(s1(1:5,1:5))'
    diag(s2(1:5,1:5))'
    diag(s3(1:5,1:5))'
    diag(s4(1:5,1:5))'
    diag(s5(1:5,1:5))'
    diag(s6(1:5,1:5))'
    ]') % print a few eigenvalues
%%

figure
imshow([... display the data modeled by the calculated vectors
    obj ...
    u1(:,1:ncs)*s1(1:ncs,1:ncs)*v1(:,1:ncs)' ...
    u2(:,1:ncs)*s2(1:ncs,1:ncs)*v2(:,1:ncs)' ...
    u3(:,1:ncs)*s3(1:ncs,1:ncs)*v3(:,1:ncs)' ...
    u4(:,1:ncs)*s4(1:ncs,1:ncs)*v4(:,1:ncs)' ...
    u5(:,1:ncs)*s5(1:ncs,1:ncs)*v5(:,1:ncs)' ...
    u6(:,1:ncs)*s6(1:ncs,1:ncs)*v6(:,1:ncs)'],[])
title({'data','obj | truncated svd | empca | weighted empca | weighted missing empca| weighted zeros empca | missing data empca'})

figure
imshow([... display the data modeled by the calculated vectors
    obj - obj ...
    obj - u1(:,1:ncs)*s1(1:ncs,1:ncs)*v1(:,1:ncs)' ...
    obj - u2(:,1:ncs)*s2(1:ncs,1:ncs)*v2(:,1:ncs)' ...
    obj - u3(:,1:ncs)*s3(1:ncs,1:ncs)*v3(:,1:ncs)' ...
    obj - u4(:,1:ncs)*s4(1:ncs,1:ncs)*v4(:,1:ncs)' ...
    obj - u5(:,1:ncs)*s5(1:ncs,1:ncs)*v5(:,1:ncs)' ...
    obj - u6(:,1:ncs)*s6(1:ncs,1:ncs)*v6(:,1:ncs)'],[])
title({'data diff','obj | truncated svd | empca | weighted empca | weighted missing empca| weighted zeros empca | missing data empca'})

figure
imshow([... display the data modeled by the calculated vectors
    obj*obj' ...
    u1*s1^2*u1' ...
    u2*s2^2*u2' ...
    u3*s3^2*u3' ...
    u4*s4^2*u4' ...
    u5*s5^2*u5' ...
    u6*s6^2*u6'],[])
title({'left covariance','obj | truncated svd | empca | weighted empca | weighted missing empca| weighted zeros empca | missing data empca'})

figure
imshow([... display the data modeled by the calculated vectors
    obj'*obj ...
    v1*s1^2*v1' ...
    v2*s2^2*v2' ...
    v3*s3^2*v3' ...
    v4*s4^2*v4' ...
    v5*s5^2*v5' ...
    v6*s6^2*v6'],[])
title({'right covariance','obj | truncated svd | empca | weighted empca | weighted missing empca| weighted zeros empca | missing data empca'})


