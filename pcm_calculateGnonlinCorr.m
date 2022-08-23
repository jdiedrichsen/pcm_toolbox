function [G,dGdtheta] = pcm_calculateGnonlinCorr(theta,M)
% function [G,dGdtheta] = pcm_calculateGnonlinCorr(theta,M)
% Calculate G and dGdtheta for nonlinear correlation models 
% built with pcm_buildCorrModel 
% In this paramaterization: 
% var(x) = exp(theta_x) 
% var(y) = exp(theta_y) 
% cov(x,y) = sqrt(var(x)*var(y))* r 
% r = (exp(2.*Z)-1)./(exp(2.*Z)+1);  % Fisher inverse 
if strcmp(M.r,'flexible') 
    z=theta(end);
    r=(exp(2.*z)-1)./(exp(2.*z)+1); % Correlation is the inverse fisherz transform 
    theta=theta(1:end-1); 
else 
    r=M.r; 
end 

offset = M.numCond * M.condEffect; 
% Determine the derivatives for the within-condition covariances 
indx1 = [offset+[1:M.numWithinParam]];  % Parameters for the first condition 
indx2 = [offset+M.numWithinParam+[1:M.numWithinParam]]; % Parameters for the second condition 
dGdtheta=bsxfun(@times,M.Gc,permute(exp(theta),[3 2 1])); % Get the derivatives in respect to the within-condition blocks

% Determine the predicted G-matrix so far
G=sum(dGdtheta(:,:,[indx1 indx2]),3); 

% Now determine the across-condition block
i1 = M.condVec==1; 
i2 = M.condVec==2; 
C=sqrt(G(i1,i1).*G(i2,i2)); % maximal covariance 

% Add the predicted covariances 
G(i1,i2)=C*r; 
G(i2,i1)=C'*r; 

% Now add the across-conditions blocks to the derivatives: 
for i=1:M.numWithinParam
    dC1=0.5.*1./C.*r.*G(i2,i2).*dGdtheta(i1,i1,indx1(i)); 
    dC1(C==0)=0; 
    dGdtheta(i1,i2,indx1(i))=dC1; 
    dGdtheta(i2,i1,indx1(i))=dC1'; 
    dC2=0.5.*1./C.*r.*G(i1,i1).*dGdtheta(i2,i2,indx2(i)); 
    dC2(C==0)=0; 
    dGdtheta(i1,i2,indx2(i))=dC2; 
    dGdtheta(i2,i1,indx2(i))=dC2'; 
end

% Now add the main within Condition co-variance
G=G+sum(dGdtheta(:,:,1:offset),3); 

% Add the derivative for the correlation parameter 
if strcmp(M.r,'flexible') 
    dGdtheta(i1,i2,M.numGparams)=C*4*exp(2*z)./(exp(2.*z)+1).^2;
    dGdtheta(i2,i1,M.numGparams)=C'*4*exp(2*z)./(exp(2.*z)+1).^2;
end
