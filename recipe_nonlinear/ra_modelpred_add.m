function [G,dGdtheta] = ra_modelpred_add(theta)
% Predicts distaces and G-matrix for additive independent
fingerParams=[1;theta(1:14)];
addParams = exp(theta(15:17));
indx1=double(tril(true(5),0));
indx1(indx1>0)=[1:15];
indx2         = indx1';
M  = [kron(ones(4,1),eye(5)) kron([0;addParams],ones(5,1))];
% Finger parameters
A=zeros(6);
for i=1:15
    d=zeros(5);
    d(indx1==i | indx2==i)=1;
    d(6,6)=0; 
    dA(:,:,i)=d; 
    A=A+dA(:,:,i)*fingerParams(i); 
end;
A(6,6)=1; 
OM = A*A'; 
for i=1:15 
    if (i>1)
        dOM = dA(:,:,i)*A'+A*dA(:,:,i)'; 
        dGdtheta(:,:,i-1)=M*dOM*M';
    end; 
end; 
G=M*OM*M';  % Second moment matrix
% Additive parameters 
for i=1:3
    add=zeros(4,1);
    add(i+1)=1;
    dM=[zeros(20,5) kron(add,ones(5,1))];
    dGdtheta(:,:,14+i)=dM*OM*M'+M*OM*dM';
    dGdtheta(:,:,14+i)=dGdtheta(:,:,14+i)*addParams(i);
end;
