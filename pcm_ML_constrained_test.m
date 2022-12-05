function [P1,P2]=ML_constrained_test
% Test constrained ML estimation -
% Speed up 
X=normrnd(0,1,5,80);
A = [1 3 0 0 0 0 0 0 0 0  
     3 2 0 0 0 0 0 0 0 0  
     0 0 4 0 0 0 0 0 0 0  
     0 0 0 4 0 0 0 0 0 0  
     0 0 0 0 4 0 0 0 0 0  
     0 0 0 0 0 4 0 0 0 0  
     0 0 0 0 0 0 4 0 0 0  
     0 0 0 0 0 0 0 5 0 0  
     0 0 0 0 0 0 0 0 5 0  
     0 0 0 0 0 0 0 0 0 5  
     6 0 0 0 0 0 0 0 0 0  
     0 6 0 0 0 0 0 0 0 0  
     0 0 6 0 0 0 0 0 0 0  
     0 0 0 6 0 0 0 0 0 0  
     0 0 0 0 6 0 0 0 0 0  
     0 0 0 0 0 6 0 0 0 0  
     0 0 0 0 0 0 6 0 0 0  
     0 0 0 0 0 0 0 6 0 0  
     0 0 0 0 0 0 0 0 6 0  
     0 0 0 0 0 0 0 0 0 6];

N=100; 
[P,Q]=size(A); 
H=max(A(:)); 
for i=1:6 
    Cc{i}=double(A==i)*i; 
    sCc{i}=sparse(Cc{i}); 
end; 
% Precompute cross-terms for M-step
CcCc = cell(H, H);
for i=1:H
    for j=i:H
        CcCc{i,j}=Cc{j}'*Cc{i}; % Calculate cross terms
        sCcCc{i,j}=sparse(CcCc{i,j}); 
    end;
end;

X=normrnd(0,1,Q,N); 
for n=1:300
    Y=A*X+normrnd(0,1,P,N);
    h1=ML_constrained(X*Y',X*X',Cc);
    [c1,v1]=ML_constrained2(Y*X',X*X',Cc,CcCc);
    [c2,v2]=ML_constrained_fast(Y*X',X*X',sCc,sCcCc); 
    h2=v2\c2; 
    h3=v1\c1;
end; 
[h1 h2 h3]  


