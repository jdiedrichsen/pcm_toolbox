function B=blockdiag(varargin);
% function B=blockdiag(A1,A2,A3,A4,...)
% returns the block diagonal of A1.....
B=[];
for i=1:length(varargin)
    row=size(varargin{i},1);
    col=size(varargin{i},2);
    B(end+1:end+row,end+1:end+col)=varargin{i};
end;
