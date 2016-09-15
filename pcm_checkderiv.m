function [d,dh,dy] = pcm_checkderiv(f, x, e, varargin);
% function [d,dh,dy] = checkderiv(f, x, e, varargin);
% checkderiv checks the derivatives in a function, by comparing them to finite
% differences approximations. The partial derivatives and the approximation
% are printed and the norm of the diffrence divided by the norm of the sum is
% returned as an indication of accuracy.
%
% usage: checkderiv(@fcn, X, e, P1, P2, ...)
%
% where x is the argument and e is the small perturbation used for the finite
% differences. P1, P2, ... are optional additional parameters which
% get passed to f. The function f should be of the type 
%
% [fX, dfX] = f(X, P1, P2, ...)
% 
% or 
% 
% [fX, dfdx , ddfdxx] = f(X, P1, P2, ...)
% 
% where fX is the function value and dfdX is a vector of partial derivatives
% and ddfdxx is the second derviative. 
%
% Joern Diedrichsen 2015 inspred by Carl Edward Rasmussen, 2001-08-01.
[y dy] = feval(f,x,varargin{:});          % get the partial derivatives dy

for j = 1:length(x)
    xu = x; 
    xu(j) = xu(j)+e; 
    y2(j,1) = feval(f,xu,varargin{:});
    xd = x; 
    xd(j) = xd(j)-e; 
    y1(j,1) = feval(f,xd,varargin{:}); 
end
dh = (y2 - y1)./(2*e);

fprintf('analyt.  numeric\n'); 
disp([dy dh])                                          % print the two vectors
d = norm(dh-dy)/norm(dh+dy);       % return norm of diff divided by norm of sum
