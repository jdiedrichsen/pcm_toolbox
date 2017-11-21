function M=pcm_addModelComp(M,MC)
% function M=pcm_addModelComp(M,MC)
% Add the model component in MC to the model specified in M 
% Adds the component automatically to the model name 
% INPUT: 
%     M: Model structure to be added to 
%     MC: Model component to be added 
% OUTPUT: 
%     M: Amended model structure 
switch (M.type)
    case {'component','componentSq'}
        if (~isfield(M,'Gc') || isempty(M.Gc))
            c=1;
        else 
            c=size(M.Gc,3)+1; 
        end;
        M.Gc(:,:,c)=MC.Gc;
    case 'feature'
        error('to be implemented');
end;

if (c>1)
    M.name = [M.name '+'];
end;
M.name = [M.name MC.name];
