function [M,CompI] = pcm_constructModelFamily(MComp,varargin)
% Constructs a model family from a set of component models
% with each model component either knocked-in or knocked-out.
% INPUT:
%       MComp:  Cell array of component models. Output will be of same type
%               as these model components
% VARARGIN:
%       'alwaysInclude': Indices of model components that are always
%               included
% OUTPUT:
%       M:      Cell array of models
%       CompI:  Component indicator for knock-in and knock-out of models
alwaysInclude = []; % If you have model components that are always included
pcm_vararginoptions(varargin,{'alwaysInclude'});

numComp = numel(MComp);
fixComp = ismember([1:numComp],alwaysInclude);
numVarComp = numComp-sum(fixComp);

Comb = dec2bin(0:(2^numVarComp-1))=='1';
i=1;
CompI= zeros(size(Comb,1),numComp);

for i=1:size(Comb,1);
    M{i}.type = MComp{1}.type;
    c=1;
    vc=1;
    M{i}.name =[];
    for j=1:numComp
        if fixComp(j) % Component is fixed (always included)
            M{i}=addComp(M{i},MComp{j},c);
            CompI(i,j)=1; 
            c=c+1;
        else
            if Comb(i,vc)
                M{i}=addComp(M{i},MComp{j},c);
                CompI(i,j)=1; 
                c=c+1;
            end;
            vc=vc+1; 
        end;
    end;
end;

function M=addComp(M,MC,c)
switch (M.type)
    case {'component','componentSq'}
        M.Gc(:,:,c)=MC.Gc;
    case 'feature'
        error('to be implemented');
end;

if (c>1)
    M.name = [M.name '+'];
end;
M.name = [M.name MC.name];
