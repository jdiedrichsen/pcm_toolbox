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

% Get number of fixed and variable components 
numComp = numel(MComp);
fixComp = ismember([1:numComp],alwaysInclude);
numVarComp = numComp-sum(fixComp);

Comb = dec2bin(0:(2^numVarComp-1))=='1';
i=1;
CompI= zeros(size(Comb,1),numComp);

for i=1:size(Comb,1);
    M{i}.type = MComp{1}.type;
    vc=1;
    M{i}.name =[];
    M{i}.Gc   =[]; 
    for j=1:numComp
        if fixComp(j) % Component is fixed (always included)
            M{i}=pcm_addModelComp(M{i},MComp{j});
            CompI(i,j)=1; 
        else
            if Comb(i,vc)
                M{i}=pcm_addModelComp(M{i},MComp{j});
                CompI(i,j)=1; 
            end;
            vc=vc+1; 
        end;
    end;
    if isempty(M{i}.Gc)
        M{i}.numGparams=0 % Empty Model 
        M{i}.name = 'null'; 
        M{i}.Gc   = zeros(size(MComp{1}.Gc)); 
    else  
        M{i}.numGparams=size(M{i}.Gc,3);  
    end; 
end;
