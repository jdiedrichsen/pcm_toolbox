function [M,CompI] = pcm_constructModelFamily(MComp,varargin)
% Constructs a model family from a set of component models
% with each model component either knocked-in or knocked-out.
% INPUT:
%       MComp:  Cell array of component models. Output will be of same type
%               as these model components
% VARARGIN:
%       'alwaysInclude': Indices of model components that are always
%               included
%       'fullModel': The model input is a full model, and should be broken
%                    into components if 'fullModel'==1
%                    only works if full model is a component model 
%                    (has MComp.Gc field)
% OUTPUT:
%       M:      Cell array of models
%       CompI:  Component indicator for knock-in and knock-out of models
alwaysInclude = []; % If you have model components that are always included
fullModel =0; 

pcm_vararginoptions(varargin,{'alwaysInclude','fullModel'});


% Check if input is full model
if fullModel==1
    % break down into components
    numComp = size(MComp.Gc,3);
    for n=1:numComp
        MTemp{n}.Gd = MComp.Gd(n,:);
        MTemp{n}.Gc = MComp.Gc(:,:,n);
        MTemp{n}.name = sprintf('Comp%d',n);
        MTemp{n}.type = MComp.type;
    end
    clear MComp;
    MComp = MTemp;
end

% Get number of fixed and variable components 
numComp = numel(MComp);
fixComp = ismember([1:numComp],alwaysInclude);
numVarComp = numComp-sum(fixComp);

% Build all combination of 0,1,2... components 
Comb = zeros(1,numVarComp); 
for i=1:numVarComp 
    A=nchoosek([1:numVarComp],i);
    n = size(A,1); 
    X=zeros(n,numVarComp); 
    for j=1:n
         X(j,A(j,:))=1;
    end; 
    Comb=[Comb;X]; 
end; 

% Now build the models with fixed components set in 
CompI= zeros(size(Comb,1),numComp);
i=1;
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
        M{i}.numGparams=0; % Empty Model 
        M{i}.name = 'null'; 
        M{i}.Gc   = zeros(size(MComp{1}.Gc)); 
    else  
        M{i}.numGparams=size(M{i}.Gc,3);  
    end; 
end;
