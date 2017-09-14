function M = pcm_optimalAlgorithm(M); 
% Determine the fastest fitting algorithm for each of the PCM 
% model - if not given. 
numModels = numel(M); 
for m=1:numModels 
    if (~isfield(M{m},'fitAlgorithm')) 
        switch (M{m}.type)
            case 'component' 
                if (isfield(M{m},'Gd')) 
                    M{m}.fitAlgorithm = 'NR_diag'; 
                else 
                    M{m}.fitAlgorithm = 'NR'; 
                end; 
            case 'feature'
                M{m}.fitAlgorithm = 'NR'; 
            case 'nonlinear' 
                M{m}.fitAlgorithm = 'NR';
            case 'freechol' 
                M{m}.fitAlgorithm = 'minimize';                
            case 'fixed' 
                M{m}.MM = pcm_diagonalize(M{m}.Gc);  %Model matrix for diagonalisation 
                M{m}.fitAlgorithm = 'pcm_NR_diag'; 
        end; 
    end; 
end;    % Over models 

