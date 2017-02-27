function M=pcm_prepFreeModel(M); 
% function M=pcm_prepFreeModel(M); 
% Prepares the free model by calculating vectorization of Cholensky matrix 
% Speeds up subsequent computation 
if (~isfield(M,'numCond'))
    error('need to provide numCond as a field'); 
end; 
M.indx         = find(tril(true(M.numCond))); 
[M.row,M.col]  = ind2sub([M.numCond M.numCond],M.indx);
M.numGparams   = length(M.indx); 