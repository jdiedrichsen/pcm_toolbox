function varargout = pcm_searchlight_run(Searchlights,ImagesIn,ImagesOut,pcmfcn,...
                                         params,varargin)
%% function varargout = pcm_searchlight_run(Searchlights,Images,pcmfcn,varargin)
% Run pcm on the surface for multiple subjects from raw time series image.
% 
% INPUTS:
%   Searchlights: Searchlight definitions for multiple subjects. 
% 
%   ImagesIn:       Input images for multiple subjects.
% 
%   ImagesOut:      Output images for pcmfcn.
% 
%   pcmfcn:         user-difined pcm model fitting function (e.g., )
%                   1st input argument should be;
%                   Y: a cell array of Ndata x Nvoxel data
% 
%   params:         all the other parameters for pcmfcn
% 
% OPTIONAL INPUTS:
%   'verbose':      if set to 1, inhibit commandline output.
% 
% 
% 
% OUTPUTS:
% 
% 
% 
% 
% 
% JD & AY 2016
% 
% Recommended spec: memory larger than ** GB?

% note: good to be compatible for both one-step (prewhitening+fitting) and 
%       two-step (prewhitening+save, then fitting) style


%% Input check (cell/structure/char)
% Options
verbose     = 0;       % Control output to commandline 
tmpDir      = {};      % Directory for tmp*.mat files (for each subject) 
columnNames = {};
vararginoptions(varargin,{'tmpDir','verbose','columnNames'});

%-------------------------------------------------------------------------%
% Checnk input
%-------------------------------------------------------------------------%
% Number of subjects
sn = numel(Searchlights);
if sn~=numel(ImagesIn);
    error('Size mismatch between Searchlight and ImagesIn!\n');
end

% Filetype for Searchlights (cell array of file names? structure array?)
% Check searchlights
if iscell(Searchlights)
    SLnames = Searchlights;
    clear Searchlights
    for s=1:sn
        Searchlights(s) = load(SLnames{s});
    end
end
Nnode1      = numel(Searchlights(1).LI); % assuming the same node size across sbjs
NsubjNode   = zeros(Nnode1,1);
for s=1:sn    
    Searchlights(s).Good = ~isnan(Searchlights(s).vORr); % non-empty node
    NsubjNode = NsubjNode + double(Searchlights(s).Good);
end
% Get node id
if ~isfield(Searchlights(1),'nodeID');
    for s=1:sn
        Searchlights(s).nodeID      = [1:Nnode1]';
        Searchlights(s).num_nodes   = Nnode1;
    end
end

%-------------------------------------------------------------------------%
% Check output
%-------------------------------------------------------------------------%
Nout = numel(ImagesOut);
if (~isempty(ImagesOut))&&(~iscell(ImagesOut))
    error('outfiles must be given in a cell array.'); 
end;
for i=1:Nout
    [d f t m]       = spm_fileparts(ImagesOut{i}(1,:));
    metricNames{i}  = fullfile(d, [f t]);
    metricColumn{i} = m(2:end);
    if (~strcmp(t,'.metric'))
        warning('Output file extension must be .metric.');
    end;
end
%=== currently only allow .metric file type for output?
% Check output size
for s=1:sn
    X{s} = rand(length(ImagesIn{1}),10);
end
Result  = feval(pcmfcn,X,params{:});
Ncol    = size(Result,2);

if isempty(tmpDir)
   % get directory from Searchlight 
end



%-------------------------------------------------------------------------%
%% Main part (how can we speedup things?)

%-------------------------------------------------------------------------% 
% Split searchlight definitions?
% local_splitSearchlights(Searchlights,length(ImagesIn{2}),tmpDir);

%-------------------------------------------------------------------------% 
% Read volumes into memory space 
% (this is memory-consuming. Any smarter way?)
for s=1:sn
    Vol{s} = spm_vol(ImagesIn{s});
end

%-------------------------------------------------------------------------% 
% loop over searchlight center nodes (this loop can be parallelized)
nodeID          = Searchlights(1).nodeID;
num_nodes       = Searchlights(1).num_nodes;
Result          = zeros(sn,Ncol,num_nodes)*NaN; % initialize data container
alpha           = 1.0;
CandidateNodes  = find(NsubjNode>=alpha*sn); % define candidate center nodes
Ncandidate      = length(CandidateNodes);
blocksize       = 1; % this needs to be adaptive
Nblock          = ceil(Ncandidate/blocksize);
isDone          = false;
isCalc          = false(Nblock,1);

if verbose==1;fprintf('Blocksize = %d\n',blocksize);end;

for i=1:Nblock;
    from    = blocksize*(i-1)+1;
    to      = min(from+blocksize-1,Ncandidate);
    
    candidates = CandidateNodes(from:to); % node indices
    
    if verbose==1
        fprintf('Running center node %d through %d...\n',candidates(1),candidates(end));
    end
    
    % get data into X
    tic;
    X = cell(sn,1); % initialize X
    for s=1:sn % can be parfor?
        LI = Searchlights(s).LI(candidates);
        if ~isempty(LI)
            % Sample images (time series data) for LI{i_centernode}
            linVox    = unique(cat(2,LI{:})); % unique voxel linear indices
            [I,J,K]   = ind2sub(Vol{s}(1).dim,linVox);
            
            X{s} = sparse(length(Vol{s}),double(max(linVox))); % create sparse matrix for speed and memory?
            for i_image=1:length(Vol{s}) % sample all candidate voxels at once (can be parfor?)
                X{s}(i_image,linVox) = spm_sample_vol(Vol{s}(i_image),double(I),double(J),double(K),0);
            end;
        end
    end;toc;
    
    tic;
    % run pcm for each candidate center node one by one
    for b=1:blocksize
        Nsubj   = NsubjNode(candidates(b));
        Y       = cell(Nsubj,1);
        c       = 1;
        for s=1:sn
            LI = Searchlights(s).LI(candidates(b));
            Y{c} = full(X{s}(:,LI{:}));
            c = c+1;
        end        
        Result(:,:,nodeID(candidates(b))) = feval(pcmfcn,Y,params{:});
    end;toc
    
end
varargout = {};

%-------------------------------------------------------------------------% 
% Write result into .metric file
local_saveMetricFile(Result,metricNames,columnNames);

end

%% Local functions

%-------------------------------------------------------------------------% 
% Split searchlight into optimal size to save memory
%-------------------------------------------------------------------------% 
function local_splitSearchlights(Searchlights,NInputFiles,tmpDir)
sn = length(Searchlights);
if length(NInputFiles)==1
    NInputFiles = repmat(NInputFiles,sn,1);
end

% Calculate memory size (=Ncenternodes*Nvoxel*Ndatapoints*Nsubjects*64bit)
for s=1:sn
    nVox                    = length([Searchlight(s).LI{Searchlights(s).Good}])/...
                                sum(Searchlights(s).Good);
    totalMemory(s)          = numel(Searchlights(s).LI)*nVox*NInputFiles(s);
end

b=1;

% cd(tmpDir);
% while (~isDone)
%     % Find the candidate search lights
%     % This starts in the x-direction and finds the slice of search lights
%     % that are within that slice, and similar for y and z.
%     % then it moves the block tightly to this location
%     isCandidate = ~isCalc;
%     for i=1:3
%         IJKc1(i)    = min(Searchlight.voxmin(isCandidate,i));
%         IJKc2(i)    = IJKc1(i)+blockSize(i);     
%         isCandidate = isCandidate & Searchlight.voxmin(:,i)>=IJKc1(i) & ...
%                                     Searchlight.voxmax(:,i)<=IJKc2(i);
%     end;
%     IJKc1  = min(Searchlight.voxmin(isCandidate,:),[],1);
%     IJKc2  = IJKc1+blockSize;
%     isCandidate=~isCalc & Searchlight.voxmin(:,1)>=IJKc1(1) & Searchlight.voxmin(:,2)>=IJKc1(2) &  ...
%                           Searchlight.voxmin(:,3)>=IJKc1(3) & Searchlight.voxmax(:,1)<=IJKc2(1) &  ...
%                           Searchlight.voxmax(:,2)<=IJKc2(2) & Searchlight.voxmax(:,3)<=IJKc2(3);
%     j=find(isCandidate);
%     if (isempty(j));
%         break;
%     end;
%     
%     % Save the substructure as a tempory file
%     T.LI    = {Searchlight.LI{j}}';
%     T.j     = j;
%     T.voxel = Searchlight.voxel(j,:);
%     if verbose==0; 
%         fprintf('block %d Corner: %d %d %d length:%d  \n',b,IJKc1,length(j));
%     end
%     save(sprintf('temp_%2.2d.mat',b),'-struct','T');
%     isCalc(j) = true;
%     isDone    = all(isCalc);
%     b = b+1;
% end
% numBlocks = b-1;
end

%-------------------------------------------------------------------------%
% Save result as .metric file
%-------------------------------------------------------------------------% 
function local_saveMetricFile(R,metricNames,column_names)
[Nsubj,Ncol,Nnode] = size(R);

%make the struct
M = struct();
for i=1:Nsubj
    M.save_names        = metricNames{i};
    M.column_name       = column_names;
    M.data              = permute(R(i,:,:),[3 2 1]);
    M.num_rows          = Nnode;
    M.num_cols          = Ncol;
    minmax_colormapping = repmat([-1 1],Ncol,1);
    M.column_color_mapping  = minmax_colormapping;
    M.encoding              = {'BINARY'};
    M.index                 = (0:(Nnode-1))';
    
    caret_savemetric(M.save_names, rmfield(M,'save_names') );   %flexible
end

end