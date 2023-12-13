function [Yestim, score, Xestim] = EstimateParametersFromGrid_Chunk(Xobs, Xgrid, Ygrid, verb)

% MRF matching with norm options and chunking of the dictionary to
% avoid RAM overload

if nargin == 3, verb = 0; end

if size(Xgrid,1) > 300000
    
    N = size(Xgrid,1);
    chunk = 100e3;
    I = floor(N/chunk);
    best = zeros(size(Xobs,1),1);
    idx = ones(size(Xobs,1),1);
    
    for i = 1:I
        [localBest, localIdx] = max(round(abs(Xobs * Xgrid((i-1)*chunk +1 : i*chunk,:)'), 10), [], 2);
        idx(best<localBest) = localIdx(best<localBest)+(i-1)*chunk;
        best(best<localBest)=localBest(best<localBest);
        
    end
    if mod(N,chunk)~=0
        [localBest, localIdx] = max(round(abs(Xobs * Xgrid(i*chunk +1 : end,:)'), 10), [], 2);
        idx(best<localBest) = localIdx(best<localBest)+(I)*chunk;
        best(best<localBest)=localBest(best<localBest);
        
    end
else
    [best, idx]    = max(round(abs(Xobs * Xgrid'), 10), [], 2);
end

score = best;
Yestim  = Ygrid(idx, :);
Xestim  = Xgrid(idx, :);