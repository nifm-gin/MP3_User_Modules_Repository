function [Yestim, score, Xestim] = EstimateParametersFromGrid(Xobs, Xgrid, Ygrid, verb)

narginchk(3, 4)
if nargin == 3, verb = 0; end

% normalization
% Xobs_normalized = Xobs ./ vecnorm(Xobs,2,2);
% Xgrid_normalized = Xgrid ./ vecnorm(Xgrid,2,2);

%% dot-product/scalar product comparison
[score, idx] = max(round(abs(Xobs * Xgrid'), 3), [], 2);

%% r2 coefficient
% score = zeros(size(Xobs_normalized, 1), 1);
% idx = ones(size(Xobs_normalized, 1), 1);
% 
% for i=1:size(Xobs_normalized, 1)
%     sig = Xobs_normalized(i, :);
% %     if ~any(isnan(sig))
% %         pause
% %     end
%     R2 = 1 - sum((sig - Xgrid_normalized).^2, 2) ./ sum((sig - mean(Xgrid_normalized, 2)).^2, 2);
%     [localScore, localIdx] = max(R2);
%     if localScore > score(i)
%         score(i) = localScore;
%         idx(i) = localIdx;
%     end
% end

%% min dist
% Dvox = pdist2(Xobs_normalized,Xgrid_normalized);%Euclidian distance
% [score,idx] = min(Dvox,[],2);


Yestim  = Ygrid(idx,:);
Xestim  = Xgrid(idx,:);