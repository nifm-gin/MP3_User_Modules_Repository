function [Yestim, score, Xestim] = EstimateParametersOnFlyPar_dfmask(Xobs, Xgrid, Ygrid, geo, Ylabels)

% MRF matching that load a Bloch dico Xgrid, dynamically convolve it by realistic
% vessels 3D distribution (input : geo), and perform matching
fn = fieldnames(geo);
goat = zeros(size(Xobs, 1), 1);

dfindex = find(strcmp('df', Ylabels), 1);
df_values = Ygrid(:, dfindex);
Yestim = zeros(size(Xobs, 1), numel(Ylabels));
Xestim = zeros(size(Xobs));

score = zeros(size(Xobs, 1), 1);

SO2index = find(strcmp('SO2', Ylabels), 1);
Vfindex = find(strcmp('Vf', Ylabels), 1);
Rindex = find(strcmp('R', Ylabels), 1);

nbDistrib = length(fn);
distribIdx = round(linspace(1, numel(fn), nbDistrib));
distribIdx = unique(distribIdx);
% distribAll = unique(round(linspace(1, 28006, 28006)));
% distribIdx = setdiff(distribAll, distribIdx);

best = cell(1,nbDistrib);
best(1:nbDistrib) = {NaN(size(Xobs, 1),1)};
idx(:) = best; % do not repeat precedent line to save 2 seconds yay

ps = parallel.Settings;
ps.SchedulerComponents(1).NumWorkers = 32;
ps.SchedulerComponents(1).NumThreads = 4;

if(isempty(gcp))
    parpool(ps.SchedulerComponents(2).NumWorkers);
end

fwhm = 30; %Hz
% mask = (min(df_values) + fwhm < df_values) & (df_values < max(df_values) - fwhm); %constraint fully defined df centered distributions (assume max distribution is 20Hz width)
mask = (0 <= df_values) & (df_values <= 30); 
% mask = (Ygrid(:, 3) >= -50) & (Ygrid(:, 3) <= 50); % df values that are between -50 and 50Hz (column 3 of Ygrid)
% mask = (Ygrid(:,3) == 1);

parfor w=1:numel(distribIdx)
%     fprintf('%' \n', w);
    k = distribIdx(w);
    Xgrid_tmp = Xgrid;

    matrix = zeros(size(unique(df_values),1), size(geo.(fn{k}).histo.bin, 2));
    
    for a=1:size(unique(df_values), 1)
        if df_values(a) >= 0
            matrix(a,:) = circshift(geo.(fn{k}).histo.bin, round(df_values(a)));
        else
            matrix(a,:) = circshift(geo.(fn{k}).histo.bin, round(numel(geo.(fn{k}).histo.bin) - abs(df_values(a))));
        end
    end
    
    Xgrid_tmp = reshape(Xgrid_tmp, size(unique(df_values), 1), [], size(Xgrid, 2));
    
    for t=1:size(Xgrid_tmp, 2)
        Xgrid_tmp(:, t, :) = matrix * squeeze(Xgrid_tmp(:, t, :));
    end
    %
    Xgrid_tmp = reshape(Xgrid_tmp, [], size(Xgrid, 2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xgrid_tmp = abs(Xgrid_tmp);
    Xgrid_tmp = Xgrid_tmp ./ vecnorm(Xgrid_tmp,2,2);
    
    % match
    Xgrid_tmp = Xgrid_tmp(mask,:);
    [best{w}, idx{w}]    = max(round(abs(Xobs * Xgrid_tmp'), 10), [], 2);
end

B = cell2mat(best);
I = cell2mat(idx);
[goat, idxGoat] = max(B, [], 2, 'linear');
iNotNan = ~isnan(goat);
[rI, cI] = ind2sub(size(I), idxGoat(iNotNan));

Ygrid = Ygrid(mask,:);
Yestim(rI, :) = Ygrid(I(idxGoat(iNotNan)), :);

for j = 1:numel(rI)
    Yestim(rI(j), SO2index) = geo.(fn{distribIdx(cI(j))}).SO2;
    Yestim(rI(j), Vfindex) = geo.(fn{distribIdx(cI(j))}).VF;
    Yestim(rI(j), Rindex) = geo.(fn{distribIdx(cI(j))}).R;
end

score(rI) = goat(rI);

end