function [Yestim, score, Xestim] = EstimateParametersOnFly(Xobs, Xgrid, Ygrid, geo, Ylabels)

% MRF matching that load a Bloch dico Xgrid, dynamically convolve it by realistic
% vessels 3D distribution (input : geo), and perform matching
fn = fieldnames(geo);
goat = zeros(size(Xobs, 1), 1);

dfindex = find(strcmp('df', Ylabels), 1);
df_values = Ygrid(:, dfindex);
Yestim = zeros(size(Xobs, 1), numel(Ylabels));
Xestim = zeros(size(Xobs));

SO2index = find(strcmp('SO2', Ylabels), 1);
Vfindex = find(strcmp('Vf', Ylabels), 1);
Rindex = find(strcmp('R', Ylabels), 1);

nb_distrib = 1000;
ks = round(linspace(1, numel(fn), nb_distrib));
ks = unique(ks);

for w=1:numel(ks)
    k = ks(w);
%     fprintf('distrib = %d \n', k);
    % convolve Xgrid by geo.(fn{k}).Histo.bin
    Xgrid_tmp = Xgrid;
    matrix = zeros(size(unique(df_values),1), size(geo.(fn{k}).histo.bin, 2));
    
    for a=1:size(unique(df_values), 1)
        if df_values(a) >= 0
            matrix(a,:) = circshift(geo.(fn{k}).histo.bin, df_values(a));
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
    [best, idx]    = max(round(abs(Xobs * Xgrid_tmp'), 10), [], 2);
    
    % save idx if better match than previous convolved dict
    compare = best>goat;
    goat(compare, :) = best(compare, :);
    
    
    Yestim(compare, :)  = Ygrid(idx(compare), :);
    Xestim(compare, :)  = Xgrid_tmp(idx(compare), :);
    Yestim(compare, SO2index) = geo.(fn{k}).SO2;
    Yestim(compare, Vfindex) = geo.(fn{k}).VF;
    Yestim(compare, Rindex) = geo.(fn{k}).R;
end
    

score = goat;