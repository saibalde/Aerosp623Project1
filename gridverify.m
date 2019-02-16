function [maxval, normalSum] = gridverify(GriFile)
%GRIDVERIFY Verify a grid file
%
% INPUT : GriFile   = name of the grid file
% OUTPUT: normalSum = [nelems x dim] array; each row contains an edge
%                     length weighted sum of outward normal directions
%                     for each elements
%         maxval    = maximum value of 2-norms of each row of normalSum

%% read grid file, and generate grid matrices

[Nodes, NB, ~, NE, ~, ~] = gridread(GriFile);
[I2E, B2E, In, Bn, ~] = gridmats(Nodes, NB, NE);

NE = cell2mat(NE);                              % flatten element groups

%% grid parameters

nElem  = size(NE,  1);
nIEdge = size(I2E, 1);
nBEdge = size(B2E, 1);

%% allocate memory for normalSum

normalSum = zeros(nElem, 2);

%% loop over interior edges

for iEdge = 1 : nIEdge
    % left face
    elem = I2E(iEdge, 1);
    vert = NE(elem, :);
    if I2E(iEdge, 2) == 1
        len = norm(Nodes(vert(2), :) - Nodes(vert(3), :));
    elseif I2E(iEdge, 2) == 2
        len = norm(Nodes(vert(3), :) - Nodes(vert(1), :));
    elseif I2E(iEdge, 2) == 3
        len = norm(Nodes(vert(1), :) - Nodes(vert(2), :));
    else
        error('Something wrong with I2E');
    end
    normalSum(elem, :) = normalSum(elem, :) + len * In(iEdge, :);
    
    % right face
    elem = I2E(iEdge, 3);
    vert = NE(elem, :);
    if I2E(iEdge, 4) == 1
        len = norm(Nodes(vert(2), :) - Nodes(vert(3), :));
    elseif I2E(iEdge, 4) == 2
        len = norm(Nodes(vert(3), :) - Nodes(vert(1), :));
    elseif I2E(iEdge, 4) == 3
        len = norm(Nodes(vert(1), :) - Nodes(vert(2), :));
    else
        error('Something wrong with I2E');
    end
    normalSum(elem, :) = normalSum(elem, :) - len * In(iEdge, :);
end

%% loop over boundary edges

for bEdge = 1 : nBEdge
    elem = B2E(bEdge, 1);
    vert = NE(elem, :);
    if B2E(bEdge, 2) == 1
        len = norm(Nodes(vert(2), :) - Nodes(vert(3), :));
    elseif B2E(bEdge, 2) == 2
        len = norm(Nodes(vert(3), :) - Nodes(vert(1), :));
    elseif B2E(bEdge, 2) == 3
        len = norm(Nodes(vert(1), :) - Nodes(vert(2), :));
    else
        error('Something wrong with B2E');
    end
    normalSum(elem, :) = normalSum(elem, :) + len * Bn(bEdge, :);
end

%% maximum norm of rows of normalSum

maxval = 0.0;
for elem = 1 : nElem
    maxval = max(maxval, norm(normalSum(elem, :)));
end

end

