function [I2E, B2E, In, Bn, Area] = gridmats(Nodes, NB, NE)
%GRIDMATS Compute matrices corresponding to the mesh
%
% INPUT : Nodes   = [nnodes x dim] array; coordinates of the nodes
%         NB      = [nbgroups x 1] cell array of boundary groups, where
%                   each cell contains an array mapping boundary faces to
%                   node indices
%         NE      = [negroups x 1] cell array of element groups, where
%                   each cell nontains an array mapping elements to node
%                   indices
% OUTPUT: I2E     = [niface x 4] array; map from interior faces to
%                   corresponding elements; each row contains four
%                   entries: [elemL, faceL, elemR, faceR], defined as
%                     - elemL = number of element on the left
%                     - faceL = local face number w.r.t. left element
%                     - elemR = number of element on the right
%                     - faceR = local face number w.r.t. right element
%         B2E     = [nbface x 3] array; map from boundary faces to
%                   correspondign elemetns; each row contains three
%                   entries: [elem, face, bgrp], defined as
%                    - elem = number of adjoining element
%                    - face = local face number w.r.t. adjoining element
%                    - bgrp = boundary group number
%         In      = [niface x 2] array; map of interior faces to normals,
%                   pointing from the left to right face
%         Bn      = [nbface x 2] array; map of boundary faces to outward
%                   normals
%         Area    = [nelems x 1] array; list of element areas

%% Grid parameters

[nNode, ~] = size(Nodes);

nBGroup = length(NB);
nBFace = zeros(nBGroup, 1);
nf = zeros(nBGroup, 1);
for iBFace = 1 : nBGroup
    [nBFace(iBFace), nf(iBFace)] = size(NB{iBFace});
end
nBFaceTot = sum(nBFace);

NE = cell2mat(NE);                              % flatten element groups
nElemTot = size(NE, 1);

%% Hash matrices

E = sparse(nNode, nNode);                       % store element number
F = sparse(nNode, nNode);                       % store local face number

%% Allocate memory

I2E  = zeros(ceil(3 * nElemTot / 2), 4);        % overallocate
In   = zeros(ceil(3 * nElemTot / 2), 2);        % overallocate
B2E  = zeros(nBFaceTot, 3);
Bn   = zeros(nBFaceTot, 2);
Area = zeros(nElemTot, 1);

%% Update I2E and In

% initialize number of interior faces to zero
nIFace = 0;

% loop over elements
for iElemTot = 1 : nElemTot
    % node numbers for elements
    ne = NE(iElemTot, :);

    % loop over edges
    for edge = 1 : 3
        % compute local node numbers
        ne1loc = mod(edge - 1, 3) + 1;          % node on the edge
        ne2loc = mod(edge,     3) + 1;          % node on the edge
        ne3loc = 6 - ne1loc - ne2loc;           % node opposite to edge

        % compute global node numbers
        ne1 = ne(ne1loc);                       % node on the edge
        ne2 = ne(ne2loc);                       % node on the edge
        ne3 = ne(ne3loc);                       % node opposite to edge

        % sort the global node numbers, lying on the edge
        nemin = min(ne1, ne2);
        nemax = max(ne1, ne2);

        % hit edge for first time; could be interior or boundary edge
        if E(nemin, nemax) == 0
            % store the element number
            E(nemin, nemax) = iElemTot;

            % store local index of opposite node
            F(nemin, nemax) = ne3loc;

        % hit edge for second time, definitly interior edge
        elseif E(nemin, nemax) > 0
            % increment number of interior edges
            nIFace = nIFace + 1;

            % recover element number from previous visit
            oldelem = E(nemin, nemax);

            % update I2E: we know oldelem < elem
            elemL = oldelem;
            elemR = iElemTot;
            faceL = F(nemin, nemax);
            faceR = ne3loc;
            I2E(nIFace, :) = [elemL faceL elemR faceR];

            % compute normal pointing into current element
            e12 = Nodes(ne2, :) - Nodes(ne1, :);
            e13 = Nodes(ne3, :) - Nodes(ne1, :);
            In(nIFace, :) = unitdir(e12, e13);

            % make element number negative
            E(nemin, nemax) = -1;
        else
            % negative element number, something is wrong!
            error('Mesh input error');
        end
    end
end
    
% clip I2E and In
I2E = I2E(1:nIFace, :);
In  = In (1:nIFace, :);

%% Update B2E and Bn

% running total of number of boundary edges already visited
edge = 0;

% loop over boundary face groups
for iBGroup = 1 : nBGroup
    % loop over boundary faces in group
    for iBFace = 1 : nBFace(iBGroup)
        % increament number of boundary edges already visited
        edge = edge + 1;

        % number of nodes in face
        ne1 = NB{iBGroup}(iBFace, 1);
        ne2 = NB{iBGroup}(iBFace, 2);
        
        % sort node indices
        nemin = min(ne1, ne2);
        nemax = max(ne1, ne2);
        
        % update B2E
        elem = E(nemin, nemax);
        face = F(nemin, nemax);
        B2E(edge, :) = [elem face iBGroup];
        
        % find global number of node opposite to edge
        neopp = NE(elem, face);
        
        % compute normal
        a = Nodes(nemax, :) - Nodes(nemin, :);
        b = Nodes(neopp, :) - Nodes(nemin, :);
        Bn(edge, :) = - unitdir(a, b);
    end
end

%% Update Area

% loop over element groups
for iElemTot = 1 : nElemTot
    % get coordinates of the vertices
    A = Nodes(NE(iElemTot, 1), :);
    B = Nodes(NE(iElemTot, 2), :);
    C = Nodes(NE(iElemTot, 3), :);

    % compute area
    Area(iElemTot) = 0.5 * abs(A(1) * (B(2) - C(2)) + ...
                               B(1) * (C(2) - A(2)) + ...
                               C(1) * (A(2) - B(2)));
end

end

function n = unitdir(a, b)
%NORMALDIR Compute unit vector orthogonal to a that makes an acute angle
% with b

n = b - dot(a, b) * a / norm(a)^2;
n = n / norm(n);

end