function gridrefine(GriFileIn, GriFileOut)
%GRIDREFINE Perform an uniform refinement of mesh
%
% INPUT : GriFileIn  = Name of input grid file
%         GriFileOut = Name of output grid file

%% read input grid

[Nodes, NB, Title, NE, Basis, Order] = gridread(GriFileIn);

%% grid parameters

nNode      = size(Nodes, 1);
nBGroup    = length(NB);
nElemGroup = length(NE);

%% hash matrix

% N(n1, n2) is the index of the newly created node between nodes n1 and n2
N = sparse(nNode, nNode);

%% running number of how many nodes in the mesh

nNodeCur = nNode;

%% refine boundary faces

% storage for new boundary groups
NB2 = cell(nBGroup, 1);

% loop over boundary groups
for iBGroup = 1 : nBGroup
    % number of faces in boundary group
    nBFace = size(NB{iBGroup}, 1);
    
    % allocate memory for twice as many faces
    NB2{iBGroup} = zeros(2 * nBFace, 2);
    
    % loop over existing faces
    for iBFace = 1 : nBFace
        % node indices corresponding to boundary edge
        n1= NB{iBGroup}(iBFace, 1);
        n2= NB{iBGroup}(iBFace, 2);
        
        % sort the node indices
        nmin = min(n1, n2);
        nmax = max(n1, n2);
        
        % create new node
        nNodeCur = nNodeCur + 1;
        pmin = Nodes(nmin, :);
        pmax = Nodes(nmax, :);
        p = 0.5 * (pmin + pmax);
        
        % conform to bump geometry
        if iBGroup == 1
            p(2) = bump(p(1));
        end
        
        % update hash matrix with new node index
        N(nmin, nmax) = nNodeCur;

        % update Nodes matrix
        Nodes(nNodeCur, :) = p;
        
        % insert the boundary edges to new boundary group
        NB2{iBGroup}(2 * iBFace - 1, :) = [nmin, nNodeCur];
        NB2{iBGroup}(2 * iBFace,     :) = [nNodeCur, nmax];
    end
end

%% refine elements

% storage for new element groups
NE2 = cell(nElemGroup, 1);

% loop over element groups
for iElemGroup = 1 : nElemGroup
    % number of elements in group
    nElem = size(NE{iElemGroup}, 1);
    
    % memory for new element group
    NE2{iElemGroup} = zeros(4 * nElem, 3);
    
    % loop over existing elements
    for iElem = 1 : nElem
        % node indices corresponding to element
        ne = NE{iElemGroup}(iElem, :);
        
        % storage for new nodes
        nw = zeros(3, 1);
        
        % loop over edges
        for edge = 1 : 3
            % local node numbers corresponding to edge
            n1loc = edge;
            n2loc = mod(edge, 3) + 1;
            
            % global node number corresponding to edge
            n1 = ne(n1loc);
            n2 = ne(n2loc);
            
            % sort global node numbers
            nmin = min(n1, n2);
            nmax = max(n1, n2);
            
            % check against if new node needs to be created
            if N(nmin, nmax) == 0
                % create new middle node
                nNodeCur = nNodeCur + 1;
                p1 = Nodes(n1, :);
                p2 = Nodes(n2, :);
                p  = 0.5 * (p1 + p2);
                
                % update hash matrix
                N(nmin, nmax) = nNodeCur;
                
                % update nodes matrix
                Nodes(nNodeCur, :) = p;
            end
                
            % copy new node index over to local new node list nw
            nw(edge) = N(nmin, nmax);
        end
        
        % insert the elements to new element group
        NE2{iElemGroup}(4 * iElem - 3, :) = [ne(1), nw(1), nw(3)];
        NE2{iElemGroup}(4 * iElem - 2, :) = [ne(2), nw(2), nw(1)];
        NE2{iElemGroup}(4 * iElem - 1, :) = [ne(3), nw(3), nw(2)];
        NE2{iElemGroup}(4 * iElem,     :) = [nw(1), nw(2), nw(3)];
    end
end

% write output grid
gridwrite(Nodes, NB2, Title, NE2, Basis, Order, GriFileOut);

end

