function [Nodes, NB, Title, NE, Basis, Order] = gridread(GriFile)
%GRIDREAD Read the grid file
%
% INPUT : GriFile = name of the .gri file
% OUTPUT: Nodes   = coordinates of the nodes
%         NB      = node indices corresponding to boundary faces
%         Title   = names of the boundary face groups
%         NE      = node indices corresponding to elements
%         Basis   = type of the elements
%         Order   = order of the elements

% open file
fileId = fopen(GriFile, 'r');

% read in number of nodes, elements and dimension
A = textscan(fileId, '%d', 3);
nNode    = A{1}(1);
nElemTot = A{1}(2);
Dim      = A{1}(3);

% verify we are using 2D mesh
if Dim ~= 2
    error('We only support 2D meshes');
end

% read the node coordinates
Nodes = zeros(nNode, 2);
for iNode = 1 : nNode
    A = textscan(fileId, '%f', 2);
    Nodes(iNode, :) = A{1};
end

% read in number of boundary groups
A = textscan(fileId, '%d', 1);
nBGroup = A{1};

nBFace = zeros(nBGroup, 1);
nf     = zeros(nBGroup, 1);
Title  = cell(nBGroup, 1);
NB     = cell(nBGroup, 1);

% read the boundary groups
for iBGroup = 1 : nBGroup
    % read number of faces in the group, nodes per face and group title
    A = textscan(fileId, '%d %d %s', 1);
    nBFace(iBGroup) = A{1};
    nf(iBGroup) = A{2};
    Title{iBGroup} = A{3}{1};
    
    % verify boundary faces are composed of two nodes each
    if nf(iBGroup) ~= 2
        error('We only support boundary faces with 2 nodes')
    end
    
    % read the indices for nodes on the boundary faces
    NB{iBGroup} = zeros(nBFace(iBGroup), 2);
    for iBFace = 1 : nBFace(iBGroup)
        A = textscan(fileId, '%d', 2);
        NB{iBGroup}(iBFace, :) = A{1};
    end
end

% read the element groups
nElemGroup = 0;
nElemCur = 0;
while nElemCur < nElemTot
    nElemGroup = nElemGroup + 1;

    % read number of elements in group, order and basis of element
    A = textscan(fileId, '%d %d %s', 1);
    nElem(nElemGroup) = A{1};
    Order(nElemGroup) = A{2};
    Basis{nElemGroup} = A{3}{1};
    
    % ensure linear, triangular element types
    if ~strcmp(Basis{nElemGroup}, 'TriLagrange')
        error('We only support TriLagrange element type');
    end
    
    if Order(nElemGroup) ~= 1
        error('We only support first order elements');
    end

    % read in the elements in the group
    NE{nElemGroup} = zeros(nElem(nElemGroup), 3);
    for iElem = 1 : nElem(nElemGroup)
        A = textscan(fileId, '%d', 3);
        NE{nElemGroup}(iElem, :) = A{1};
    end
    
    % update number of elements read in
    nElemCur = nElemCur + nElem(nElemGroup);
end

% ensure that the number of elements is correct
if (nElemTot ~= nElemCur)
    error('Number of elements does not match!');
end

% close file
fclose(fileId);

end

