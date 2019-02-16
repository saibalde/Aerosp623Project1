function gridwrite(Nodes, NB, Title, NE, Basis, Order, GriFile)
%GRIDWRITE Write the grid file
%
% INPUT : Nodes   = coordinates of the nodes
%         NB      = node indices corresponding to boundary faces
%         Title   = names of the boundary face groups
%         NE      = node indices corresponding to elements
%         Basis   = type of the elements
%         Order   = order of the elements
%         GriFile = name of the .gri files

[nNode, Dim] = size(Nodes);
nBGroup = length(NB);
nElemGroup = length(NE);

nElem = zeros(nElemGroup, 1);
nn = zeros(nElemGroup, 1);
for iElemGroup = 1 : nElemGroup
    [nElem(iElemGroup), nn(iElemGroup)] = size(NE{iElemGroup});
end
nElemTot = sum(nElem);

% open file
fileId = fopen(GriFile, 'w');

% write number of nodes, number of elements, mesh dimension
fprintf(fileId, "%d %d %d\n", nNode, nElemTot, Dim);

% write node coordinates
for iNode = 1 : nNode
    for i = 1 : Dim
        fprintf(fileId, "% 12.5e ", Nodes(iNode, i));
    end
    fprintf(fileId, "\n");
end

% write number of boundary groups
fprintf(fileId, "%d\n", nBGroup);

% write boundary groups
for iBGroup = 1 : nBGroup
    % write boundary group metadata
    [nBFace, nf] = size(NB{iBGroup});
    fprintf(fileId, "%d %d %s\n", nBFace, nf, Title{iBGroup});
    
    % write boundary group details
    for iBFace = 1 : nBFace
        for i = 1 : nf
            fprintf(fileId, "%d ", NB{iBGroup}(iBFace, i));
        end
        fprintf(fileId, "\n");
    end
end

% write element groups
for iElemGroup = 1 : nElemGroup
    % write element group metadata
    fprintf(fileId, "%d %d %s\n", nElem(iElemGroup), Order(iElemGroup), ...
        Basis{iElemGroup});
    
    % write element group details
    for iElem = 1 : nElem(iElemGroup)
        for i = 1 : nn(iElemGroup)
            fprintf(fileId, "%d ", NE{iElemGroup}(iElem, i));
        end
        fprintf(fileId, "\n");
    end
end

% close file
fclose(fileId);

end

