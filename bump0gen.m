function bump0gen()

%% runtime path modification

addpath('mesh2d');
libpath;

%% specify boundary nodes

% To ensure the mesh is correctly written into bump0.gri file, choose the
% initial boundary vertices to be dense enough to ensure no extra boundary
% edges are generated. Otherwise the numbering of vertices gets messed
% up, which leads to problems in specifying the boundary faces

% Also, set the bottom face y-coordinates to zero for now

boundary_nodes = [-1.50 0.00; % bottom left corner
                  -1.25 0.00;
                  -1.00 0.00;
                  -0.75 0.00;
                  -0.50 0.00;
                  -0.25 0.00;
                  -0.15 0.00;
                  -0.05 0.00;
                   0.00 0.00;
                   0.05 0.00;
                   0.15 0.00;
                   0.25 0.00;
                   0.50 0.00;
                   0.75 0.00;
                   1.00 0.00;
                   1.25 0.00;
                   1.50 0.00; % bottom right corner
                   1.50 0.20;
                   1.50 0.40;
                   1.50 0.60;
                   1.50 0.80; % top right corner
                   1.25 0.80;
                   1.00 0.80;
                   0.75 0.80;
                   0.50 0.80;
                   0.25 0.80;
                   0.00 0.80;
                  -0.25 0.80;
                  -0.50 0.80;
                  -0.75 0.80;
                  -1.00 0.80;
                  -1.25 0.80;
                  -1.50 0.80; % top left corner
                  -1.50 0.60;
                  -1.50 0.40;
                  -1.50 0.20]; % one above bottom left corner

% fix the y-coordinates of the bottom face: find the indices for which the
% second column entries are zero, then fix these entries using the bump
% formula

inds = find(boundary_nodes(:, 2) == 0);
boundary_nodes(inds, 2) = bump(boundary_nodes(inds, 1));

%% generate mesh using mesh2d

opts.kind = 'delfront';
opts.rho2 = 1.025;
hfun      = 0.25;

[vert, ~, tria, ~] = refine2(boundary_nodes, [], [], opts, hfun);

%% separate out the bottom, right, top and left boundary faces

nBot = 16;
nRht = 4;
nTop = 12;
nLft = 4;

n1 = nBot;
n2 = nBot + nRht;
n3 = nBot + nRht + nTop;
n4 = nBot + nRht + nTop + nLft;

bot = zeros(nBot, 2); bot(:, 1) =      1 : n1; bot(:, 2) = bot(:, 1) + 1;
rht = zeros(nRht, 2); rht(:, 1) = n1 + 1 : n2; rht(:, 2) = rht(:, 1) + 1;
top = zeros(nTop, 2); top(:, 1) = n2 + 1 : n3; top(:, 2) = top(:, 1) + 1;
lft = zeros(nLft, 2); lft(:, 1) = n3 + 1 : n4; lft(:, 2) = lft(:, 1) + 1; lft(nLft, 2) = 1;

%% mesh details, in terms gridwrite.m understands

Nodes = vert;
NB    = {bot, rht, top, lft};
Title = {'Bottom', 'Right', 'Top', 'Left'};
NE    = {tria};
Basis = {'TriLagrange'};
Order = [1];

%% write the grid

gridwrite(Nodes, NB, Title, NE, Basis, Order, 'bump0.gri');

end
