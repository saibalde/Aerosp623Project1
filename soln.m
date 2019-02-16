% soln.m
% 
% A sequence of commands to generate the results for the entire project

% generate initial mesh
bump0gen;

% compute I2E, B2E, In, Bn, Area for test mesh
[Nodes1, B2N1, Title1, E2N1, Basis1, Order1] = gridread('test.gri');
[I2E1, B2E1, In1, Bn1, Area1] = gridmats(Nodes1, B2N1, E2N1);

% compute I2E, B2E, In, Bn, Area for bump0 mesh
[Nodes2, B2N2, Title2, E2N2, Basis2, Order2] = gridread('bump0.gri');
[I2E2, B2E2, In2, Bn2, Area2] = gridmats(Nodes2, B2N2, E2N2);

% verify test.gri and bump0.gri
VerVal1 = gridverify('test.gri');
VerVal2 = gridverify('bump0.gri');

% uniform refinements of bump0.gri
gridrefine('bump0.gri', 'bump1.gri');
gridrefine('bump1.gri', 'bump2.gri');
gridrefine('bump2.gri', 'bump3.gri');
gridrefine('bump3.gri', 'bump4.gri');

% plot refinements
plotgri('bump0.gri');
plotgri('bump1.gri');
plotgri('bump2.gri');
% plotgri('bump3.gri'); % commented out as takes fairly long
% plotgri('bump4.gri'); % commented out as takes fairly long