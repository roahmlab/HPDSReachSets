%%% this is the axis aligned reachable set computation for systems 
%%% with constant input. Note that in this instance, we assume the initial  
%%% set (defined as a zonotope) using X0_center and X0_gen to be 
%%% axis-aligned to the Z-Eigenvectors

close all;
clear;

%% Initial setup of problem
dim = 3;
T = 0.05;

% set options for reachability analysis:
time = 0:0.01:T;

% Construct orthonormal columns from svd
[orthCol, ~, ~] = svd(randn(dim, dim));

% Construct orthogonally decomposable tensor
odecoTensor = zeros(dim, dim, dim);
lambda = [-0.5, -0.2, 0.1]; % edit according to dimension of system
for j = 1:size(orthCol, 2)
    kronVector = kron(kron(orthCol(:, j), orthCol(:, j)), orthCol(:, j));
    odecoTensor = odecoTensor+lambda(j)*reshape(kronVector, dim, dim, dim);
end

%creating random constant input
b = randn(dim,1);

%Calculating btilde - btilde will be equal to b in transformed coordinate
%system
btilde = b;

% creating initial condition sets -  ensure that 0 does not belong to set
% in any of the coordinates, these sets are axis-aligned
tX0_center = [ 10; 40; 50 ]; % edit according to dimension of system
tX0_gen = 6 * eye( dim );

% the initial sets in the original frame of reference
X0_center = orthCol * tX0_center;
X0_gen = orthCol * tX0_gen;

% tVec represents eigenvectors in transformed coordinate system- i.e.,
% standard basis
tVec = eye(dim);

%% Algorithm 3 of paper
for i = 1:length(time)
    cen = zeros(dim,1);
    gen = [];
    R{ i } = 0;
    for m = 1:size( tX0_gen, 2 )
        upp=solveHypergeom(3,btilde(m),lambda(m),tX0_center(m)+tX0_gen(m, m),time(i));
        low = solveHypergeom(3,btilde(m),lambda(m),tX0_center(m)-tX0_gen(m, m),time(i));
        cen(m) = (upp+low)/2;
        gen = [gen,(upp-low)*tVec(:,m)/2];
    end
    R{i} = zonotope( cen, gen );
end

% flipping to different coordinate system
S = R;
for i = 1:length(R)
    c = R{i}.center;
    g = R{i}.generators;
    S{i} = zonotope( orthCol * c, orthCol * g );
    
end

%% Setup for CORA
% the dynamical system xdot = Ax^(k-1)
homoPoly =@(x,u) [...
    [1 0 0]*(reshape(reshape(odecoTensor, [9, 3])*[x(1); x(2); x(3)], 3, 3)*[x(1); x(2); x(3)]+b);
    [0 1 0]*(reshape(reshape(odecoTensor, [9, 3])*[x(1); x(2); x(3)], 3, 3)*[x(1); x(2); x(3)]+b);
    [0 0 1]*(reshape(reshape(odecoTensor, [9, 3])*[x(1); x(2); x(3)], 3, 3)*[x(1); x(2); x(3)]+b);
    ];  % edit according to dimension of system

tank = nonlinearSys('tank',homoPoly);

% simulating examples - proxy for the ground truth
hparams.R0 = zonotope( [ X0_center,  X0_gen ] );
hparams.U = zonotope([0]);
hparams.tStart = 0;
hparams.tFinal = T;

options.timeStep = 0.01;
options.taylorTerms=5; % number of taylor terms for reachable sets
options.zonotopeOrder= 10; % zonotope order... increase this for more complicated systems.
options.maxError = 1000*ones(dim, 1); % our zonotopes shouldn't be "splitting", so this term doesn't matter for now
options.verbose = false;
options.uTrans = 0; % we won't be using any inputs, as traj. params specify trajectories
options.advancedLinErrorComp = 0;
options.tensorOrder = 3;
options.reductionInterval = inf;
options.reductionTechnique = 'girard';
options.alg = 'poly';
options.intermediateOrder = 50;
options.errorOrder = 1;

simOpt.points = 100;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;
hsys = nonlinearSys(homoPoly);
% simulating the randomly initialized points
simRes = simulateRandom(hsys, hparams, simOpt);

% Generating reachable sets using CORA 
R_cora = reach(tank, hparams, options);

%% Plotting the results - 3 dimensional systems; modify as needed
dims = {[1 2] [2 3] [3 1]};
dim_labels = {'x', 'y', 'z'}; % Labels for dimensions

for k = 1:length(dims)
    projDim = dims{k};
    dim_label_x = dim_labels{projDim(1)};
    dim_label_y = dim_labels{projDim(2)};
    
    figure; hold on;

    % plot reachable sets
    for i = 1:length(time)
        plot(S{i}, projDim, 'c', 'Filled', false, 'EdgeColor', 'red','DisplayName', 'Proposed Algorithm');
    end

    % plot initial set
    plot(hparams.R0, projDim, 'w', 'Filled', false, 'EdgeColor', 'blue', 'DisplayName', 'Initial set');

    % plot simulation results
    plot(simRes, projDim, 'k', 'DisplayName', 'Simulation');
    plot(R_cora, projDim, 'Filled', false,'EdgeColor', 'green', 'DisplayName', 'RS-Cora');

    % Label plot
    xlabel(dim_label_x);
    ylabel(dim_label_y);
    title(['Projection onto ' dim_label_x ' vs ' dim_label_y]);
    %legend();
end
