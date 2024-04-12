%%% this is the non-axis aligned reachable set computation. Note that in
%%% this instance we assume the initial set (defined as a zonotope
%%% using tX0_center and tX0_gen is not axis aligned and compute
%%% a polytopic outerapproximation to it before doing analysis.
%% Initial problem setup
dim = 3;
T = 0.03;
epsilon = 1e-2; %% this epsilon determines how large of an overapproximation to the original initial set we construct using box aligned objects

% set options for reachability analysis:
time = 0:0.01:T;

% fix the random generator
rng(123)

% Construct orthonormal columns from svd
[orthCol, ~, ~] = svd(randn(dim, dim));

% Construct orthogonally deomposable tensor using definition
odecoTensor = zeros(dim, dim, dim);
lambda = [-0.5, -0.2, 0.1]; % change while changing dimensions
for j = 1:size(orthCol, 2)
    kronVector = kron(kron(orthCol(:, j), orthCol(:, j)), orthCol(:, j));
    odecoTensor = odecoTensor+lambda(j)*reshape(kronVector, dim, dim, dim);
end

% creating initial condition sets
X0_center = [ 10; 20; 30 ]; % change while changing dimensions
X0_gen = 3 * eye( dim );
Zono_X0 = zonotope( X0_center, X0_gen );

tX0_center = orthCol' * X0_center;
tX0_gen = orthCol' * X0_gen;
Zono_tX0 = zonotope( tX0_center, tX0_gen );

% obtaining box decompositions of non axis-aligned set 
outerintervals_tX0 = multi_outer( mptPolytope( Zono_tX0 ), ...
    mptPolytope( Zono_tX0 ), epsilon );

% Explicit formula original system (note the "a" argument is the
% coordinates of the projection of the initial condition onto the orthCol'
% basis set)
homoPolyExp = @(t, a) sum((1 - lambda .* a .* t).^(-1) .* a .* orthCol, 2);

% the dynamical system (xdot)
homoPoly = @(x, u) reshape(reshape(odecoTensor, [numel(x)^2, numel(x)]) * x, numel(x), numel(x))*x;

% Explicit formula in the transformed reference frame, the orthcol terms
% found in homoPolyExp become I here because the axes are aligned to the
% Z-Eigenvectors
thomoPoly = @(a, u) sum((lambda' .* a.^2) .* eye(numel(a)), 2);

% the dynamical system in the transformed frame (xdot) 
thomoPolyExp =@(t,a) sum((1 - lambda' .* a .* t).^(-1) .* a .* eye(numel(a)), 2);

%% Algorithm 2 of the paper
for i = 1:length(time)
    
    for k = 1:size( outerintervals_tX0, 2 )
        % iterating over each box of the decomposition
        zonoInt = zonotope( outerintervals_tX0( :, k ) );
        outerZono_tX0_gen = zonoInt.generators;
        outerZono_tX0_center = zonoInt.center;
        R{ i, k } = 0;
        cen = zeros(3,1);
        gen = [];
  
        for j = 1:size( outerZono_tX0_gen, 2 )
            upp = thomoPolyExp( time( i ), outerZono_tX0_center+ outerZono_tX0_gen( :, j ) );
            low = thomoPolyExp( time( i ), outerZono_tX0_center-outerZono_tX0_gen( :, j ) );
            cen(j) = (upp(j)+low(j))/2;
            gen = [gen,(upp-low)/2];
        end
        R{ i, k } = zonotope( cen, gen );
    end
end

% flipping RS to original coordinate system
S = R;
for i = 1:size(R,1)
    for k = 1:size(R,2)
        c = R{i, k}.center;
        g = R{i, k}.generators;
        S{i,k} = zonotope( orthCol * c, orthCol * g );
    end
end


%% Problem setup for implementing CORA
% figuring out the initial zonotope set
params.R0 = zonotope( tX0_center,  tX0_gen );
params.U = zonotope([0]);
params.tStart = 0;
params.tFinal = T;

sys = nonlinearSys(thomoPoly );

% defining initial conditions
hparams.R0 = zonotope( [ X0_center,  X0_gen ] );
hparams.U = zonotope([0]);
hparams.tStart = 0;
hparams.tFinal = T;

% simulation parameters
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
simRes = simulateRandom(hsys, hparams, simOpt);

% Generating reachable sets using CORA
R_cora = reach(hsys, hparams, options);
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
        curr = mptPolytope(S{i,1});
        for o = 1:size(R,2)
        curr = or(curr,mptPolytope(S{i,o}));
        end
        plot(curr, projDim, 'c', 'Filled', false, 'EdgeColor', 'red','DisplayName', 'RS-my algo');
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

