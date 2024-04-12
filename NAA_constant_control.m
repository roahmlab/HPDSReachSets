%%% this is the non-axis aligned reachable set computation for 
%%% systems with constant input. Note that in this instance, 
%%% we assume the initial set (defined as a zonotope) using tX0_center and 
%%% tX0_gen is not axis aligned and compute a polytopic outerapproximation 
%%% to it before doing analysis.

close all;
clear;

%% Initial setup of problem
dim = 3;
T = 0.05;
epsilon = 1;

% set options for reachability analysis:
time = 0:0.01:T;

% Construct orthonormal columns from svd
[orthCol, ~, ~] = svd(randn(dim, dim));

% Construct orthogonally decomposable tensor
odecoTensor = zeros(dim, dim, dim);
lambda = [-0.5, -0.2, 0.1];
for j = 1:size(orthCol, 2)
    kronVector = kron(kron(orthCol(:, j), orthCol(:, j)), orthCol(:, j));
    odecoTensor = odecoTensor+lambda(j)*reshape(kronVector, 3, 3, 3);
end

%creating random constant input
b = randn(3,1);

%Calculating btilde - btilde will be equal to b in transformed coordinate
%system
btilde = b;

% creating initial condition sets -  ensure that 0 does not belong to set
% in any of the coordinates
X0_center = [ 10; 20; 30 ];
X0_gen = 5*eye( 3 );
Zono_X0 = zonotope( X0_center, X0_gen );

% tVec represents eigenvectors in transformed coordinate system- i.e.,
% standard basis
tVec = eye(3);

tX0_center = orthCol' * X0_center;
tX0_gen = orthCol' * X0_gen;
Zono_tX0 = zonotope( tX0_center, tX0_gen );

outerintervals_tX0 = multi_outer( mptPolytope( Zono_tX0 ), ...
    mptPolytope( Zono_tX0 ), epsilon );

%% Algorithm 4 of paper
for i = 1:length(time)
    
    for k = 1:size( outerintervals_tX0, 2 )
        
        zonoInt = zonotope( outerintervals_tX0( :, k ) );
        outerZono_tX0_gen = zonoInt.generators;
        outerZono_tX0_center = zonoInt.center;
        R{ i, k } = 0;
        cen = zeros(3,1);
        gen = [];
  
        for m = 1:size( outerZono_tX0_gen, 2 )
            upp=solveHypergeom(3,btilde(m),lambda(m),outerZono_tX0_center(m)+outerZono_tX0_gen(m, m),time(i)); 
            mustBeReal(upp); % if the problem terminates here, change the initial conditions: singularity present in RS
            low = solveHypergeom(3,btilde(m),lambda(m),outerZono_tX0_center(m)-outerZono_tX0_gen(m, m),time(i)); 
            mustBeReal(low);
            cen(m) = (upp+low)/2;
            gen = [gen,(upp-low)*tVec(:,m)/2];
        end
        R{ i, k } = zonotope( cen, gen );
    end
end


% flipping to different coordinate system
S = R;
for i = 1:size(R,1)
    for k = 1:size(R,2)
        c = R{i, k}.center;
        g = R{i, k}.generators;
        S{i,k} = zonotope( orthCol * c, orthCol * g );
    end
end

%% Setup for CORA
homoPoly =@(x,u) [...
    [1 0 0]*(reshape(reshape(odecoTensor, [9, 3])*[x(1); x(2); x(3)], 3, 3)*[x(1); x(2); x(3)]+b);
    [0 1 0]*(reshape(reshape(odecoTensor, [9, 3])*[x(1); x(2); x(3)], 3, 3)*[x(1); x(2); x(3)]+b);
    [0 0 1]*(reshape(reshape(odecoTensor, [9, 3])*[x(1); x(2); x(3)], 3, 3)*[x(1); x(2); x(3)]+b);
    ];

tank = nonlinearSys('tank',homoPoly);

% simulating examples - proxy for the ground truth
hparams.R0 = zonotope( [ X0_center,  X0_gen ] );
hparams.U = zonotope([0]);
hparams.tStart = 0;
hparams.tFinal = T-0.01;

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

simOpt.points = 75;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;
hsys = nonlinearSys(homoPoly);

% simulating the randomly initialized points
simRes = simulateRandom(hsys, hparams, simOpt);

% Generating reachable sets using CORA 
R_cora = reach(tank, hparams, options);

dims = {[1 2] [2 3] [3 1]};
dim_labels = {'x', 'y', 'z'}; % Labels for dimensions

%% Plotting the results - 3 dimensional systems; modify as needed
for k = 1:length(dims)
    projDim = dims{k};
    dim_label_x = dim_labels{projDim(1)};
    dim_label_y = dim_labels{projDim(2)};
    
    figure; hold on;

    % plot reachable sets
    for i = 1:length(time)-1
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
