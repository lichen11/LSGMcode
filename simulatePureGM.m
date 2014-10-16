% add path of other files
addpath './algorithms/'

%%%% compile methods
%%%% run the compile method once.
%cd './algorithms/graphm-0.52'
%system('make clean')
%system('./graphm_install')
%cd '../../'
%%% compile for glag: first cd into the glag directory
%mex vector_th_alpha_beta.cpp


% simulation parameters
num_runs = 20;
num_exp = 1;

% matrix for vector of N for each experiment
num_blocks = 2;
% each column of Ns is one experiment ie. N = Ns(:,i)
Ns = repmat([50, 100, 150, 200], [num_blocks,1]);
%Ns = repmat([250, 300, 400, 500], [2,1]);
num_params = size(Ns,2);
% number of seed vertices
%ms = 10*ones(num_params,1);
ms = zeros(num_params,1);
% correlation between graphs
corrln = .9;
% make lambda matrix
lam=  .3*eye(num_blocks)+.3*ones(num_blocks);

save_file_name = 'gm_small_2blocks_corrln_9_glag.mat';

numdim = rank(lam);
%numdim = 10;
num_clust = size(Ns,1);
acc = zeros(num_exp, num_params, num_runs);
runtime = zeros(num_exp, num_params, num_runs);

parfor r = 1:num_runs
    %r  
    % initalize
    acc_ = zeros(num_exp, num_params);
    runtime_ = zeros(num_exp, num_params);
    for i = 1:num_params
	    % make experiments reproducable
	    rng(r*num_params+i);

	    % extract parameters
	    N = Ns(:,i);
	    m = ms(i);
	    nonseeds = m+1:m+sum(N);
	
	    % generate correlated graphs
	    [A, B, shuffle] = sampleGraphs(m, N, corrln, lam);
	    ex = 1;
	
%	    % Frank-wolfe
%	    start = tic;
%        match = graphmatchell2(A, B, m);
%	    runtime_(ex, i) = toc(start);
%	    acc_(ex, i) = mean(shuffle(nonseeds)==match(nonseeds));
%	    ex = ex+1;
%	    
%	    % random guesses
%	    start = tic;
%	    match =  graphmAlg(A, B, m, 'rand');
%	    runtime_(ex,i) = toc(start);
%	    acc_(ex, i) = mean(shuffle(nonseeds)==match(nonseeds));
%	    ex = ex+1;
%	
%	    % Umeyaha
%	    start = tic;
%	    match =  graphmAlg(A, B, m, 'U');
%	    runtime_(ex, i) = toc(start);
%	    acc_(ex, i) = mean(shuffle(nonseeds)==match(nonseeds));
%	    ex = ex+1;
%	
%	    % RANK
%	    start = tic;
%	    match =  graphmAlg(A, B, m, 'RANK');
%	    runtime_(ex, i) = toc(start);
%	    acc_(ex, i) = mean(shuffle(nonseeds)==match(nonseeds));
%	    ex = ex+1;
%	    
%	    % QCV
%	    start = tic;
%	    match =  graphmAlg(A, B, m, 'QCV');
%	    runtime_(ex, i) = toc(start);
%	    acc_(ex, i) = mean(shuffle(nonseeds)==match(nonseeds));
%	    ex = ex+1;
%	
%	    % PATH
%	    start = tic;
%	    match =  graphmAlg(A, B, m, 'PATH');
%	    runtime_(ex,i) = toc(start);
%	    acc_(ex, i) = mean(shuffle(nonseeds)==match(nonseeds));
%	    ex = ex+1;
	
	    % GLAG 
	    start = tic;
	    match =  graphmGLAG(A, B, m);
	    runtime_(ex, i) = toc(start);
	    acc_(ex, i) = mean(shuffle(nonseeds)==match(nonseeds));
	    ex = ex+1;
    end 
    runtime(:,:,r) = runtime_;
    acc(:,:,r) = acc_;
end

save(save_file_name, 'ms', 'Ns', 'corrln', 'lam', 'acc', 'runtime');


