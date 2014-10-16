% add path of other files
addpath './algorithms/'

%% compile methods (eg. PATH)
%cd './algorithms/graphm-0.52'
%system('make clean');
%system('./graphm_install');
%cd '../../'

%% compile GLAG
%cd './algorithms/GLAG'
%mex vector_th_alpha_beta.cpp;
%cd '../../'


% simulation parameters
num_runs = 100;
num_exp = 8;

% matrix for vector of N for each experiment ie. N = Ns(:,i)
num_blocks = 5;
Ns = repmat([50, 100, 150, 200], [num_blocks,1]);
%Ns = repmat([250, 300, 400, 500], [2,1]);
num_params = size(Ns,2);
% number of seed vertices
ms = 10*ones(num_params,1);
max_clust_sizes = [50, 100, 150, 200]*1.1;

% correlation between graphs
corrln = .6;
% make lambda matrix
lam=  .3*eye(num_blocks)+.3*ones(num_blocks);

save_file_name = 'lsgmr_small_5blocks.mat';

numdim = rank(lam);
%numdim = 10;
acc = zeros(num_exp, num_params, num_runs);
runtime = zeros(num_exp, num_params, num_runs);
parfor r = 1:num_runs
  r
  
  % initalize
  acc_ = zeros(num_exp, num_params);
  runtime_ = zeros(num_exp, num_params);
  for i = 1:num_params
	% make experiments reproducable
	rng(r*num_params+i);
	
	% extract
	N = Ns(:,i);
	m = ms(i);
	max_clust = max_clust_sizes(i);
	nonseeds = m+1:m+sum(N);
	
	% generate correlated graphs
	[A, B, shuffle] = sampleGraphs(m, N, corrln, lam);
	ex = 1;
	
	% lsgm
	start = tic;
	match = BigGMr( A,B,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, @seedgraphmatchell2);
	runtime_(ex,i)  = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;
	
	start = tic;
	match = BigGMr( A,B,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, @graphmatchell2);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;

	% lsgmm
	% choices for method: I U RANK QCV rand PATH s
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'rand');
	match = BigGMr( A,B,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;
	
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'U');
	match = BigGMr( A,B,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;
	
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'RANK');
	match = BigGMr( A,B,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;
	
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'QCV');
	match = BigGMr( A,B,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;
	
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'PATH');
	match = BigGMr( A,B,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;
	
	% GLAG algorithm
	start = tic;
	match = BigGMr( A,B,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, @graphmGLAG);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;
	
	% remeber to change num_exp when adding an experiment
  end
  runtime(:,:,r) = runtime_;
  acc(:,:,r) = acc_;
end

save(save_file_name, 'ms', 'Ns', 'corrln', 'lam', 'acc', 'runtime');
