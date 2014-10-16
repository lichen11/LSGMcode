matlabpool(12);

% add path of other files
addpath './algorithms/'
addpath './preprocessing/'

%% compile methods (eg. PATH)
%cd './algorithms/graphm-0.52'
%system('make clean');
%system('./graphm_install');
%cd '../../'
%
%% compile GLAG
%cd './algorithms/GLAG'
%mex vector_th_alpha_beta.cpp;
%cd '../../'

%% preprocess brain graph
%start = tic;
%down = 4;
%directed = false;
%thresh = 0;
%
%load('data/KKI-21_KKI2009-08_big_graph.mat');
%A = downsample(fibergraph,down,directed);
%clear fibergraph;
%load('data/KKI-21_KKI2009-26_big_graph.mat');
%B = downsample(fibergraph,down,directed);
%clear fibergraph;
%
%v1 = findLCC(A,directed);
%v2 = findLCC(B,directed);
%v = intersect(v1,v2);
%A = double(A(v,v)>thresh);%A(v,v);%
%B = double(B(v,v)>thresh);%B(v,v);%
%
%disp( ['done downsampling: ', num2str(toc(start))] );

% load preprocessed brain graph
load('KKI-08-29.mat');


% name of save file
save_file_name = 'lsgm_data_08-29_2500.mat';


% simulation parameters
num_runs = 20;
num_exp = 7;%

% max cluster size
max_clust = 800;%150;%500%
%% number of seeds
num_seeds_vec = 2500;%[200, 1000, 2000];
num_params = length(num_seeds_vec);
%m = num_seeds;
%% index vector for nonseed vertices
%nonseeds = m+1:N;

% number of vertices
N = size(A,1);

% maximum dimension considered for the embedding
numdim = 100;

acc     = zeros(num_exp, num_params, num_runs);
runtime = zeros(num_exp, num_params, num_runs);
parfor r = 1:num_runs
  % temporary vectors for storing acc and runtime
  acc_		= zeros(num_exp, num_params);
  runtime_	= zeros(num_exp, num_params);
  
  for i = 1:num_params
	% make experiments reproducable
	rng(r*num_params+i);
	
	%load parameters
	% index vector for nonseed vertices
	num_seeds = num_seeds_vec(i);
	m = num_seeds;
	nonseeds = m+1:N;
	
	[r num_seeds]
	
	% random seeding
	Aperm = randperm(N);
	seedinds = Aperm(1:num_seeds);
	Bperm = randperm(N);
	% move vertices in B to match seeds in A
	for seed_ind = 1:num_seeds
		seed = seedinds(seed_ind);
		Bperm(Bperm==seed) = Bperm(seed_ind);
		Bperm(seed_ind) = seed;
	end
	Bend = Bperm(num_seeds+1:end);
	% create shuffled matrices
	AA = A(Aperm,Aperm);
	BB = B(Bperm,Bperm);
	
	% run experiments
	ex = 1;
	
	% lsgm
	start = tic;
	rng(r*num_params+i);
	[match oracle_acc] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, @seedgraphmatchell2);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
	ex = ex+1;

	start = tic;
	rng(r*num_params+i);
	[match oracle_acc] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, @graphmatchell2);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
	ex = ex+1;
	
	% lsgmm
	% choices for method: I U RANK QCV rand PATH s
	start = tic;
	rng(r*num_params+i);
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'rand');
	[match oracle_acc] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
	ex = ex+1;
	
	start = tic;
	rng(r*num_params+i);
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'U');
	[match oracle_acc] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
	ex = ex+1;
	
	start = tic;
	rng(r*num_params+i);
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'RANK');
	[match oracle_acc] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
	ex = ex+1;
	
	start = tic;
	rng(r*num_params+i);
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'QCV');
	[match oracle_acc] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
	ex = ex+1;
	
%	start = tic;
%	rng(r*num_params+i);
%	gmAlg = @(A,B,s) graphmAlg(A,B,s,'PATH');
%	[match oracle_acc] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, gmAlg);
%	runtime_(ex,i) = toc(start);
%	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
%	ex = ex+1;
%	
%	% GLAG algorithm
%	start = tic;
%	rng(r*num_params+i);
%	[match oracle_acc] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbed, @kmeansAlgr, @graphmGLAG);
%	runtime_(ex,i) = toc(start);
%	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
%	ex = ex+1;	
%	% remeber to change num_exp when adding an experiment
%	
	% oracle accuracy
	acc_(ex,i) = oracle_acc;
	ex = ex+1;
	
  end
  runtime(:,:,r) = runtime_;
  acc(:,:,r) = acc_;
end

save(save_file_name, 'num_seeds_vec', 'max_clust', 'acc', 'runtime');
exit
