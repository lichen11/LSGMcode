% add path of other files
addpath './algorithms/'


% simulation parameters
num_runs = 200;
num_exp = 2;

%K = 4:8;
num_blocks = 8;
num_cores = 1:4;
max_clust = 200;
num_params = length(num_cores);%max_clust_vec);%size(Ns,2);
N = 200*ones(num_blocks,1);
% number of seed vertices
%ms = 2*ones(num_blocks,1);
m = 20;
%nonseeds = m+1:m+sum(N);
% correlation between graphs
corrln = .9;
% make lambda matrix
%lam=  .3*eye(num_blocks)+.3*ones(num_blocks);

save_file_name = './scalability/desktop_simulate_corr_9_cores_1_4_paired_200_sim_200.mat';%'lsgm-sim-corr_6-cores_4_8_200_paired.mat';

%matlabpool close;
matlabpool(2);
acc = zeros(num_exp, num_params, num_runs);
runtime = zeros(num_exp, num_params, num_runs);
for r = 1:num_runs
  % initalize
  	acc_ = zeros(num_exp, num_params);
  	runtime_ = zeros(num_exp, num_params);
	nonseeds = m+1:m+sum(N);
	
	lam=  .3*eye(num_blocks)+.3*ones(num_blocks);
	numdim = rank(lam);
%	max_clust = max_clust_vec(i);
	
	% generate correlated graphs
	rng(r);	
	[A, B, shuffle] = sampleGraphs(m, N, corrln, lam);
  for i = 1:num_params
	
	num_core = num_cores(i);
	
	% set up parfor
	matlabpool close
	matlabpool(num_core);
	fprintf('\nnumber of cores = %d\n', num_core);

	ex = 1;
	
	% lsgm
	start = tic;
    [match clust_labels] = BigGMforScala( A,B,m, numdim, max_clust);%, @spectralEmbed, @kmeans0, @seedgraphmatchell2); %@kmeansAlg, @seedgraphmatchell2);
	%[match clust_labels] = BigGM( A,B,m, numdim, max_clust, @spectralEmbed, @kmeans0, @seedgraphmatchell2); %@kmeansAlg, @seedgraphmatchell2);
	runtime_(ex,i)  = toc(start);
	acc_(ex,i) = mean(shuffle(nonseeds)==match(nonseeds));
	ex = ex+1;

	
	% oracle accuracy
	acc_(ex,i) = mean(clust_labels(nonseeds,1)==clust_labels(shuffle(nonseeds),2));
	ex = ex+1;

  end
  runtime(:,:,r) = runtime_;
  acc(:,:,r) = acc_;
end

save(save_file_name, 'num_cores', 'm', 'N', 'max_clust', 'corrln', 'lam', 'acc', 'runtime');



%exit
