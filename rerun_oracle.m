

ex = 7

parfor r = 1:num_runs
  % temporary vectors for storing acc and runtime
%  acc_		= acc(7,r);%, num_params);
  
%  for i = 1:num_params
	% make experiments reproducable
	rng(r);%*num_params+i);
	
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
	
	% run lsgm
	rng(r);
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'rand');
	[match clust_labels] = BigGMr( AA,BB,m, numdim, max_clust, @spectralEmbedElbow, @kmeansAlgr, gmAlg);
	acc(ex,r) = mean(clust_labels(nonseeds,1)==clust_labels(nonseeds,2));
	
%  acc(:,r) = acc_;
end

save(save_file_name, 'num_seeds', 'max_clust', 'acc', 'runtime');
