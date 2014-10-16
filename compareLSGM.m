function [ acc, runtime] = compareLSGM(num_runs, N, m, corrln, lam)
% Computes the accuracy and runtime of LSGM, SGM, U, RANK, QCV, rand, PATH
%Input: numsim: number of simulations
%   N: number of vertices
%   m: number of seeds
%   corrln: correlation
%   lam: block probability matrix
%Output: acc accuracy of 6 methods
%   runtime: computes the run time of 6 methods

%numdim = rank(lam);
numdim = 10;
max_clust_size = sum(N)/length(lam);
acc = zeros(6, num_runs);
runtime = zeros(6, num_runs);
for r = 1:num_runs
	% make experiments reproducable
	rng(r);
   
	% generate correlated graphs
	[A, B, shuffle] = sampleGraphs(m, N, corrln, lam);
	ex = 1;
	
	% lsgm
	start = tic;
	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, @seedgraphmatchell2);
	runtime(ex,r)  = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
	ex = ex+1;
	
%	start = tic;
%	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, @seedgraphmatchell2test);
%	runtime(ex,r)  = toc(start);
%	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
%	ex = ex+1;

	start = tic;
	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, @graphmatchell2);
	runtime(ex,r) = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
	ex = ex+1;

	% lsgmm
	% choices for method: I U RANK QCV rand PATH s
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'U');
	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, gmAlg);
	runtime(ex,r) = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
	ex = ex+1;
	
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'RANK');
	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, gmAlg);
	runtime(ex,r) = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
	ex = ex+1;
	
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'QCV');
	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, gmAlg);
	runtime(ex,r) = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
	ex = ex+1;
	
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'rand');
	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, gmAlg);
	runtime(ex,r) = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
	ex = ex+1;
	
	start = tic;
	gmAlg = @(A,B,s) graphmAlg(A,B,s,'PATH');
	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, gmAlg);
	runtime(ex,r) = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
	ex = ex+1;
end



end

