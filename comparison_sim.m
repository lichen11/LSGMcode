num_runs = 25;
% we use the same parameters in Figure 2 in the LSGM paper
% vector of vertices in each block
K = 900;
N = 30 * ones(K,1);
% number of seed vertices
m_vec = 10:10:50;
% correlation between graphs
corrln = .9;
%alpha=.5;
%lam=alpha*[.5,.3,.4;.3,.8,.6;.4,.6,.3]+.5*(1-alpha)*ones(3);
%lam = [0.6 0.3 0.2; 0.3 0.7 0.3; 0.2 0.3 0.7];
lam = rand(K);
lam=  triu(lam) + triu(lam,1)';
acc_tensor = zeros(7, num_runs, length(m_vec));
runtime_tensor = zeros(7, num_runs, length(m_vec));
for i = 1:length(m_vec)
    m = m_vec(i);
    [acc, runtime] = compareLSGM(num_runs, N, m, corrln, lam);
    acc_tensor(:,:,i) = acc;
    runtime_tensor(:,:,i) = runtime;
end

save('result.mat', 'acc_tensor', 'runtime_tensor', 'lam', 'K', 'N', 'm_vec', 'corrln', 'num_runs');