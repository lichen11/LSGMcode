function [ match, clust_labels ] = BigGM( A, B, s, numdim, numclust, embedAlg, clustAlg, graphMatchAlg)


start = tic;
sumn = length(A)-s;

% perform embedding
startt = tic;
[XA XB] = embedAlg(A, B, numdim);
fprintf( 'done projection: %f\n', toc(startt) );

% compute procrusties othogonal projection (on the seed vertices)
startt = tic;
[~,~,TRANSFORM]=procrustes(XA(1:s,:),XB(1:s,:));
TRANSFORM.c=ones(sumn+s,1)*TRANSFORM.c(1,:);
XB = TRANSFORM.b * XB * TRANSFORM.T + TRANSFORM.c;
fprintf( 'done procrusties: %f\n', toc(startt) );


% cluster using the embedding
startt = tic;
XAXB=[XA;XB];
nonseedsA = s+1:s+sumn;
nonseedsB = s+sumn+ s+1:2*(s+sumn);
[IDX, centroid, Dis] = clustAlg(XAXB, numclust);
fprintf( 'done clustering: %f\n', toc(startt) );

% fixing cluster sizes to be equal in both graphs
startt = tic;
IDXA=IDX(nonseedsA,:);
DisA=Dis(nonseedsA,:);
IDXB=IDX(nonseedsB,:);
DisB=Dis(nonseedsB,:);

clear IDX Dis

clustsizesA=zeros(numclust,1);
clustsizesB=zeros(numclust,1);
for i=1:numclust
   iiAi=find(IDXA==i);
   iiBi=find(IDXB==i);
   clustsizesA(i)=length(iiAi);
   clustsizesB(i)=length(iiBi);
end

clustsizes=round( (clustsizesA+clustsizesB)/2  );

% add or subtract 1 from largest cluster sizes
[val,ind] = sort(clustsizes, 1, 'descend');
temp = sumn -sum(clustsizes);
mask = ind(1:abs(temp));
clustsizes(mask) = clustsizes(mask) +sign(temp);

% display some cluster sizes
%sizes_of_some_clusters = clustsizes(1:min(end,20))'

clust_labels = zeros(s+sumn,2,numclust);

% save clusters(subgraphs) to match (so it can be parallelized)
pieceA_ = cell(numclust,1);
pieceB_ = cell(numclust,1);
gA_ = cell(numclust,1);
gB_ = cell(numclust,1);

for i=1:numclust
	if clustsizes(i) == 0
		continue
	end
	% sort by distance
	[~, TA] = sort(DisA(:,i));
	[~, TB] = sort(DisB(:,i));
	% use closest vertices
    gA=TA(1:clustsizes(i));
    gB=TB(1:clustsizes(i));
    % make these vertices unselectable for next time
	DisA(gA,:)=inf;
	DisB(gB,:)=inf;
    
    gA = gA+s;
    gB = gB+s;
    
    % sav cluster labels
    clust_labels(gA,1) = numclust;
    clust_labels(gB,2) = numclust;
    
    gAaug=[1:s gA'];
    gBaug=[1:s gB'];
    
    pieceA=A(gAaug,gAaug);
    pieceB=B(gBaug,gBaug);
    
    pieceA_{i} = pieceA;
    pieceB_{i} = pieceB;
    
    gA_{i} = gA;
    gB_{i} = gB;
end

clear A B DisA DisB
% perform graph matching in parallel
match = zeros(s+sumn,numclust);
parfor i = 1:numclust
	% load subgraph adjacency matrix
	gA = gA_{i};
	gB = gB_{i};
	pieceA = pieceA_{i};
	pieceB = pieceB_{i};
	
	% if graph is empty
	if length(gA)==0
		continue;
	end
	
    % perform graph match
	ord = graphMatchAlg(pieceA, pieceB, s);
	ord = ord(s+1:end)-s
	
	% save results
	temp = zeros(s+sumn,1);
	temp(gA) = gB(ord);
    match(:,i) = temp;
end
% combine results
match(1:s,1) = 1:s;
match = sum(match,2)';

fprintf( 'done matching: %f\n', toc(startt) );
fprintf( 'total time: %f\n', toc(start) );

end
