function [A, B, findme] = sampleGraphs(s, n, crln, lam)
% creates 2 correlated k block model graphs

numblocks=length(n);
sumn=sum(n);

% create block membership vector
indicate=[];
for i=1:numblocks
    indicate=[indicate; i*ones(n(i),1)];
end
% randomly generate seed memberships
for i=1:s
    temp=0;
    a=rand;
    for i=1:numblocks
        if a>temp
            j=i;
        end
        temp=temp+n(i)/sumn;
    end
    indicate=[j; indicate];
end
% shuffle the nonseeds
mix=[1:s randperm(sumn) ]';
indicate=indicate(mix);

nsize = zeros(numblocks,1);
for i=1:numblocks
	nsize(i) = sum(indicate==i);
end

% generate A and B vectors
A=zeros(sumn+s,sumn+s);
B=zeros(sumn+s,sumn+s);

for ii=1:numblocks
	maski = indicate==ii;
	for jj = 1:numblocks
		maskj = indicate==jj;
		A(maski,maskj) = rand(nsize(ii),nsize(jj)) < lam(ii,jj);
		B(maski,maskj) = rand(nsize(ii),nsize(jj)) < (  (1-crln)*lam(ii,jj)+A(maski,maskj)*crln  );
	end
end
% symmetrize A and B
A = triu(A,1);
B = triu(B,1);
for i=1:sumn+s
	A(i+find(A(i,i+1:end)),i) = 1;
	B(i+find(B(i,i+1:end)),i) = 1;
end
% shuffle A
findme=[ 1:s randperm(sumn)+s ];
A=A(findme,findme);
