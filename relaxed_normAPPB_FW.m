function [P,Pp]=relaxed_normAPPB_FW(A,B,varargin)

verbose = 0;

AtA = A'*A;
BBt = B*B';
C=AtA;
D=BBt;

p=size(A,1);

f1 = @(P) norm(A*P-P*B,'fro')^2;

tol=5e-2;
tol2=1e-6;

P=ones(p)/p;

 if ~isempty(varargin)
     if (size(varargin{1},1) > 1)
        P=varargin{1}; 
     else
        tol2=varargin{1};    
     end
     
 end
 

f=f1(P);
var=1;

while (f>tol) && (var > tol2)
    fold=f;

    grad = AtA*P -A'*P*B - A*P*B' + P*BBt;
    
    corr=lapjv(grad,0.01);
    Ps=perm2mat(corr);

    %Ps=hungarian_mex(grad*1e3);
    %Ps =Ps';
    
    C = A*(P-Ps) + (Ps-P)*B; 
    D = A*Ps-Ps*B;
    aopt = -0.5*trace(C*D'+D*C')/trace(C*C');

     
    Ps4 = aopt*P + (1-aopt)*Ps;
    
    f=f1(Ps4);
    P=Ps4;
    
    if (f>fold)
        fprintf('it shouldnt');
        f=fold;
        f2 = f1(Ps);
        
        if (f2 < fold)
            f=f2;
            P=Ps;
        end
        
    end

    var=abs(f-fold);
    if (verbose) fprintf('f: %1.5f  var %.15f\n',f,var); end

end

    corr=lapjv(-P,0.01);
    Pp=perm2mat(corr);


end
