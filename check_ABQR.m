%% function for checking properties

function check_ABQR(A,B,Q,R)

if min(eig(R)) > 0
    display('R is positive definite')
elseif min(eig(R)) == 0
    display('R is positive semidefinite -- DLQR.M will fail!')
else
    display('R is not positive semidefinite  -- DLQR.M will fail!')
end

if min(eig(Q)) > 0
    display('Q is positive definite')
elseif min(eig(Q)) == 0
    display('Q is positive semidefinite')
    % check detectability
    [ABAR,BBAR,CBAR,T,K] = obsvf(A,B,sqrtm(Q));
    nobs = size(A,2) - sum(K);
    if max(abs(eig(ABAR(1:nobs,1:nobs)))) == 0
        display('(Q^0.5,A) is constructible')
    elseif max(abs(eig(ABAR(~K,~K)))) < 1
        display('(Q^0.5,A) is detectable')
    else
        display('(Q^0.5,A) is not detectable -- DLQR.M will fail!')
    end
    
    
else
    display('Q is not positive semidefinite -- DQLR.M will fail!')
end

r = rank(ctrb(A,B));
if r == 2
    display('(A,B) is reachable')
else
    display('(A,B) is not reachable')
    % check stabilizability
    [ABAR,BBAR,CBAR,T,K] = ctrbf(A,B,sqrtm(Q));
    nrea = size(A,2) - sum(K);
    if max(abs(eig(ABAR(1:nrea,1:nrea)))) == 0
        display('(A,B) is controllable')
    elseif max(abs(eig(ABAR(~K,~K)))) < 1
        display('(A,B) is stabilizable')
    else
        display('(A,B) is not stabilizable -- DLQR.M will fail!')
    end
end

end