function r = nmf2dCluster(x, fs, method, maxiter,d,Tau,Phi)




if method=='SAGE'
    [W, H, ~] = is_nmf2D_mu(cochleagram(gammatoneFast(x,128,[50 fs/2])),maxiter,d,Tau,Phi);% the MU algorithm
elseif method=='MU'
    [W, H, ~] = is_nmf2D_em(cochleagram(gammatoneFast(x,128,[50 fs/2])),maxiter,d,Tau,Phi);% the MU algorithm
end

for idx = 1:d
    Rec(:,:,idx) = isp_nmf2d_rec(W, H, Tau,Phi, idx);
end
for k = 1:size(Rec,3)
    mask = logical(Rec(:,:,k)==max(Rec,[],3));
    r(k,:) = synthesisFast(x,mask,[50 fs/2]);
end
