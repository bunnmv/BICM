
function u_hat = block_dec_sdd(code, r)
    
    corr = -r * code.vocabulario';
    [~, indx] = max(corr,[],2);
    c_hat = code.vocabulario(indx,:);
    u_hat = c_hat(:,1:code.k)>0;

end
