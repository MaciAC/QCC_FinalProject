function H = Q_entropy_2D(M)
    dim = length(M);
    [v_B, val] = eig(M);
    A = zeros(dim);
    for i = [1:dim]
        p_B = val(i,i);
        if p_B ~= 0
            A = A + p_B*log2(p_B)*v_B(:,i)*v_B(:,i)';
        end
    end
    H = -trace(A);
end