function H = Q_entropy_2D(M)
    [v_B, val] = eig(M);
    p_B = [val(1,1) val(2,2)];
    A = zeros(2);
    for i = [1,2]
        if p_B(i) ~= 0
            A = A + p_B(i)*log2(p_B(i))*v_B(:,i)*v_B(:,i)';
        end
    end
    H = -trace(A);
end