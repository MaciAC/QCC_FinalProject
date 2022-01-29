function H = Q_entropy_2D(M)
    [v_B, val] = eig(M);
    p_B = [val(1,1) val(2,2)];
    H = -trace( p_B(1)*log2(p_B(1))*v_B(:,1)*v_B(:,1)' + p_B(2)*log2(p_B(2))*v_B(:,2)*v_B(:,2)');
end