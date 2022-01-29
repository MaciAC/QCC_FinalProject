function m = PartialTrace(M,option)

    if option == 1
        m(1,1) = M(1,1) + M(1,3);
        m(1,2) = M(1,2) + M(3,4);
        m(2,1) = M(2,1) + M(4,3);
        m(2,2) = M(2,2) + M(4,4);
    elseif option == 2
        m(1,1) = M(1,1) + M(1,3);
        m(1,2) = M(1,2) + M(3,4);
        m(2,1) = M(2,1) + M(4,3);
        m(2,2) = M(2,2) + M(4,4);
    end
end
