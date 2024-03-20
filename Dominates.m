function b=Dominates(x,y)

    if isstruct(x)
        x=x.Cost;
    end

    if isstruct(y)
        y=y.Cost;
    end

    % size(x)
    % size(y)
    b=all(x<=y) && any(x<y);


end