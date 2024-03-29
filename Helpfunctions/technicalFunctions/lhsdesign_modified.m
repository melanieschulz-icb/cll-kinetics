function [X_scaled, X_normalized] = lhsdesign_modified(n, min_ranges_p, max_ranges_p)
    p=length(min_ranges_p);
    [M,N]=size(min_ranges_p);
    if M<N
        min_ranges_p=min_ranges_p';
    end
    
    [M,N]=size(max_ranges_p);
    if M<N
        max_ranges_p=max_ranges_p';
    end
    
    slope=max_ranges_p-min_ranges_p;
    offset=min_ranges_p;
    
    SLOPE=ones(n,p);
    OFFSET=ones(n,p);
    
    for i=1:p
        SLOPE(:,i)=ones(n,1).*slope(i);
        OFFSET(:,i)=ones(n,1).*offset(i);
    end
    X_normalized=lhsdesign(n,p);
    X_scaled=SLOPE.*X_normalized+OFFSET;
    %X_scaled=X_scaled+OFFSET;
end