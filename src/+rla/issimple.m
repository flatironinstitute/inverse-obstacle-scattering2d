function out = issimple(xs,ys,nover)
    if(nargin == 2) 
       noveruse = 5;
    elseif (nargin == 3)
       noveruse = ceil(nover);
    end

    nn = numel(xs);
    x_ft = fft(xs);
    y_ft = fft(ys);
    x_big_ft = zeros([noveruse*nn,1]);
    y_big_ft = zeros([noveruse*nn,1]);
    x_big_ft(1:(nn/2)) = x_ft(1:(nn/2));
    y_big_ft(1:(nn/2)) = y_ft(1:(nn/2));
    x_big_ft(noveruse*nn -nn/2 +1:end) = x_ft((nn/2+1):end);
    y_big_ft(noveruse*nn -nn/2 +1:end) = y_ft((nn/2+1):end);
    xbig = real(ifft(x_big_ft))*(noveruse);
    ybig = real(ifft(y_big_ft))*(noveruse);
    %plot(xs,ys,'b.',xbig,ybig,'r.')
    
    P = polyshape(xbig,ybig,'KeepCollinearPoints',true);
%     fprintf('NumRegions=%d\n',P.NumRegions)
    if P.NumRegions>1
        out = false;
    else
        out = true;
    end
    
    return
    
