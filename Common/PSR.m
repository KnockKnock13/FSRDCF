function [ psr ] = PSR( map )
%caculating the psr of the correlation output
%   psr = (max(map)-mean(side lobe))/sigma(side lobe)
%   the side lobe of the correlation output is the pixels exclude the 11x11
%   window around the peak
    [m,n]=size(map);
    g = max(map(:));
    [x,y] = find(map==g);
    count=0;
    for i=-5:5
        for j=-5:5
            k = x+i;
            p = y+j;
            if(k>0 && p>0 && k<=m && p<=n)
                map(k,p)=0;
                count =count +1;
            end
        end
    end
    mu = sum(map(:))/(numel(map)-count);
    map =power(map-mu,2);
    sigma = (sum(map(:))-count*mu^2)/(numel(map)-count);
    psr = (g-mu)/sqrt(sigma);
end

