function [angleout,Centroid]=orientation(x,y)
            
            % Adapted from RegionProps function in Matlab
            
            Centroid=[mean(x) mean(y)];
            x = x - mean(x);
            y = -(y - mean(y)); % This is negative for the
            % orientation calculation (measured in the
            % counter-clockwise direction).

            N = length(x);

            % Calculate normalized second central moments for the region. 1/12 is
            % the normalized second central moment of a pixel with unit length.
            uxx = sum(x.^2)/N + 1/12;
            uyy = sum(y.^2)/N + 1/12;
            uxy = sum(x.*y)/N;

            % Calculate orientation.
            if (uyy > uxx)
                num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
                den = 2*uxy;
            else
                num = 2*uxy;
                den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
            end
            if (num == 0) && (den == 0)
                angleout = 0;
            else
                angleout = (180/pi) * atan(num/den);
            end