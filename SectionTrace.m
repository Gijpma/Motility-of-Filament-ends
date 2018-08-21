function [points] = SectionTrace(X,Y, param)
%% sections the traces into stopped and moving sections

  if length(X)>1
            error=10;
            points(1,:)=[X(1) Y(1) 1];
            points(2,:)=[X(end) Y(end) length(X)];
            % this loop gradually breaks up the line connecting the start and
            % end point of the position trace by searching for the points on
            % the graph that are furthest away from this line, and subsequently
            % drawing lines through these far away points.
            while error>param.errormax
                clear temppoint
                temppoint(1,:)=[X(1) Y(1) 1];  % always start at start point
                error=0;
                for i=1:length(points(:,1))-1
                    % first find the slope and intercept of the existing line
                    a=(points(i,2)-points(i+1,2))/(points(i,1)-points(i+1,1));
                    b=points(i,2)-a*points(i,1);
                    if isnan(a)      % Isnan occurs when line is vertical. at
                        % this point error defaults to 0
                        error=0;
                    else
                        %finds the maximum distances to the line given by ax+b
                        [maxout,maxspot]=max(Y(points(i,3):points(i+1,3))-a*X(points(i,3):points(i+1,3))-b);
                        [minout,minspot]=min(Y(points(i,3):points(i+1,3))-a*X(points(i,3):points(i+1,3))-b);
                        if max([abs(minout) abs(maxout)])>param.errormax
                            error=max([abs(minout) abs(maxout)]);   %error defined as the max point distance
                        end
                        % this handles the three cases of the maximum and
                        % minimum relative locations and their possible
                        % existence/nonexistence (error below treshold)
                        if abs(maxout)>=param.errormax & abs(minout)>=param.errormax
                            k=points(i,3)+min([maxspot minspot])-1;
                            temppoint(end+1,:)=[X(k) Y(k) k];
                            k=points(i,3)+max([maxspot minspot])-1;
                            temppoint(end+1,:)=[X(k) Y(k) k];
                        elseif abs(maxout)<=param.errormax & abs(minout)>=param.errormax
                            k=points(i,3)+minspot-1;
                            temppoint(end+1,:)=[X(k) Y(k) k];
                        elseif abs(maxout)>=param.errormax & abs(minout)<=param.errormax
                            k=points(i,3)+maxspot-1;
                            temppoint(end+1,:)=[X(k) Y(k) k];
                        end
                    end
                    % finish with the end point
                    temppoint(end+1,:)=[X(points(i+1,3)) Y(points(i+1,3)) points(i+1,3)];
                    
                end
                points=temppoint;
                clear temppoint
                if length(points)>150 % if this happens something probably went wrong
                    pause
                end
                
            end
            
        else
            points=[0 0];
            
        end
        