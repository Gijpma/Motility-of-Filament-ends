function [hist,TotalDist,Trace,Kill]=TraceCreate(ImageUse,MaxIter,firstframe,lastframe,param)
%% Generates initial trace based on trace image. minimal processing

    if sum(size(ImageUse))<20   %indicates stopped filament end
        hist(1)=firstframe;
        hist(2)=lastframe;
        TotalDist(1)=0;
        TotalDist(2)=0;
        Trace(1:2,1)=round(param.Square/2);
        Trace(1:2,2)=round(param.Square/2);
        Kill=0;
    else
    % generate skeleton image
    Skeleton=(bwmorph(ImageUse,'thin',Inf));
    % now prune all the braches. Skeleton is gradually decreased in size by 
    % removing endpoints untill only 2 endpoints remain.
    % subsequently the longest two branches are preserved as
    % they are presumed to be part of the main trace
    BranchImage=Skeleton;           
    BranchImage=padarray(BranchImage,[1 1],0,'both');  %pad to allow for edge detected endpoints
    [x,y]=find(bwmorph(BranchImage,'endpoints'));
    NumEndpoint=sum(x);
    iterations=1;
    if NumEndpoint>2
        for i=1:length(x)
            BranchImage=Skeleton;
            BranchImage=padarray(BranchImage,[1 1],0,'both');   %pad to allow for edge detected endpoints
            endx(i,1)=x(i);             % store current end points
            endy(i,1)=y(i);
            BranchImage(x(i),y(i))=0;   % Remove current end points
            pixelfound=true;            
            j=1;
            %% NOT CLEAR WHY I DID THIS
            while pixelfound
                [xout,yout]=find(BranchImage(endx(i,j)-1:endx(i,j)+1,endy(i,j)-1:endy(i,j)+1));    
                if length(xout)==1
                    j=j+1;
                    endx(i,j)=endx(i,j-1)-2+xout;
                    endy(i,j)=endy(i,j-1)-2+yout;
                    BranchImage(endx(i,j),endy(i,j))=0;
                else
                    pixelfound=false;
                end
            end
            branchlength(i)=j;
        end    
   
        % removes the shortest branches
        [~,b]=sort(branchlength);
        for i=1:length(b)-2
            for j=1:length(endx(b(i),:))
                if endx(b(i),j)>1 & endy(b(i),j)>1
                Skeleton(endx(b(i),j)-1,endy(b(i),j)-1)=0;
                end
            end
        end
    end
    % finds the route of the filament over the skeleton, starting at one
    % endpoint and moving to the other
    Trace=find(Skeleton);
    
    [EndX,EndY] = find(bwmorph(Skeleton,'endpoints'));
    %if there are endpoints, and the longest branches weren't too long, and 
    % the trace is at least three long, start at lowest frame value
    if ~isempty(EndX) & length(Trace)>3 & iterations<max([10 round(MaxIter*length(Trace))]) 
        [~,Startpoint]=min(ImageUse(sub2ind(size(Skeleton),EndX,EndY)));
        [a,b]=ind2sub(size(Skeleton),Trace);
        coordinates=[a b];      % coordinates of the skeleton trace points
        % removes the startpoint
        coordinates(find(coordinates(:,1)==EndX(Startpoint) & coordinates(:,2)==EndY(Startpoint)),:)=[];
        clear coordinatesout
        coordinatesout(1,:)=[EndX(Startpoint) EndY(Startpoint)]'; %output coordinates
        % walks over trace, by finding the next closest point
        for walk=2:length(Trace)
            [~,closest]=min((coordinates(:,1)-coordinatesout(walk-1,1)).^2 + (coordinates(:,2)-coordinatesout(walk-1,2)).^2);
            coordinatesout(walk,:)=coordinates(closest,:);
            coordinates(closest,:)=[];
        end
        %  smooth the trace
        Trace=[smooth(coordinatesout(:,1),4) smooth(coordinatesout(:,2),4)];
        clear intensity
        %using the smoothed trace, calculate the historyimage intensity
        %values of the fractional pixels
        for i=1:length(Trace)
            T=double(ImageUse);
            % make sure no points of zero intensity are used in averaging
            [a, b]=find((T(floor(Trace(i,1)):ceil(Trace(i,1)),floor(Trace(i,2)):ceil(Trace(i,2))))==0);
            if ~isempty(a)
                for j=1:length(a)
                T(floor(Trace(i,1))-1+a(j),floor(Trace(i,2))-1+b(j))=max(max((T(floor(Trace(i,1)):ceil(Trace(i,1)),floor(Trace(i,2)):ceil(Trace(i,2))))));
                end
            end
            intensity(i)=T(floor(Trace(i,1)),floor(Trace(i,2)))*(1-mod(Trace(i,1),1))*(1-mod(Trace(i,2),1))+...
                T(ceil(Trace(i,1)),floor(Trace(i,2)))*mod(Trace(i,1),1)*(1-mod(Trace(i,2),1))+...
                T(ceil(Trace(i,1)),ceil(Trace(i,2)))*mod(Trace(i,1),1)*mod(Trace(i,2),1)+...
                T(floor(Trace(i,1)),ceil(Trace(i,2)))*(1-mod(Trace(i,1),1))*mod(Trace(i,2),1);
        end
        % calculate the distance between the points on the trace
        TotalDist(1)=0;
        dist(1)=0;
        for i=1:length(Trace)
            if i<length(Trace)
            dist(i+1)=sqrt((Trace(i+1,1)-Trace(i,1))^2+(Trace(i+1,2)-Trace(i,2))^2);
            end
            if i==1
                TotalDist(i)=0;
            else
                TotalDist(i)=TotalDist(i-1)+dist(i);
            end
            hist(i)=intensity(i);
        end
        Kill=0;
        if firstframe<hist(1)-20
            hist(2:end+1)=hist;
            TotalDist(2:end+1)=TotalDist;
            Trace(2:end+1,1:2)=Trace;
            hist(1)=firstframe;
        end
        if lastframe>hist(end)+20
            hist(end+1)=lastframe;
            TotalDist(end+1)=TotalDist(end);
            Trace(end+1,1:2)=Trace(end,1:2);
            
        end
    else
        % values thrown out to indicate failed trace
        Trace=0;
        hist=0;
        TotalDist=0;
        Kill=1;
    end
    end
    