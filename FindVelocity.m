function [Points, Smooth,Speed,Posx,Posy,fail] = FindVelocity(Points,Frames, Dist, Trace, param, lastframe, firstframe)
%% finds velocity from cleaned and sectioned traces.
Posx=0;
Posy=0;
CutOffVelocity=param.CutOffVelocity/param.VelocityFactor;
% calculates the slope of sections of the curve
for j=1:length(Points(:,1))-1
    Points(j,4)=(Points(j+1,2)-Points(j,2))...
        /(Points(j+1,1)-Points(j,1));
    Points(j,7)=(Points(j+1,1)-Points(j,1));
    if Points(j,1)==0
        Points(j,4)=0;
    end
end


if length(Points(:,1))>1
    cutoff=find(Points(:,4)>50,1);
    if cutoff
        %                 Points(cutoff+1:end,:)=[];
        Points(cutoff,:)=[];
        Points(cutoff,4)=0;
    end
    
    if range(Points(:,3))>0
        Points(:,5)=(Points(:,4)>CutOffVelocity);                     %moving sections
        Points(find(Points(:,7)<(1/CutOffVelocity) & Points(:,5)==0),5)=-1;   %sections too short to distinguish movement
        Points(find(Points(:,1)==0),5)=-1;
        % merge two adjacent stopped or go regions
        doubles=find(diff(Points(1:end-1,5))==0);
        
        while doubles
            %             
            Points(doubles+1,:)=[];
            for j=1:length(Points(:,1))-1
                Points(j,4)=(Points(j+1,2)-Points(j,2))...
                    /(Points(j+1,1)-Points(j,1));
            end
            stops=find(Points(1:end-1,4)<CutOffVelocity);
            
            doubles=find(diff(Points(1:end-1,5))==0 & Points(1:end-2,1)>0);
        end
        
        % Trace Smoothing
        %% testing smoothing prior to interpolation
        
        Go=find(Points(1:end-1,5)==1);
        if Go
            for Goi=1:length(Go)
                
                from=(find(Frames==Points(Go(Goi),1),1)+1);
                to=(find(Frames==Points(Go(Goi)+1,1),1,'last')-1);
                
                Points(Go(Goi),1)=Frames(from);
                Points(Go(Goi)+1,1)=Frames(to);
                Frames(from:to)=smooth(Frames(from:to),4);
            end
        end
        
        
        %% end testing smoothing prior to interpolation
        start=Points(1,1);
        Go=find(Points(1:end-1,5)==1);
        Stop=find(Points(1:end-1,5)==0);
        Points(1:end-1,6)=Points(2:end,1)-Points(1:end-1,1);
        
        
        % resampling of data to get one point per frame in moving sections
        
        
        clear Speed
        Speed=[];
         Smooth=[Frames(1) Dist(1)];
        if Go
            for Goi=1:length(Go)
                if Points(Go(Goi),6)>param.MinPresence 
                    from=ceil(Points(Go(Goi),3))+1;
                    to=floor(Points(Go(Goi)+1,3))-1;
                    %%
                    Pos1=Trace(from:to,1);
                    Pos2=Trace(from:to,2);
                    
                    x=[floor(Frames(from)):min(ceil(Frames(to)),max(floor(Frames)))];
                    
                    y=interp1(Frames(from:to),Dist(from:to),[floor(Frames(from)):min(ceil(Frames(to)),max(floor(Frames)))],'pchip');
                    Posx([floor(Frames(from)):min(ceil(Frames(to)),max(floor(Frames)))])=interp1(Frames(from:to),Pos1(1:end),[floor(Frames(from)):min(ceil(Frames(to)),max(floor(Frames)))],'pchip');
                    Posy([floor(Frames(from)):min(ceil(Frames(to)),max(floor(Frames)))])=interp1(Frames(from:to),Pos2(1:end),[floor(Frames(from)):min(ceil(Frames(to)),max(floor(Frames)))],'pchip');
                    
                    %%
                    if length(y) <= param.smooth
                        Speed(x)=mean(param.VelocityFactor*(diff(y))); %% fix the smoothing !!
                        insmoothx=[Smooth(:,1); x(1:end)'];
                        insmoothy=[Smooth(:,2); smooth(y,2)];
                    else
                        
                        if (y(param.smooth)-y(1))>0
                            a=polyfit([1:param.smooth],y(1:param.smooth),1);
                            y(1:param.smooth)=a(2)+a(1)*[1:param.smooth];
                        else
                            y(1:param.smooth)=0;       % unlikely to occur, but would cause error
                            % otherwise
                        end
                        if (y(end)-y(end-param.smooth+1))>0
                            a=polyfit([1:param.smooth],y(end-param.smooth+1:end),1);
                            y(end-param.smooth+1:end)=a(2)+a(1)*[1:param.smooth];
                           
                        else
                            y(1:param.smooth)=0;
                        end
                        
                        Speed(x(1:end-1))=param.VelocityFactor*(diff(smooth(y,param.smooth)));
                        
                        insmoothx=[Smooth(:,1); x'];
                        insmoothy=[Smooth(:,2); smooth(y,param.smooth)];
                    end
                    Smooth=[insmoothx insmoothy];
                    
                else
                    from=ceil(Points(Go(Goi),3))+1;
                    to=floor(Points(Go(Goi)+1,3))-1;
                    x=[floor(Frames(from)):min(ceil(Frames(to)),max(floor(Frames)))];
                    Speed(x(1:end-1))=0;
                end
            end
        end
        
        if Stop
            
            for stopi=1:length(Stop)
                from=round(Points(Stop(stopi),1));
                to=round(Points(Stop(stopi)+1,1));
                if ~from==0 & ~to==0
                    
                Speed(from:to)=-10;
                Posx(from:to)=Trace(Points(Stop(stopi),3),1);
                Posy(from:to)=Trace(Points(Stop(stopi),3),2);
                end
            end
        end
        
        if ~isempty(Speed)
            if Speed(end)<lastframe 
                if 1/(lastframe-length(Speed))<CutOffVelocity
                Posx(end+1:lastframe)=Posx(end);
                Posy(end+1:lastframe)=Posy(end);
                Speed(end+1:lastframe)=-10;
                end
            end
            if find(Speed,1)>firstframe 
                if 1/(find(Speed,1)-firstframe)<CutOffVelocity
                Posx(firstframe:find(Speed,1)-1)=Posx(find(Speed,1));
                Posy(firstframe:find(Speed,1)-1)=Posy(find(Speed,1));
                Speed(firstframe:find(Speed,1)-1)=-10;
                end
            end
            fail=0;
        else
            fail=4;
        end
        Smooth(end+1,:)=[Frames(end) Dist(end)];
        
    else
        fail=5;
    end
else
    fail=6;
end
