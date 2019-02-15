function [Trace, TraceDist,Position] = CleanTrace(Trace,TraceDist,PosTrace,Shift,Coord,Start)
% cleans the trace by removing trace errors
Presence=length(find(Coord(:,1)));

Position(:,1)=PosTrace(:,1)+Shift(2);
Position(:,2)=PosTrace(:,2)+Shift(1);
if Trace(1)>min(Trace(2:end))
    TraceDist(1)=[];
    Position(1,:)=[];
    Trace(1)=[];
end
% remove points that show negative frame movement
removepoint=[];
for i=1:length(Trace)-1
    if max(Trace(1:i))>Trace(i+1)+0.001
        removepoint(end+1)=i+1;
        
    end
end

if removepoint
    if length(removepoint)>max([0.2*length(Trace) 10])  % indicates bad trace
        TraceDist=[0 0];
        Position=[0 0];
        Trace=[0 0];
    else
        TraceDist(removepoint)=[];              
            Position(removepoint,:)=[];
            Trace(removepoint)=[];
    end
end


% remove vertical parts at start and end of trace
removefront=find(abs(Trace(2:end)-Trace(1))<0.0001,1,'last');

if removefront
    Trace(1:removefront)=[];
    TraceDist(1:removefront)=[];
    Position(1:removefront,:)=[];
end
removeback=find(abs(Trace(1:end-1)-Trace(end))<0.0001,1)+1;
if removeback
    Trace(removeback:end)=[];
    TraceDist(removeback:end)=[];
    Position(removeback:end,:)=[];
end

%remove vertical segments
removevertical=find(abs(diff(Trace))<0.001 & Trace(1:end-1)>0);
Trace(removevertical+1)=[];
TraceDist(removevertical+1)=[];
Position(removevertical+1,:)=[];

if Trace
    % remove jumps
    Flaws=find(diff(TraceDist)>6);
    if ~isempty(Flaws)
        for i=1:length(Flaws)
            TraceDist(Flaws(i)+1:end)=TraceDist(Flaws(i)+1:end)-(TraceDist(Flaws(i)+1)-TraceDist(Flaws(i)));
        end
    end
end

if Trace
    if Trace(1)>Start+50
        Trace=[Start Trace];
        TraceDist=[TraceDist(1) TraceDist];
        Position=[Position(1,:); Position];
        
    end
    TraceDist=TraceDist-TraceDist(1);
end

