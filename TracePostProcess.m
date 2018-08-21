function [TraceOut, FinalSpeeds, FinalPosx, FinalPosy, medianlength, FinalLength, Lengths]=TracePostProcess(FilEnd, param, frames)

%% processes the generated trace by first cleaning the trace, then sectioning
% followed by calculation of velocities
for j=1:length(FilEnd)
     if FilEnd(j).reject==0
    [TraceOut(j).Final,TraceOut(j).FinalDist,TraceOut(j).FinalTrace]=...
        CleanTrace(FilEnd(j).Frames,FilEnd(j).Dist,FilEnd(j).Trace,FilEnd(j).historyshift,FilEnd(j).Coord,FilEnd(j).firstframe);
    
    TraceOut(j).FinalPoints=SectionTrace(TraceOut(j).Final, TraceOut(j).FinalDist, param);
     end
end
FinalSpeeds=zeros(length(FilEnd),frames);
for i=1:length(TraceOut)
    if length(TraceOut(i).Final)>=2  & FilEnd(i).reject==0
        try
            [TraceOut(i).FinalPoints, TraceOut(i).FinalSmooth, Speed,Posx,Posy, FilEnd(i).reject] = FindVelocity...
                (TraceOut(i).FinalPoints,TraceOut(i).Final, TraceOut(i).FinalDist, TraceOut(i).FinalTrace, param, FilEnd(i).lastframe, FilEnd(i).firstframe);
        catch
            TraceOut(i).FinalPoints=0; TraceOut(i).FinalSmooth=0; Speed=0;Posx=0;Posy=0;
             FilEnd(i).reject=3;
            fail=i
        end
        FinalSpeeds(i,1:length(Speed))=Speed;
        FinalPosx(i,1:length(Posx))=Posx;
        FinalPosy(i,1:length(Posy))=Posy;
        medianlength(i)=median(FilEnd(i).FilLength);
        
    end
    
end

for i=1:length(FinalSpeeds(:,1))
A=find(abs(FinalSpeeds(i,FilEnd(i).firstframe:FilEnd(i).lastframe))>0);
if A
FinalLength(i,FilEnd(i).firstframe+A-1)=FilEnd(i).FilLength(A);
else
    FinalLength(i,FilEnd(i).firstframe+A-1)=0;
end
end

for i=1:length(FinalLength(1,:))
Lengths(i)=mean(nonzeros(FinalLength(:,i)));
end
