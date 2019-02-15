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
FinalSpeeds=sparse(length(FilEnd),frames);
FinalPosx=sparse(length(FilEnd),frames);
FinalPosy=sparse(length(FilEnd),frames);
for i=1:length(TraceOut)
    if length(TraceOut(i).Final)>=4  & FilEnd(i).reject==0
        try

            [TraceOut(i).FinalPoints, TraceOut(i).FinalSmooth, Speed,Posx,Posy, FilEnd(i).reject] = FindVelocity...
                (TraceOut(i).FinalPoints,TraceOut(i).Final, TraceOut(i).FinalDist, TraceOut(i).FinalTrace, param, FilEnd(i).lastframe, FilEnd(i).firstframe);
        catch
            TraceOut(i).FinalPoints=0; TraceOut(i).FinalSmooth=0; Speed=0;Posx=0;Posy=0;
            FilEnd(i).reject=3;
            fail=i
        end
        FinalSpeeds(i,find(Speed))=nonzeros(Speed);
        FinalPosx(i,find(Posx))=nonzeros(Posx);
        FinalPosy(i,find(Posy))=nonzeros(Posy);
        medianlength(i)=median(nonzeros(FilEnd(i).FilLength(:,1)));
        
    end
    
end

for i=1:length(FinalSpeeds(:,1))
    FinalLength(i,FilEnd(i).FilLength(:,2))=FilEnd(i).FilLength(:,1);
    
end

for i=1:length(FinalLength(1,:))
    Lengths(i)=mean(nonzeros(FinalLength(:,i)));
end
