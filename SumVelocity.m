function [Velocity, fmot, fnumber]=SumVelocity(Speeds)
for i=1:length(Speeds(1,:))
    Velocity(i)=mean(Speeds(find(Speeds(:,i)>0),i));
    if isnan(Velocity(i))
        Velocity(i)=0;
    end
    fmot(i)=1-length(find(Speeds(:,i)<0))/(length(find(Speeds(:,i)>0 | Speeds(:,i)<0)));
    fnumber(i)=length(find(Speeds(:,i)>0));
end