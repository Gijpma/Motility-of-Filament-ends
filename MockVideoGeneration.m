clearvars;
close all;

% --- This part creates traces with bent random bending, then fits them to
% splines, which it then stores.



thickness=1;

height = 480;
width = 720;
resolutionfactor=10/(30*0.124);
vmax = [5 5 1]*resolutionfactor;%(1/(0.124))/30;% in microns per second, with 0.124 micron/px
vmin=[1 1 0.2]*resolutionfactor;%(1/(0.124))/30;
SpeedChangeStart=300;
SpeedChangeEnd=500;
fmotmax = [1 0.5 1];
fmotmin= [1 0.5 0.5];
fmotChangeStart=300;
fmotChangeEnd=500;
StopAndGoRate=0.003;
tracelength = 10000;
time_trace = 1:tracelength+2;

bend =[0.5 0.5 0.5]/3;
filamentcount= 25;
max_frame = 1000;
min_fil_length = [25 25 300];
max_fil_length = [100 100 300];
% fmot = 0.5; % between [0, 1.0], at 0.5, 50% of filaments starts moving
% fmot_change = 0; % between [-0.5, 0.5]. E.g. for 0.5 all filaments will move by the end
% fmot_change_rate = 0; % between [1,100]. 1 changes fastest, 5 changes slower, and so on
speckle = [0 0.005 0.005]; % speckle noise on the whole image (0.005 standard)
orientation = [1 1 1];
tic
for Complex=3
    for Regime=3
        foldername = 'D:\Gijs.Ijpma\My papers\In Progress 2018 - Filament Tracking\TestFiles\2\';
        output = ['Vmax', num2str(vmax(Regime)), '_Vmin', num2str(vmin(Regime)),'_bend', num2str(bend(Complex)),...
            '_fmotmax',num2str(fmotmax(Regime)),...
            '_fmotmin',num2str(fmotmin(Regime)),...
            '_speckle',num2str(speckle(Complex)),...
            '_orientation',num2str(orientation(Complex)),...
            'Fil_min',num2str(min_fil_length(Complex)),...
            'Fil_max',num2str(max_fil_length(Complex))];
%         output='Scenario4B_v2.avi';
        output
        for i=1:filamentcount
            xstart(i)=rand*width;
            ystart(i)=rand*height;
            traces(i,1,1)=xstart(i);
            traces(i,1,2)=ystart(i);
            theta1=orientation(Complex)*(rand*2*pi);
            traces(i,2,1)=xstart(i)+1.*cos(theta1);
            traces(i,2,2)=ystart(i)+1.*sin(theta1);
            
            for t = 1:tracelength
                theta2=atan2((traces(i,t+1,2)-traces(i,t+0,2)),(traces(i,t+1,1)-traces(i,t+0,1)));
                newtheta=theta2+(0.5-rand)*bend(Complex);
                traces(i,t+2,1)=traces(i,t+1,1)+.1.*cos(newtheta);
                traces(i,t+2,2)=traces(i,t+1,2)+.1.*sin(newtheta);
            end
            
            s1 = spline(time_trace,traces(i,:,1));
            s2 = spline(time_trace,traces(i,:,2));
            splines{1,i} = s1;
            splines{2,i} = s2;
            [num2str(i) ' out of ' num2str(filamentcount) ' traces is done'];
        end
        
        %%
        % --- Here we actually let filaments walk on these traces. Here we vary the
        % size and the speed of the filaments
  
        
        
        toc
        
        filament_length = 10*(rand(filamentcount,1)*(max_fil_length(Complex) - min_fil_length(Complex)) + min_fil_length(Complex));
        speed = zeros(filamentcount, max_frame);
        fmot = zeros(max_frame,1);
        n = 2;
        
        for i = 1:filamentcount
            time{i,1} = 1:filament_length(i,1);
            speed(i,1:SpeedChangeStart) = vmax(Regime);
            speed(i,SpeedChangeStart+1:SpeedChangeEnd) = vmax(Regime)-((vmax(Regime)-vmin(Regime))/...
                (SpeedChangeEnd-SpeedChangeStart))*[1:SpeedChangeEnd-SpeedChangeStart];
            speed(i,SpeedChangeEnd:end) = vmin(Regime);
        end
        fmot(1:fmotChangeStart) = fmotmax(Regime);
        fmot(fmotChangeStart+1:fmotChangeEnd) = fmotmax(Regime)-((fmotmax(Regime)-fmotmin(Regime))/...
            (fmotChangeEnd-fmotChangeStart))*[1:fmotChangeEnd-fmotChangeStart];
        fmot(fmotChangeEnd:end) = fmotmin(Regime);
        
        ismoving(:,1) = round(rand(filamentcount,1) + fmotmax(Regime) - 0.5);
        ismoving(:,2) = 0;
        
        aviobj = VideoWriter([foldername output '.avi']);
        open(aviobj)
        M.colormap = [];
        Mwrite.colormap = [];
        
        frameoutlarge=uint8(zeros(10*480,10*720));
        FlipRemain=0;
        for tt = 1:max_frame
           
            ismoving(:,2)=ismoving(:,2)+1;
            CurrentStopped=filamentcount-sum(ismoving(:,1));
            StoppedFil=round(filamentcount*(1-fmot(tt)));
            if tt<100
                StopAndGo=0;
            else
                StopAndGo=StopAndGoRate;
            end
            FlipRaw=filamentcount*(1-fmot(tt))*StopAndGo+FlipRemain;
            FlipNumber=floor(FlipRaw);
            FlipRemain=FlipRaw-FlipNumber;
            ToStop=round((StoppedFil-CurrentStopped));
            if ToStop<0 & FlipNumber>abs(ToStop)
                ToGo=FlipNumber;
                ToStop=FlipNumber-ToStop;
            elseif ToStop<0
                ToGo=abs(ToStop);
                ToStop=0;
            elseif ToStop>=0 & FlipNumber>abs(ToStop)
                ToGo=FlipNumber-ToStop;
                ToStop=FlipNumber;
            else
                ToGo=0;
            end
            [~,sortedstop]=sort((1-ismoving(:,1)).*ismoving(:,2),'descend');
            if ToGo>0
                ismoving(sortedstop(1:ToGo),2)=0;
                ismoving(sortedstop(find(ismoving(:,1)==0,ToGo)),1)=1;
            end
            [~,sortedgo]=sort(ismoving(:,1).*ismoving(:,2),'descend');
            if ToStop>0
                ismoving(sortedgo(find(ismoving(:,1)==1,ToStop)),2)=0;
                ismoving(sortedgo(find(ismoving(:,1)==1,ToStop)),1)=0;
            end
            
            frameoutlarge(:,:)=0;
            
            for i = 1:filamentcount
                
                x=round(10*(ppval(splines{1,i},time{i,tt}(1):0.05:time{i,tt}(end))));
                y=round(10*(ppval(splines{2,i},time{i,tt}(1):0.05:time{i,tt}(end))));
                use=sub2ind(size(frameoutlarge),mod(y,10*480)+1,mod(x,10*720)+1);
                xout(i,tt)=(ppval(splines{1,i},time{i,tt}(1)));
                frameoutlarge(use)=255;
                speed(i,tt) = round(ismoving(i,1)) * speed(i,tt);
                time{i,tt+1} = time{i,tt} + speed(i,tt);
                
            end
            
            
            frameout=10*imresize(frameoutlarge,[480 720]);
            
            
            frameout=imgaussfilt(frameout,1.2)*1.2;
            
            frameout = imnoise(frameout, 'gaussian',speckle(Complex),speckle(Complex));
            
            writeVideo(aviobj,frameout);
%             imwrite(frameout, [foldername output '.tiff'], 'writemode', 'append')
            moving(:,tt)=ismoving(:,1);
        end
        
        
        % close
        close(aviobj);
        
        fprintf('DONE\n')
        toc
    end
end 
    
    
