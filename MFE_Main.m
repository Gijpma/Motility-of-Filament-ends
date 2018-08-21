%% Motility of Filament Ends software
% user can select single or multiple files for analysis. parameter settings
% are set in MFE_Process.m

clearvars
close all

%% setting parameters

tic
%% end setting parameters

% Loading video file
[filelist,dirout]=uigetfile('*.avi','multiselect','on');
if iscell(filelist)
    if ~isempty(dirout)
        parfor VideoNow=1:length(filelist)
            tic
            file=filelist{VideoNow};
            try
                MFE_Process(dirout,file,VideoNow)
            catch
                ['error in ' file]
            end
            toc
        end
    end
else
    file=filelist;
    try
        MFE_Process(dirout,file,1)
    catch
        ['error in ' file]
    end
    
end
toc