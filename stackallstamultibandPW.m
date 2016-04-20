%%
if 1
    clear, close all
    % now load all those daylongvariables files and put it together
    d=dir('daylongxcvariables*.mat');
    xcsumall=[];
    timeall =[];
    for nn=1:length(d)
        load(d(nn).name);
        xcsumall=[xcsumall, xcsum];
        timeall =[timeall, time];
        clear xcsum time
    end
    %
    cutoffxc=2.;
    [~,b]=find(xcsumall>cutoffxc);
    [~,peaksind]=findpeaks(xcsumall,'MINPEAKHEIGHT',cutoffxc);
    disp([num2str(length(peaksind)),' events at ',num2str(cutoffxc)]);
    
    % Try cutoffxc=2.0;
    %stationF9C doesn't start until the 17th
    peaksind(1:43)=[];
    % station SW died on the 21st
    peaksind(23:end)=[];
    % station F9A has NaNs for the event on 21-Jan-2009 06:06:22
    peaksind(22)=[];
    
    
    % station F9NE has really large Z on event 6 (20-Jan-2009 03:16:28)
    % if I leave this in, then it throws off the other indices:
    peaksind(6)=[];
    
    % QC clear out some others
    
    % if the event 6 from above is removed, use this
    peaksind([1:3,6,10,14,16:18])=[]; %otherwise:
    %     peaksind([1:3,7,11,15,17:19])=[];
    peaksind(2)=[];
    
    %htese are very low amplitude
    % if the event 6 from above is F9NE removed, use this
    peaksind([4,8])=[];% otherwise:
    %     peaksind([4,9])=[];
    
    
    % station F9SW has opposite polarity for this one (21-Jan-2009 03:52:58)
    peaksind(7:8)=[];
    % peaksind(8:9)=[];
    % peaksind(6)=[];
    
    % finally, this event has the wrong polarity on F9A Z for BP400-60, but I
    % think we should keep it in the mix peaksind(5)=[];
    
    
    
    %     % Use these for cutoffxc=2.3
    %         %stationF9C doesn't start until the 17th
    %     peaksind(1:33)=[];
    %         % station SW died on the 21st
    %     peaksind(17:39)=[];
    %
    %     % station F9A has NaNs for the event on 21-Jan-2009 06:06:22
    %     peaksind(16)=[];
    
    disp([num2str(length(peaksind)),' events on all stations']);
    
    % Use these for cutoffxc=2.5
    % station SW
    %     peaksind(31:45)=[];
    %     % %stationNW is missing some as well
    %     peaksind(10:22)=[];
    %     %stationF9C
    %     peaksind(1:9)=[];
    
    % % Use these for cutoffxc=2.6 %6 events
    % % station SW
    % peaksind(24:36)=[];
    % %stationNW is missing some as well
    % peaksind(8:17)=[];
    % %stationF9C
    % peaksind(1:7)=[];
    
    % Use these for cutoffxc=2.7 % 4 events
    % peaksind(17)=[];
    % peaksind(1:8)=[];
    % %station SW is missing a lot of days
    % peaksind(9:17)=[];
    % %
    % % %stationNW is missing some as well
    % peaksind(1:4)=[];
    %
    figure
    plot(timeall,xcsumall)
    hold on;
    plot(timeall(peaksind),xcsumall(peaksind),'or')
    
    datetick('x')
    
    
    %
    % close all
    
    
    rootdir='/Users/gpwaite/Data/Guate/Fuego2012/FuegoVLP2012';
    
    %make sure the data are consistent with padded GFs
    td=2^16; % need 655.36 seconds to get the tilt
    % fd=td/2;
    
    hcall=[1/60 1/10  1/10  1/10 1/10 1/10];
    lcall=[1/400 1/400 1/120 1/60 1/30 1/90];
    skipmat=zeros(20,11);
    
    for lcn=6%1:length(lcall);
        close all
        lc=lcall(lcn);
        hc=hcall(lcn);
        % load the data
        
        % Here we will get the time domain seismograms from the mat files using
        % GISMO. This should allow for easier testing of many different events.
        % staall=[{'SW1'};{'SW2'};{'S1'};{'SE1'};{'SE2'};{'N1'};{'NE1'};{'NE2'};{'NW1'};{'NW2'};];
        
        %These are espc: F9SW, F9NW, F904, F9B, F9NE,
        % F9E and F9D have bad tilt signals
        %   staall=[{'F9SW'};{'F9NW'};{'F900'};{'F9E'};{'F9D'};{'F9C'};{'F9B'};{'F9A'};{'F9NE'};{'F9SE'};];
        staall=[{'F9SW'};{'F9NW'};{'F900'};{'F9C'};{'F9B'};{'F9A'};{'F9NE'};{'F9SE'};];
        %     staall={'F900'}
        for stasub=1:length(staall);
            
            sta=staall(stasub);
            % this is the order that the data are in in the GFs
            allchn=[{'HHE'};{'HHN'};{'HHZ'}];
            inst=[{'F9_espc'};{'F9_espc'};{'F9_cmg40t'};{'F9_cmg40t'};{'F9_espc'};{'F9_cmg40t'};{'F9_espc'};{'F9_cmg40t'};];
            instsub=inst(stasub);
            net='XT';
            loc='01';
            
            stn=char(sta);
            
            
            loadsource=datasource('file','/Users/gpwaite/Data/Guate/Fuego2009/Fuego_2009_mat/%s/%s%s%04d%02d%02d.mat','station','station','channel','year','month','day');
            % now search through all the data and find events to stack
            scnle = scnlobject(stn,char(allchn(1)),net,loc);
            scnln = scnlobject(stn,char(allchn(2)),net,loc);
            scnlz = scnlobject(stn,char(allchn(3)),net,loc);
            
            if strcmp(char(instsub),'F9_cmg40t');
                
                fo=filterobject('B',[lc, hc],2);
                %             fot=filterobject('B',[1/500 1/80],2);
            elseif strcmp(char(instsub),'F9_espc');
                
                fo=filterobject('B',[lc, hc],2);
                %             fot=filterobject('B',[1/500 1/80],2);
            else
                error('wrong instrument');
                
            end
            
            tmpVLPE=zeros(td,length(peaksind));
            tmpVLPN=tmpVLPE;
            tmpVLPZ=tmpVLPE;
            %         tmpULPE=tmpVLPE;
            %         tmpULPN=tmpVLPE;
            %         tmpULPZ=tmpVLPE;
            
            % used 3000 seconds for td=2^15
            % use 3164 seconds for td=2^16
            prewin=3164; %seconds
            %         prewin=3328
            %         prewin=3000
            prewinsamp=prewin*100 -3250-16384; % an extra 164seconds
            
            for n=1:length(peaksind)
                tmpstartday=timeall(peaksind(n))-prewin/86400; %tmpstartm is in days
                tmpstopday=tmpstartday+.75;
                
                disp(['working on ',datestr(tmpstartday)]);
                
                %grab a whole day, or this will be a mess
                woe=waveform(loadsource,scnle,tmpstartday,tmpstopday);
                sps=get(woe,'freq');
                if ~isempty(woe)
                    
                    woe=set(woe,'data',double(woe)-mean(double(woe)));
                    % get the VLP
                    wof=filtfilt(fo,woe);
                    tmpVLP=double(wof);
                    tmpVLPE(:,n)=tmpVLP(prewinsamp+1:prewinsamp+td);
                    % and the ULP
                    %                 woft=filtfilt(fot,woe);
                    %                 tmpULP=double(woft);
                    %                 tmpULPE(:,n)=tmpULP(prewinsamp+1:prewinsamp+td);
                    clear tmpVLP
                    
                    won=waveform(loadsource,scnln,tmpstartday,tmpstopday);
                    won=set(won,'data',double(won)-mean(double(won)));
                    wof=filtfilt(fo,won);
                    tmpVLP=double(wof);
                    tmpVLPN(:,n)=tmpVLP(prewinsamp+1:prewinsamp+td);
                    % and the ULP
                    %                 woft=filtfilt(fot,won);
                    %                 tmpULP=double(woft);
                    %                 tmpULPN(:,n)=tmpULP(prewinsamp+1:prewinsamp+td);
                    clear tmpVLP
                    
                    woz=waveform(loadsource,scnlz,tmpstartday,tmpstopday);
                    woz=set(woz,'data',double(woz)-mean(double(woz)));
                    wof=filtfilt(fo,woz);
                    tmpVLP=double(wof);
                    tmpVLPZ(:,n)=tmpVLP(prewinsamp+1:prewinsamp+td);
                    % and the ULP
                    %                 woft=filtfilt(fot,woz);
                    %                 tmpULP=double(woft);
                    %                 tmpULPZ(:,n)=tmpULP(prewinsamp+1:prewinsamp+td);
                    clear tmpVLP
                    
                else
                    disp(['skipping this one ',datestr(tmpstartday),',  number ',num2str(n)])
                    tmpVLPZ(:,n)=zeros(td,1);
                    tmpVLPN(:,n)=zeros(td,1);
                    tmpVLPE(:,n)=zeros(td,1);
                    sps=100;
                    skipmat(n,stasub)=1;
                    error(['skip this one ',datestr(timeall(peaksind(n))),',  number ',num2str(n)])
                end
                %     imagesc(tmpVLPZ'); drawnow
                %     pause
            end
            
            
            figure
            suptitle2([stn,': BP ',num2str(1/lc),'-',num2str(1/hc)])
            subplot(311)
            imagesc((1:td)/sps,1:length(peaksind),tmpVLPE')
            ylabel('east')
            subplot(312)
            imagesc((1:td)/sps,1:length(peaksind),tmpVLPN')
            ylabel('north')
            subplot(313)
            imagesc((1:td)/sps,1:length(peaksind),tmpVLPZ')
            ylabel('vert')
            xlabel('time')
            figfile=[stn,'_BP',num2str(1/lc),'-',num2str(1/hc),'.ps'];
            print('-dpsc',figfile)
            
            drawnow
            
            %         if lc<1/500
            %             figure
            %             suptitle2([stn,': BP 500-80'])
            %
            %             subplot(311)
            %             imagesc((1:td)/sps,1:length(peaksind),tmpULPE')
            %             subplot(312)
            %             imagesc((1:td)/sps,1:length(peaksind),tmpULPN')
            %             subplot(313)
            %             imagesc((1:td)/sps,1:length(peaksind),tmpULPZ')
            %
            %         end
            stafile=[stn,'_BP',num2str(1/lc),'-',num2str(1/hc),'.mat'];
            save(stafile,'tmpVLPE','tmpVLPN','tmpVLPZ','lc','hc');
            %             pause
        end
        %     clear  tmpVLPE tmpVLPN tmpVLPZ
    end
    
end
%% stack
% finally ready to stack
% use F900 HHZ to normalize everything
clear

dt50=1/50;
sps=50;

td=2^16;
tddec=td/2;
fd=tddec;
taperwin=tukeywin(td,.5); %

hcall=[1/60  1/10   1/10   1/10  1/10  1/10];
lcall=[1/400 1/400, 1/120, 1/60, 1/30, 1/90];
for mmm=6%1:length(lcall)
    lc=lcall(mmm);
    hc=hcall(mmm);
    
    % get noise subset for degrees of freedom
    noiselen=100*sps;
    if lc<1/400
        noisestart=1;
    else
        noisestart=sps*101;
    end
    noiseend=noisestart+noiselen;
    
    %
    
    load(sprintf('F900_BP%d-%d.mat',1/lc,1/hc))
    normvec=range(tmpVLPZ);
    normfrac=normvec/max(normvec);
    
    % now load and stack each of the other stations
    staall=[{'F9SW'};{'F9NW'};{'F900'};{'F9C'};{'F9B'};{'F9A'};{'F9NE'};{'F9SE'};];
    staall=[{'F9SW'};{'F9NW'};{'F900'};{'F9C'};{'F9B'};{'F9A'};{'F9NE'};];
    
    %     staall=[{'F9SW'};{'F9NW'};{'F900'};{'F9C'};{'F9B'};{'F9A'}];
    %     goodind(1,:)=[ones(1,3),zeros(1,9)]; %F9SW
    %     goodind(2,:)=[ones(1,5),zeros(1,3),ones(1,4)]; %F9NW
    %     goodind(3,:)=ones(1,12); %F900
    %     goodind(4,:)=[ones(1,7),zeros(1,3),ones(1,2)]; %F9C
    %     goodind(5,:)=[ones(1,7),zeros(1,3),ones(1,2)]; %F9B
    %     goodind(6,:)=[ones(1,3),zeros(1,5),ones(1,2),zeros(1,2)]; %F9A
    %     goodind(7,:)=ones(1,12); %F9NE
    %     goodind(8,:)=[ones(1,8),zeros(1,4)]; %F9SE
    
    
    
    nt=1;
    for stasub=1:length(staall);
        datatd(nt,:)=zeros(1,tddec);
        datatd(nt+1,:)=zeros(1,tddec);
        datatd(nt+2,:)=zeros(1,tddec);
        
        sta=staall(stasub);
        stn=char(sta);
        disp([' working on ',stn])
        load(sprintf('%s_BP%d-%d.mat',stn,1/lc,1/hc));
        %         tmpind=logical(goodind(stasub,:));
        tmpind=1:size(tmpVLPE,2);
        tmpE=tmpVLPE(:,tmpind);
        tmpN=tmpVLPN(:,tmpind);
        tmpZ=tmpVLPZ(:,tmpind);
        
        if strcmp('F9B',stn)
            disp('swapping horizontal channels for F9B')
            tmpE=tmpVLPN(:,tmpind);
            tmpN=tmpVLPE(:,tmpind);
            
        end
        if strcmp('F9NE',stn)
            disp('removing the vertical on F9NE')
            tmpZ=tmpVLPZ(:,tmpind)*0;
            
        end
        
        
        %     normvecZ=range(tmpVLPZ(:,tmpind));
        %     % compute the fraction of the sum of each
        %     normvecZsum=sum(normvecZ);
        %     normfracZ=normvecZ/normvecZsum;
        % normalize each according to the size of F900HHZ
        % now normalize to the whole sequence from F900
        normfracsub=normfrac(tmpind);
        
        for n=1:size(tmpZ,2)
            tmpEnew(:,n)=tmpE(:,n)/normfracsub(n);
            tmpNnew(:,n)=tmpN(:,n)/normfracsub(n);
            tmpZnew(:,n)=tmpZ(:,n)/normfracsub(n);
            
            tmp=decimate(tmpE(:,n),2);
            tmpE50(:,n)=tmp/normfracsub(n); clear tmp
            tmp=decimate(tmpN(:,n),2);
            tmpN50(:,n)=tmp/normfracsub(n); clear tmp
            tmp=decimate(tmpZ(:,n),2);
            tmpZ50(:,n)=tmp/normfracsub(n); clear tmp
        end
        stackE=sum(tmpEnew,2)/size(tmpE,2);
        stackN=sum(tmpNnew,2)/size(tmpN,2);
        stackZ=sum(tmpZnew,2)/size(tmpZ,2);
        
        
        % have not been deconvolved
        stackE=stackE.*taperwin;
        stackN=stackN.*taperwin;
        stackZ=stackZ.*taperwin;
        
        stackE50=decimate(stackE,2);
        stackN50=decimate(stackN,2);
        stackZ50=decimate(stackZ,2);
        
        
        %         % Here is where we use phase weighted stack
        %                 [pwsE,tvecE] = pwstack(tmpE50',size(tmpE50,2),dt50);
        %                 [pwsN,tvecN] = pwstack(tmpN50',size(tmpN50,2),dt50);
        %                 [pwsZ,tvecZ] = pwstack(tmpZ50',size(tmpZ50,2),dt50);
        [pwsE] = pws_simple(tmpE50',dt50,round((10)/dt50));
        [pwsN] = pws_simple(tmpN50',dt50,round((10)/dt50));
        [pwsZ] = pws_simple(tmpZ50',dt50,round((10)/dt50));
        
        % taper these babies because the ends are not dealt with properly
        % by the tukey window smoothing
        tkw=tukeywin(length(pwsE),.5); %with .5 its a hanning window
        pwsE=tkw.*pwsE';
        pwsN=tkw.*pwsN';
        pwsZ=tkw.*pwsZ';
        
        datatd(nt,:) =stackE50';
        datatdpws(nt,:)=pwsE;
        %         DoF(nt,:)=dof(pwsE(noisestart:noiseend))*fd/noiselen;
        % FFT
        tmpfft = fft(datatdpws(nt,:));
        % put in matrix
        data(nt,:)=tmpfft(1:fd); clear tmpfft
        %NORTH
        datatd(nt+1,:) =stackN50';
        datatdpws(nt+1,:)=pwsN;
        %         DoF(nt+1,:)=dof(pwsN(noisestart:noiseend))*fd/noiselen;
        
        % FFT
        tmpfft = fft(datatdpws(nt+1,:));
        % put in matrix
        data(nt+1,:)=tmpfft(1:fd); clear tmpfft
        %VERT
        datatd(nt+2,:)= stackZ50';
        datatdpws(nt+2,:)=pwsZ;
        %         DoF(nt+2,:)=dof(pwsZ(noisestart:noiseend))*fd/noiselen;
        
        
        % FFT
        tmpfft = fft(datatdpws(nt+2,:));
        % put in matrix
        data(nt+2,:)=tmpfft(1:fd); clear tmpfft
        nt=nt+3;
        
        figure
        suptitle2([char(sta),' BP ',num2str(1/lc), ' - ',num2str(1/hc)])
        tvec=(0:tddec-1)*dt50;
        subplot(311)
        plot(tvec,stackE50)
        hold on
        subplot(311)
        ylabel('E')
        plot(tvec,pwsE)
        
        subplot(312)
        plot(tvec,stackN50)
        hold on
        subplot(312)
        ylabel('N')
        plot(tvec,pwsN)
        
        subplot(313)
        plot(tvec,stackZ50)
        hold on
        subplot(313)
        plot(tvec,pwsZ)
        ylabel('Z')
        drawnow
        legend('TD stack','PWS')
        
        
    end
    %     DoF
    savefile=sprintf('PWstack2009BP%d_%d.mat',1/lc,1/hc);
    save(savefile,'data','datatd','datatdpws','staall','sps')
    %     figure
    %     plot(datatd'); title(savefile)
end
%%

% Adding a new section that pltos the ffts for each band
clear, close all
staall=[{'F9SW'};{'F9NW'};{'F900'};{'F9C'};{'F9B'};{'F9A'};{'F9NE'};];

datdir='/Users/gpwaite/Data/Guate/Fuego2009';
figdir='~/Documents/Papers/Fuego-ULP-VLP';
dt50=1/50;
sps=50;

td=2^16;
fd=td/2;
fvec=linspace(0,sps,fd);

hcall=[1/60  1/10   1/10   1/10  1/10  1/10];
lcall=[1/400 1/400, 1/120, 1/60, 1/30, 1/90];
for mmm=1:length(lcall)
    lc=lcall(mmm);
    hc=hcall(mmm);
    
    loadfile=sprintf('PWstack2009BP%d_%d.mat',1/lc,1/hc);
    load([datdir,filesep,loadfile]);
    nt=1;
    for stasub=1:length(staall);
        
        mxabs(mmm,stasub)=max(max(abs(data(nt:nt+2,:))));
        nt=nt+3;
        
    end
    
end
maxmxabs=max(mxabs);
%
for mmm=1:length(lcall)
    lc=lcall(mmm);
    hc=hcall(mmm);
    
    loadfile=sprintf('PWstack2009BP%d_%d.mat',1/lc,1/hc);
    load([datdir,filesep,loadfile]);
    
    nt=1;
    for stasub=1:length(staall);
        sta=staall(stasub);
        stn=char(sta);
        disp([' working on ',stn])
        
        
        
        figure(stasub)
        %         suptitle2([char(sta),' BP ',num2str(1/lc), ' - ',num2str(1/hc)])
        subplot(311)
        semilogx(fvec,abs(data(nt,:))/maxmxabs(stasub))
        hold on
        ylabel('E')
        xlim([0 0.2])
        ylim([0 1])
        
        subplot(312)
        semilogx(fvec,abs(data(nt+1,:))/maxmxabs(stasub))
        hold on
        ylabel('N')
        xlim([0 0.2])
        ylim([0 1])
        
        
        subplot(313)
        semilogx(fvec,abs(data(nt+2,:))/maxmxabs(stasub))
        hold on
        ylabel('Z')
        xlim([0 0.2])
        ylim([0 1])
        
        % add some text to the last plot
        if mmm==length(lcall)
            suptitle2(stn);
            legend(num2str(1./lcall'))
            drawnow
            print('-dpsc',[figdir,filesep,'station',stn,'_spectra_PWS_6bands.ps'])
            
        end
        
        
        
        if strcmp(stn,'F900')
            figure(99)
            subplot(3,3,1)
            semilogx(fvec,abs(data(nt,:))/maxmxabs(stasub))
            hold on
            ylabel(stn);
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            subplot(3,3,2)
            semilogx(fvec,abs(data(nt+1,:))/maxmxabs(stasub))
            hold on
            %             ylabel('N')
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            
            subplot(3,3,3)
            semilogx(fvec,abs(data(nt+2,:))/maxmxabs(stasub))
            hold on
            %             ylabel('Z')
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            
        elseif strcmp(stn,'F9NW')
            figure(99)
            subplot(3,3,4)
            semilogx(fvec,abs(data(nt,:))/maxmxabs(stasub))
            hold on
            ylabel(stn);
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            subplot(3,3,5)
            semilogx(fvec,abs(data(nt+1,:))/maxmxabs(stasub))
            hold on
            %             ylabel('N')
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            
            subplot(3,3,6)
            semilogx(fvec,abs(data(nt+2,:))/maxmxabs(stasub))
            hold on
            %             ylabel('Z')
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            
            %         elseif strcmp(stn,'F9B')
            %             figure(99)
            %             subplot(4,3,7)
            %             semilogx(fvec,abs(data(nt,:))/maxmxabs(stasub))
            %             hold on
            %             ylabel(stn);
            %             xlim([0 1])
            %
            %             subplot(4,3,8)
            %             semilogx(fvec,abs(data(nt+1,:))/maxmxabs(stasub))
            %             hold on
            %             %             ylabel('N')
            %             xlim([0 1])
            %
            %
            %             subplot(4,3,9)
            %             semilogx(fvec,abs(data(nt+2,:))/maxmxabs(stasub))
            %             hold on
            %             %             ylabel('Z')
            %             xlim([0 1])
            %
            
        elseif strcmp(stn,'F9SW')
            figure(99)
            subplot(3,3,7)
            semilogx(fvec,abs(data(nt,:))/maxmxabs(stasub))
            hold on
            ylabel(stn);
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            subplot(3,3,8)
            semilogx(fvec,abs(data(nt+1,:))/maxmxabs(stasub))
            hold on
            %             ylabel('N')
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            
            subplot(3,3,9)
            semilogx(fvec,abs(data(nt+2,:))/maxmxabs(stasub))
            hold on
            %             ylabel('Z')
            xlim([0 0.2])
            ylim([0 1])
            set(gca,'YTickLabel','');
            set(gca,'Xtick',[0.01 0.1 1]);
            tmpx=get(gca,'XTick');
            set(gca,'XtickLabel',num2str(tmpx'))
            
            
            
            % add some text to the last plot
            if mmm==length(lcall)
                legend(num2str(1./lcall'))
                drawnow
                print('-dpsc',[figdir,filesep,'example_spectra_PWS_6bands.ps'])
                
            end
            
            
        end
        nt=nt+3;
    end
end