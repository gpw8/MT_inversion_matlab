% Invert VLP data from either 2009 or 2012
% 25 sps
% This version properly accounts for tilt 11 April 2014
% it is written for GFs with 60 channels. The first 30 channels are for the
% 10 stations in 2009 and the next 30 channels are for the 10 stations from
% 2012.
%    X     Y    name
% 3360	3600	F9SW  1
% 4560	5400	F9NW  2
% 4920	5000	F900  3
% 4960	5600	F9E   4
% 5040	5680	F9D   5
% 5120	5720	F9C   6
% 5120	5640	F9B   7
% 5120	5520	F9A   8
% 5480	5320	F9NE  9
% 9640	2280	F9SE 10
% 3360	3560	SW1  11
% 3040	3280	SW2  12
% 3640	1560	S1   13
% 6440	3680	SE1  14
% 6840	3400	SE2  15
% 5000	5280	N    16
% 5040	5040	NE1  17
% 5480	5280	NE2  18
% 4920	5040	NW1  19
% 4360	5400	NW2  20

% In this code, the seismograms are not deconvolved. Instead the Green fcns
% are convolved with the seismometer responses. This enables inclusion of
% tilt data

clear
close all

% if you want to get data from the database, use loadfresh=1
loadfresh=0;
plot_flag=0; % use plot_flag=2 for more plots, incl. the GFs. This is for testing only
usetilt = 1;
outputfiles_flag=0;
if outputfiles_flag
    disp('writing out all the STFs will slow this down')
end

%high corner of filter
HCall=[1/60, .1 .1 .1 .1 .1];
% decide with low corner filter to use. These 4 are a good sample
LCall=1./[400 400 120 60 30 90];
tp=0.5; disp(['using tp = ',num2str(tp)]);
% HCall=.1;
% LCall=1/400;
% This works well for 2009 data
chnwgtstr='22chaPWS';
year_flag=2009;
%2009
stasub=[1:3,6:8];

%choose event 1, 2, or 3
event=4;
if event==1
    % this one is really good: Event 2009.1
    datastart=datenum(2009,1,19,16,9,25)-10/(24*60)-40.96/86400;
    for n=1:length(LCall)
        tmp=sprintf('event1_%03d-%03d.mat',1/LCall(n),1/HCall(n));
        eventfile(n,:)={tmp};
    end
elseif event==2
    % this is a little noisier: Event 2009.2
    datastart=datenum(2009,1,19,18,45,5)-10/(24*60)-40.96/86400;
elseif event==3
    % this is also very good: Event 2009.3
    datastart=datenum(2009,1,19,21,10,10)-10/(24*60)-40.96/86400;
    for n=1:length(LCall)
        tmp=sprintf('event3_%03d.mat',1/LCall(n));
        eventfile(n,:)={tmp};
    end
elseif event==4
    for n=1:length(LCall)
        datastart=datenum(2009,1,1); % dummy value
        tmp=sprintf('PWstack2009BP%d_%d.mat',1/LCall(n),1/HCall(n));
        eventfile(n,:)={tmp};
    end
    stasub=[1:3,6:9];
    chnwgtstr='22chaPWS';
    %Need to switch
    
end

dataend = datastart+20/(24*60);


ratmult=1;
% factor to scale the tilt gfs to something closer to the trans gfs

nummom=6; % this does not work yet for forces only

% This is the directory where the GFs are located. They are in individual
% directories labeled with node x-y-z (like 117-099-122) and labeled like
% greenf.1.ft and greenf.1.rot.ft for the translation and rotation GFs
% respectivley. Actually, the rotation GFs are simply the curl of the
% displacement field so the rotation is half of this.
rootdir='/Users/gpwaite/Data/Guate/Fuego2012/FuegoVLP2012/fdmtilt';

%make sure the data are consistent with padded GFs
td=2^15;
fd=td/2;

if loadfresh
    % use this for td=16383
    smallwinstart=16162;
    smallwinend = 48929;
    
end
% this simply cuts the ends off the trace before computing errors
timetrim=[10,645];

% 2012
%stasub=[11:15,17:20];

% determine weighting here

% 27 weight (everything except N1)
weight=ones(1,27);


% 21cha
if strcmp(chnwgtstr,'17cha')
    weight=ones(1,18);
    weightrot=ones(1,18);
    weight(13)=1e-20;
    weight(14)=1e-20;
    
    weightrot(13)=1e-20;
    weightrot(14)=1e-20;
    
    % all the z compoents of tilt shoudl be zero
    %     weightrot(3)=1e-20;
    %     weightrot(6)=1e-20;
    %     weightrot(9)=1e-20;
    %     weightrot(12)=1e-20;
    %     weightrot(15)=1e-20;
    %     weightrot(18)=1e-20;
    
    % 23cha
elseif strcmp(chnwgtstr,'23cha')
    weight(3)=1e-20;
    weight(10)=1e-20;
    weight(11)=1e-20;
    weight(22)=1e-20;
    % 25 cha
elseif strcmp(chnwgtstr,'25cha')
    weight(10)=1e-20;
    weight(22)=1e-20;
elseif strcmp(chnwgtstr,'27cha')
elseif strcmp(chnwgtstr,'22chaPWS')
    weight=ones(1,21);
    
    % don't use F9B E, which is a dead channel
    weight(13)=1e-20;
    
    % what if we throw out NE Z, which is crazy
    weight(21)=1e-20;
    
    weightrot=weight;
else
    error([chnwgtstr,' does not match anything'])
end
nonzerocha=floor(sum(weight));

%% load the stf used to make the GFs or just create it on the fly
% this will be used to convolve with the data prior to inversion
%
% assume 50 sps
tvecstf=(0:td-1)*.02;
stftd=zeros(td,1);
aa=find(tvecstf<tp);
for n=1:length(aa)
    stftd(aa(n)) = (1-cos(2*pi*tvecstf(aa(n))/tp))/2;
end

tmp=fft(stftd);
stf=tmp(1:fd,1);
fvecstf=linspace(0,25,fd);
if plot_flag==2
    figure; subplot(211); plot(tvecstf,stftd);
    subplot(212); plot(fvecstf,abs(stf));
end

%% start the loop over different lower corners
% since we filter the data, it is easiest to do this outside the part where
% the data is loaded.
tic
for locor=1:length(LCall)
    LC=LCall(locor);
    HC=HCall(locor);
    
    % load the data
    
    % Here we will get the time domain seismograms from the mat files using
    % GISMO. This should allow for easier testing of many different events.
    % staall=[{'SW1'};{'SW2'};{'S1'};{'SE1'};{'SE2'};{'N1'};{'NE1'};{'NE2'};{'NW1'};{'NW2'};];
    
    %These are espc: F9SW, F9NW, F904, F9B, F9NE,
    staall=[{'F9SW'};{'F9NW'};{'F900'};{'F9E'};{'F9D'};{'F9C'};{'F9B'};{'F9A'};{'F9NE'};{'F9SE'};{'SW1'};{'SW2'};{'S1'};{'SE1'};{'SE2'};{'N'};{'NE1'};{'NE2'};{'NW1'};{'NW2'}];
    sta=staall(stasub);
    % this is the order that the data are in in the GFs
    allchn=[{'HHE'};{'HHN'};{'HHZ'}];
    inst=[{'F9_espc'};{'F9_espc'};{'F9_cmg40t'};{'F9_cmg40t'};{'F9_cmg40t'};{'F9_cmg40t'};{'F9_espc'};{'F9_cmg40t'};{'F9_espc'};{'F9_cmg40t'};{'F12_cmg40t'};{'F12_cmg40t'};{'F12_cmg40t'};{'F12_cmg40t'};{'F12_cmg40t'};{'trilliumcompact'};{'trilliumcompact'};{'F12_cmg40t'};{'trilliumcompact'};{'F12_cmg40t'};];
    instsub=inst(stasub);
    if year_flag==2009
        net='XT';
    elseif year_flag==2012
        net='XJ';
    end
    loc='01';
    tic
    nt=0;
    if loadfresh
        for n=1:length(sta)
            stn=char(sta(n));
            % remember that this station had a bad east channel, so we rotated it
            % 90 degrees
            if strcmp('F9B',stn)
                allchn=[{'HHN'};{'HHE'};{'HHZ'}];
                disp('swapping horizontal channels for F9B')
                
            else
                allchn=[{'HHE'};{'HHN'};{'HHZ'}];
                
            end
            for m=1:length(allchn)
                nt=nt+1;
                chn=char(allchn(m));
                if year_flag==2009
                    loadsource=datasource('file','/Users/gpwaite/Data/Guate/Fuego2009/Fuego_2009_mat/%s/%s%s%04d%02d%02d.mat','station','station','channel','year','month','day');
                elseif year_flag==2012
                    loadsource=datasource('file','/Users/gpwaite/Data/Guate/Fuego2012/Fuego_2012_mat/%s/%s%s%04d%02d%02d.mat','station','station','channel','year','month','day');
                end
                scnl = scnlobject(stn,chn,net,loc);
                wo=waveform(loadsource,scnl,datastart,dataend);
                wo=set(wo,'data',double(wo)-mean(double(wo)));
                if isempty(wo)
                    error(['no data for ',scnl])
                end
                
                % filter
                spsorig=get(wo,'freq');
                fo=filterobject('B',[LC HC],2);
                wof=filtfilt(fo,wo);
                %             disp(['filtering from ',num2str(1/LC),'-',num2str(1/HC),'s period'])
                
                %         extract smaller window and decimate
                decifac=2;
                sps=spsorig/decifac;
                tmpraw=double(wof);
                tmpraw50=decimate(tmpraw,decifac);
                tmpraw=tmpraw50(smallwinstart:smallwinend);
                
                dt=1/sps;
                
                % apply cosine taper
                % .5 is fine as long as the window is well centered
                taperwin=tukeywin(td,.5); %
                
                % have not been deconvolved
                tmpraw=tmpraw.*taperwin;
                
                % this will be used later
                
                combdatatd(nt,:)=tmpraw;
                if plot_flag==2
                    figure(1)
                    subplot(211)
                    plot(combdatatd(nt,:));
                    subplot(212)
                    plot(tmpraw);
                    suptitle2([get(scnl,'station'),':',get(scnl,'channel')])
                    drawnow
                    pause;
                end
                
                if (length(tmpraw)~=td), error('data not trimmed correctly'); end
                
                % data should be in counts, make sure the range is large
                if (range(double(wo))<1e5), disp('check data to be sure they are in counts'); end
                
                datatd(nt,:)=tmpraw;
                % FFT
                tmpfft = fft(tmpraw);
                
                % put in matrix
                data(nt,:)=tmpfft(1:fd);
                clear wo wof tmpfft tmpVLPraw
                
                %             disp(['finished with ',stn,':',chn,' in ',num2str(toc)])
            end
        end
        save(char(eventfile(locor)),'data','sps','combdatatd','datatd')
    else
        load(char(eventfile(locor)))
        sps=50;
        combdatatd=datatd;
    end
    
    
    
    % stations are in this order for the regular for 2012
    % 3360        3560	SW1 40T 800 V/m/s
    % 3040        3280	SW2 40T 800 V/m/s
    % 3640        1560	S1  40T 800 V/m/s
    % 6440        3680	SE1 40T 800 V/m/s
    % 6840        3400	SE2 40T 800 V/m/s
    % 5000        5280	N1  Trillium Compact
    % 5040        5040	NE1 Trillium Compact (S/N 181) 749.1 V/m/s
    % 5480        5280	NE2 40T 800 V/m/s (Z-724, NS-992, EW-814)
    % 4920        5040	NW1 Trillium Compact
    % 4360        5400	NW2 40T 800 V/m/s
    %
    % generate a matrix of translation transfer functions
    % use all the stations so the number matches the number of green fcns
    chncnt=0;
    % freqvec=make_freq(td,1/sps);
    freqvec=linspace(0,sps,td);
    for nsta=1:length(inst)
        calib=get_inst_resp(char(inst(nsta)));
        %     [tr,re,f,t]=response(calib,1,freqvec,[]);
        % digitizer and sensitivity are applied outside the response
        %     tr=response(calib,1,freqvec,[])*calib(2);
        % what if I just do this independently of the goofy ff used in the
        % response function
        npole =real(calib( 3,1));
        nzero =real(calib(33,1));
        poles =calib( 4:npole+ 3,1);
        zeroes=calib(34:nzero+33,1);
        bb    =poly(zeroes);             % convert zeros to polynomial coefficients
        aa    =poly(poles);              % convert poles to polynomial coefficients
        
        % check that normalization is correct
        testnorm=abs(polyval(poly(poles),2*pi*1i)/polyval(poly(zeroes),2*pi*1i));
        
        norm=calib(1)*calib(2);
        % everything should be in rad/s for freqs (Laplace transform)
        w=2*pi*freqvec;
        tr=freqs(bb,aa,w)*norm;
        
        %     here we use 3 channels
        for dum=1:3
            chncnt=chncnt+1;
            transtransmat(chncnt,:)=tr(1:fd);
        end
    end
    
    
    %
    if tp<1
        
        allGFdirs=dir([rootdir,'/*-*-*']);
    else
        allGFdirs=dir([rootdir,'/*-*-*-tp',num2str(tp)]);
        
    end
    gind=[];
    cnt=1;
    for nd=1:length(allGFdirs)
        
        if exist([rootdir,'/',allGFdirs(nd).name,'/greenf.7.ft'],'file')
            gind(cnt)=nd;
            cnt=cnt+1;
        end
    end
    allGFs=allGFdirs(gind);
    if size(allGFs,1)< size(allGFdirs,1), disp('missing some greenf files');end
    
    
    % load the GFs, first the translation
    
    % chasub is a vector of indices for the GFs that actually correspond to the
    % stations and channels
    chasub=sort([3*(stasub-1)+1,3*(stasub-1)+2,3*(stasub-1)+3]);
    
    if nummom==9
        allmom=1:9;
    elseif nummom==6
        allmom=1:6;
    elseif nummom==3
        allmom=7:9;
    end
    E1=zeros(1,length(allGFs));
    E2=E1;
    
    %%% if you only want to compute a single solution with event 1
    %         best_504_LC4=[252,264,271,215,208];
    %         for ng=best_504_LC4(locor)
    
    % this is E2 best fit for 120-600: 123-105-140 of 0.3683
    %         for ng=252
    % this is E1 best fit for 120-600: 117-107-140 of 0.021455
    % for ng=47
    tic
    %
    %         best_PWstack=[404,211,137,57,188]
    % after selection from best gamma within 10% of the min E2
    %                 best_PWstack=[469 469 293 54 54];
%     best_PWstack=[469 469 293 54 137];
% %     
%  disp(['results are from 22chaPWS_6moment_output_expandedZ__flippedrotx2016_04_08_16_20_49_2009PWS.mat'])
%  disp('and refined with Matozas gamma metric');
%  best_PWstack=[459 190 82 45 103 82]; % this is from Matoza's gamma metric
% % 
% % % 
%     for ng=best_PWstack(locor)
        
        % %509 is 129-111-142
        % testhigh=[509]
        %             for ng=testhigh(1)
        %
        %%% otherwise, do all
    disp(['inverting for ',num2str(size(allGFs,1)),' nodes'])
    for ng=1:length(allGFs)
        
        for nm=allmom
            % In the forward models (with topo_v17) I used a larger scalar
            % seismic moment for the off-diagonal components (1 vs 1/sqrt(2)
            % for the diagonal) because I wanted to keep them symmetric. This
            % corrects that, although it seems to make very little difference
            if nm==4 || nm==5 || nm==6
                scalarmomentnormfac=1/sqrt(2);
                
            else
                scalarmomentnormfac=1;
            end
            % use smaller tp for trans
            %         file=sprintf('/Users/gwaite/FuegoVLP2012/fdm/%s/greenf.%d.ft',allGFs(ng).name,nm);
            gfile=sprintf('%s/%s/greenf.%d.ft',rootdir,allGFs(ng).name,nm);
            %         disp(['loading ',gfile]);
            fid=fopen(gfile,'r');
            datahead=textscan(fid,'%d %d %f %d',1);
            a=1:9:88;
            b=a+1;
            c=b+1;
            actual=[a,b,c];
            for nc=1:datahead{2}
                datacell=textscan(fid,'%f %f',datahead{1});
                
                datar=datacell{1};
                datai=datacell{2};
                if plot_flag==2
                    figure(99); subplot(311); plot(abs(complex(datar,datai))); hold on
                end
                % pad with zeros here
                tmp=zeros(datahead{4},1);
                tmp(1:datahead{1},1)=complex(datar,datai);
                tmp1g=ifft(tmp,'symmetric');
                tmp1=tmp1g*scalarmomentnormfac;
                if plot_flag==2
                    figure(99); subplot(312); plot(tmp1)
                end
                % check to see if it is reversed
                if rms(tmp1(1:datahead{1}))>rms(tmp1(datahead{1}+1:datahead{4})), error('not reversed');end
                tmptd=zeros(td,1);
                tmptd(1:length(tmp1),1)=flipud(tmp1);
                Gtd(nc,nm,:)=flipud(tmp1);
                if plot_flag==2
                    subplot(313); plot(tmptd)
                    suptitle2(num2str(nc))
                    if find(nc==actual)
                        pause(1)
                    end
                    close(99)
                end
                tmpfd=fft(tmptd);
                G(nc,nm,:)=tmpfd(1:fd);
                clear datacell
            end
            fclose(fid);
            clear datahead tmp
        end
        
        % next the load the rotations
        % in the first column, the rotation is w.r.t. the x axis, so should be
        % related to tilt recorded on the x component
        %in the second column, the rotation is w.r.t. the y axis, so should be
        %relatedto the y (N) tilt
        % in the last column we have rotation w.r.t. z, which we do not measure
        for nm=allmom
            % In the forward models (with topo_v17) I used a larger scalar
            % seismic moment for the off-diagonal components (1 vs 1/sqrt(2)
            % for the diagonal) because I wanted to keep them symmetric. This
            % corrects that, although it seems to make very little difference
            if nm==4 || nm==5 || nm==6
                scalarmomentnormfac=1/sqrt(2);
            else
                scalarmomentnormfac=1;
            end
            
            gfile=sprintf('%s/%s/greenf.%d.rot.ft',rootdir,allGFs(ng).name,nm);
            fid=fopen(gfile,'r');
            
            datahead=textscan(fid,'%d %d %f %d',1);
            for nc=1:datahead{2}
                datacell=textscan(fid,'%f %f',datahead{1});
                
                datar=datacell{1};
                datai=datacell{2};
                if plot_flag==2
                    figure(99);
                    fvec=linspace(0,0.5/datahead{3},datahead{1});
                    subplot(311); plot(fvec,abs(complex(datar,datai))); hold on
                end
                % pad with zeros here
                tmp=zeros(datahead{4},1);
                tmp(1:datahead{1},1)=complex(datar,datai);
                tmp1=ifft(tmp,'symmetric');
                
                if plot_flag==2
                    figure(99); subplot(312); tvec=(0:double(datahead{4}-1))*datahead{3}; plot(tvec,tmp1)
                end
                % check to see if it is reversed
                if rms(tmp1(1:datahead{1}))>rms(tmp1(datahead{1}+1:datahead{4})), disp('not reversed');end
                
                % these synthetics are the curl of the displacement, we need to
                % divide by 2 to get to rotation in radians
                tmp1=tmp1/2;
                %             tmp1=1e6*tmp1/2;
                %             % multiply the tilt GFs by a scale factor
                %             tmp1=tmp1*ratmult*scalarmomentnormfac;
                
                tmptd=zeros(td,1);
                tmptd(1:length(tmp1),1)=flipud(tmp1);
                
                % filter these
                [B,A]=butter(2,(1/6)/25,'low');
                tmptd=filtfilt(B,A,tmptd);
                
                Grottd(nc,nm,:)=tmptd(1:length(tmp1));
                if plot_flag==2
                    subplot(313); plot(tmptd)
                    pause(1)
                end
                tmpfd=fft(tmptd);
                Grotuncorr(nc,nm,:)=tmpfd(1:fd);
                clear datacell
            end
            fclose(fid);
            clear datahead tmp
        end
        % here is where we sort out the tilt GF so that it makes sense
        % the sign of the tilt is opposite to the way it is described in
        % Lyons et al 2012, in which positive N tilt is down to the N.
        % Here, positive N tilt would come from negative rotx, but the sign
        % convention for roty (E tilt) is consistent with Lyons et al 2012
        if usetilt
            Grot=Grotuncorr;
            for nc=1:3:size(Grotuncorr,1)
                Grot(nc,:,:) =  Grotuncorr(nc+1,:,:); % roty is E-W tilt so it goes first
                Grot(nc+1,:,:)=-Grotuncorr(nc,:,:); % rotx is N-S tilt and has the opposite sign!
                Grot(nc+2,:,:)= zeros(1,nummom,fd); % rotz should be zero
            end
        end
        
        %     % compute tilt GFs from displacement GFs
        %     rotx = ( Uz(l,m,n) - Uz(l,m-1,n  ) ) - ( Uy(l,m,n) - Uy(l,m,  n-1) )
        %     rotx = rotx / dh
        %     roty = ( Ux(l,m,n) - Ux(l,  m,n-1) ) - ( Uz(l,m,n) - Uz(l-1,m,n  ) )
        %     roty = roty / dh
        %     rotz = ( Uy(l,m,n) - Uy(l-1,m,  n) ) - ( Ux(l,m,n) - Ux(l,  m-1,n) )
        %     rotz = rotz / dh
        %
        %     % in G stations are in 1, 4, 7, etc
        %     % in 2, 5, 8, etc are surface positions 40 m E of the station
        %     % in 3, 6, 9, etc. are surface positions 40 m N of the station
        %     stasub=1:3:size(Gtd,1);
        %     Gtdsta=Gtd(stasub,:,:);
        %     stasubE=2:3:size(Gtd,1);
        %     GtdstaE=Gtd(stasub,:,:);
        %
        %     stasubN=3:3:size(Gtd,1);
        %     GtdstaN=Gtd(stasub,:,:);
        %
        %     for nm=allmom
        %         stacnt=1;
        %         for nc=1:3:size(Gtd,1)
        %             ewtilttd(stacnt,nm,:)=(Gtd(nc,nm,:)-Gtd(nc+1,nm,:))/40;
        %             nstilttd(stacnt,nm,:)=(Gtd(nc,nm,:)-Gtd(nc+2,nm,:))/40;
        %             stacnt=stacnt+1;
        %         end
        %     end
        %
        %normalize the moments for each
        % this should be done before convolving with the instrument responses
        normwgt=zeros(nummom,1);
        %     for nm=allmom
        %         tmp(:,:)=G(:,nm,:);
        %         normwgt(nm)=max(max(abs(tmp)));
        %         tmp=tmp/normwgt(nm);
        %         G(:,nm,:)=tmp;
        %         clear tmp
        %     end
        %
        %%%%%%%%%%%%%
        % Convolve the GFs with the appropriate transfer functions
        % there are 90 channels total in the G matrix, so use them all
        for ncha=1:size(transtransmat,1)
            tmptrans=transtransmat(ncha,:);
            for nm=1:nummom
                tmpG(1,:)=G(ncha,nm,:);
                %             size(tmpG)
                %             size(tmptrans)
                %stabilize as in deconvolution? with water level
                %             for nf=1:2000
                %                 G(ncha,nm,:)=tmpG(nf).*tmptrans(nf);
                %                 G(ncha,nm,:)=tmpG(nf).*(tmptrans(nf).*conj(tmptrans(nf))+1e-8)./conj(tmptrans(nf));
                %             end
                G(ncha,nm,:)=tmpG.*tmptrans;
            end
            clear tmptrans
            
        end
        clear tmpG
        if usetilt
            % create the tilt transfer function
            % here we use the translation transfer function transtransmat
            % multiplied by a vector representing the tilt transfer
            % constant
%             disp('need to check tilt transfer function')
            fvecfd=linspace(0,sps,fd);
            constant=-9.8./(2*pi*fvecfd).^2; %% this is from Maeda et al., GJI, 2011
            constant(1)=constant(2);%-1e-9; %First frequency is zero so change it to something small
            for ncha=1:size(transtransmat,1)
                tmptrans=transtransmat(ncha,:);
                % we want zero response on the verticals
                if ~mod(nc,3), tmptrans=0; end
                for nm=1:nummom
                    tmpGrot(1,:)=Grot(ncha,nm,:);
                    Grot(ncha,nm,:)=tmpGrot.*tmptrans.*constant;
                end
                clear tmptrans
                
            end
            clear tmpGrot
        end
        
        
        % normalize the GFs and save the norm weights for unnorming the STF
        % later
        for nm=allmom
            tmpG(:,:)=G(:,nm,:);
            if usetilt
                tmpGrot(:,:)=Grot(:,nm,:);
                % add these together
                tmpGcomb=tmpG+tmpGrot;
                normwgt(nm)=max(max(abs(tmpGcomb)));
                tmpGcomb=tmpGcomb/normwgt(nm);
                G(:,nm,:)=tmpGcomb;
                clear tmpGcomb tmpG tmpGrot
            else
                normwgt(nm)=max(max(abs(tmpG)));
                G(:,nm,:)=tmpG/normwgt(nm);
                clear tmpG
            end
        end
        
        
        
        %%%%%
        
        % save a copy of the GFs for generating synthetics
        % this should have normalization applied, but not the
        % instrument responses or weighting or else we won't get synthetics for
        % zero-weighted data
        Gorig=G;
        if usetilt
            Grotorig=Grot;
        end
        
        %weight GFs
        for nt=1:length(chasub)
            G(chasub(nt),:,:)=G(chasub(nt),:,:)*weight(nt);
        end
        if nt~=length(weight), error('check weighting vector'); end
        
        % convolve with instrument responses
        % since the data are not deconvolved, we convolve the GFs with the
        % instrument responses
        % However, we will not do that with the tilt
        
        % check the responses?
        if 0
            f=linspace(0,25,fd);
            for ncha=1:size(transtransmat,1)
                transfer=transtransmat(ncha,:);
                mag=abs(transfer); phase=angle(transfer);
                figure
                hold off
                clf
                subplot(2,1,1)
                loglog(f,mag);xlabel('frequency (Hz)');ylabel('amplitude');
                grid;title('transfer function')
                subplot(2,1,2)
                semilogx(f,unwrap(phase));xlabel('frequency (Hz)');ylabel('phase');
                grid;title('transfer function')
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % invert and generate FD synthetics
        % it isn't necessary to invert over the whole range of frequencies
        % since the sample freq is 50 and the largest signals are well under 1
        % Hz. It would be quicker to restrict this to just the first 100
        % frequencies
        
        m=zeros(nummom,fd);
        % skip the first frequency because it never works
        %for nf=2:fd
        % skip the frequencies above 1 Hz
        for nf=2:657
            % G is either translation only or translation plus rotation GFs
            tmpG=zeros(length(chasub),nummom);
            tmpG(:,:)=G(chasub,:,nf);
            tmpdata=data(:,nf);
            
            % this will do LS inversion
            if sum(sum(isnan(tmpG)))>0
                disp(['nans in G for ',allGFs(ng).name])
                m(:,nf)=zeros(nummom,1);
            elseif rank(tmpG)<nummom
                disp(['rank is too low for freq number ',num2str(nf)])
                m(:,nf)=zeros(nummom,1);
            else %invert
                % multiply data with stf since the GFs are convolved with it
                m(:,nf)=tmpG\(tmpdata*stf(nf)); %#ok<*SAGROW>
                % using Matlab's \ function is equivalent to:
                % m(:,nf)=(tmpG'*tmpG)\(tmpG'*tmpdata*stf(nf));
            end
            clear tmpG tmpdata
        end
        % un-norm the moments
        for nm=1:nummom
            m(nm,:)=m(nm,:)/normwgt(nm);
        end
        
        synthfd=zeros(size(data,1),fd);
        if usetilt
            synthfd=zeros(size(data,1),fd);
        end
        for nf=1:fd
            %         tmpG(:,:)=Gorig(chasub,:,nf);
            tmpG(:,:)=G(chasub,:,nf);
            % we need to divide by the stf since the GF is convolved with it
            % moment components are multiplied by the weighting vector used to
            % invert
            %         if ~usetilt
            synthfd(:,nf)=tmpG*(m(:,nf).*normwgt)/stf(nf);
            %         else
            %             %             tmpGrot(:,:)=Grotorig(chasub,:,nf);
            %             tmpGrot(:,:)=Grot(chasub,:,nf);
            %             totalG=[tmpG+tmpGrot];
            %             abc=max(max(abs(tmpG)));
            %             def=max(max(abs(tmpGrot)));
            %             ratio(nf)=abc/def;
            %             synthfd(:,nf)=totalG*(m(:,nf).*normwgt)/stf(nf);
            %         end
            clear tmpG
            
        end
        % generate synthetics
        for nc=1:size(synthfd,1)
            tmp=zeros(1,td);
            tmp(1:fd)=synthfd(nc,:);
            synthtd(nc,:)=ifft(tmp,'symmetric');
        end
        
        timetrimsamp=timetrim*sps;
        tts=timetrimsamp;
        sysub=synthtd(:,tts(1):tts(2));
        if usetilt
            %  dasub=[datatd(:,tts(1):tts(2))+datatdt(:,tts(1):tts(2))];
            dasub= combdatatd(:,tts(1):tts(2));
        else
            dasub= datatd(:,tts(1):tts(2));
        end
        
        for ns=1:length(sta)
            for nc=1:length(allchn)
                nt=(ns-1)*3+nc;
                num1(nt)=sum( ( dasub(nt,:) - sysub(nt,:) ).^2 * weight(nt));
                num2(nc)=sum( ( dasub(nt,:) - sysub(nt,:) ).^2 * weight(nt));
                den1(nt)=sum( dasub(nt,:).^2 );
                den2(nc)=sum( dasub(nt,:).^2 );
            end
            E2tmp(ns)=sum(num2)/sum(den2);
        end
        E1(ng)=sum(num1)/sum(den1);
        
        E2(ng)=sum(E2tmp)/length(sta);
        E1all(locor,ng)=E1(ng);
        E2all(locor,ng)=E2(ng);
        
        if ~mod((100*ng/length(E1)),10)
            disp([num2str(100*ng/length(E1)),'% done with LC ',num2str(1/LC),' HC ',num2str(1/HC)])
            ttt=toc;
            disp([num2str(ttt/60),' minutes elapsed'])
            disp([num2str(((ttt*length(E1)/ng)-ttt)/60),' minutes left'])
        end
        
        if outputfiles_flag
            stfoutfile=sprintf('%s/%s/STF_%03d-%03d.dat',rootdir,allGFs(ng).name,round(1/LC),round(1/HC));
            fid=fopen(stfoutfile,'w');
            fprintf(fid,'%d %d %f %d\n',td,size(m,1),1/sps,fd);
            for nm=1:size(m,1)
                tmpm=zeros(1,td);
                tmpm(1:fd)=m(nm,:);
                mtd(nm,:)=ifft(tmpm,'symmetric');
                fprintf(fid,'%f\n',mtd(nm,:));
            end
            fclose(fid);
            %         stfdir=sprintf('%s/%s',rootdir,allGFs(ng).name);
        end
        
    end
    % if usetilt && 0
    %     ratiotilttransdata=max(max(abs(datarot')))/max(max(abs(data')))
    %     ratiotilttransGF=max(max(max(abs(Grot))))/max(max(max(abs(G))))
    %     ratiotilttransGF=max(max(max(abs(Grotorig))))/max(max(max(abs(Gorig))))
    % end
    %
    [~,bbb]=find(E1>0);
    [a,b]=min(E1(bbb));
    disp(['min E1 at ',allGFs(bbb(b)).name,' of ',num2str(a)]);
    disp(['x=',num2str(40*str2double(allGFs(b).name(1:3))),' y=',num2str(40*str2double(allGFs(b).name(5:7)))])
    
    [aaa,bbb]=find(E2>0);
    [a,b]=min(E2(bbb));
    disp(['min E2 at ',allGFs(bbb(b)).name,' of ',num2str(a)]);
    disp(['x=',num2str(40*str2double(allGFs(b).name(1:3))),' y=',num2str(40*str2double(allGFs(b).name(5:7)))])
    
    disp('the best fit VLP source from Lyons and Waite, 2011 is 119-107-134')
    
    % save some output
    dte = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    savefile=sprintf('%s_%dmoment_output_expandedZ_flippedrotx%s.mat',chnwgtstr,nummom,dte);
    
    % Compute PDF
    
    %     PDFall(locor,:) = get_PDF(E2,mean(DoF,'omitnan'),6);
    
    
    toc
end
if event==1
    savefile=sprintf('%s_%dmoment_output_expandedZ_%s_Event2009.1.mat',chnwgtstr,nummom,dte);
elseif event==2
    savefile=sprintf('%s_%dmoment_output_expandedZ_%s_Event2009.2.mat',chnwgtstr,nummom,dte);
elseif event==3
    savefile=sprintf('%s_%dmoment_output_expandedZ_%s_Event2009.3.mat',chnwgtstr,nummom,dte);
elseif event==4
    savefile=sprintf('%s_%dmoment_output_expandedZ__flippedrotx%s_2009PWS.mat',chnwgtstr,nummom,dte);
end
if 1
    save(savefile,'allGFs','E1all','E2all','LCall','HCall');
end
figure;
subplot(211)
plot(E1all'); title('E1')
legend(num2str(1./LCall'))

hold on
[a,b]=min(E1all,[],2);
plot(b,a,'*')
subplot(212)
plot(E2all'); title('E2')
legend(num2str(1./LCall'))

hold on
[a,b]=min(E2all,[],2);
plot(b,a,'*')


stfdir=sprintf('%s/%s',rootdir,allGFs(ng).name);

%     if 0
%         %     file=sprintf('%s/%s/STF.dat',rootdir,allGFs(ng).name);
%         stfoutfile=sprintf('%s/%s/STF_%03d.dat',rootdir,allGFs(ng).name,round(1/LC));
%
%         fid=fopen(stfoutfile,'w');
%         fprintf(fid,'%d %d %f %d',td,size(m,1),dt,fd);
%         for nm=1:size(m,1)
%             tmpm=zeros(1,td);
%             tmpm(1:fd)=m(nm,:);
%             mtd(nm,:)=ifft(tmpm,'symmetric');
%             fprintf(fid,'%f\n',mtd(nm,:));
%         end
%         fclose(fid);
%     end
% %%
% figure
% tvec=(0:td-1)*dt;
% j=0;
% for ns=1:length(sta)
%     for nc=1:length(cha)
%         nt=(ns-1)*3+nc;
%         j=j+1;
%         ht(j)=subplot(ns+1,3,nt);
%         plot(tvec,datatd(nt,:)); hold on
%         plot(tvec,synthtd(nt,:),'r--')
%     end
% end

%% plot residuals
if plot_flag
    tit=[datestr(datastart+10/(24*60)+40.96/86400),' ',allGFs(ng).name];
    %     for nt=1:size(combdata,1)
    %     tmp1=zeros(1,td);
    %     tmp1(1:fd)=combdata(nt,:);
    %     dataraw(nt,:)=ifft(tmp1,'symmetric');
    %     end
    if usetilt
        
        plotresiduals(combdatatd,synthtd(1:21,:),weight,sta,tit,270,380,1/sps,0);
        print('-dpsc',['res_vel',num2str(sum(weight)),'cha_',num2str(round(1/LC)),'-',num2str(round(1/HC)),'s_',allGFs(ng).name,'.ps'])
        
        %integrate the data
        for nt=1:size(combdatatd,1)
            intdatatd(nt,:) =cumtrapz(combdatatd(nt,:))/50;
            intsynthtd(nt,:)=cumtrapz(synthtd(nt,:))/50;
        end
        plotresiduals(intdatatd,intsynthtd,weight,sta,tit,1,655,1/sps,0);
        %         print('-dpsc','residualsint.ps')
        print('-dpsc',['res_disp',num2str(sum(weight)),'cha_',num2str(round(1/LC)),'-',num2str(round(1/HC)),'s_',allGFs(ng).name,'.ps']);
        
        tmp=pwd;
        cd (stfdir)
        
        stfoutfile=sprintf('%s/%s/STF_%03d-%03d.dat',rootdir,allGFs(ng).name,round(1/LC),round(1/HC));
        stfmfilename(stfoutfile,1,1,655,tit)
        cd(tmp)
    else
        plotresiduals(datatd,synthtd(1:18,:),weight,sta,tit,283,373,1/sps,0);
        print('-dpsc','residuals.ps')
        tmp=pwd;
        cd (stfdir)
        stfm(1,283,373,tit)
        cd(tmp)
    end
end


%% plot error volumes


clear tmp
if 0
    for ng=1:length(E2)
        x(ng)=str2double(allGFs(ng).name(1:3));
        y(ng)=str2double(allGFs(ng).name(5:7));
        z(ng)=str2double(allGFs(ng).name(9:11));
    end
    xmin=min(x);
    ymin=min(y);
    zmin=min(z);
    xmax=max(x);
    ymax=max(y);
    zmax=max(z);
    
    for n=1:length(E2)
        xi=x(n)-xmin+1;
        yi=y(n)-ymin+1;
        zi=z(n)-zmin+1;
        v(xi,yi,zi)=E2(n);
    end
    
    %     xs=xmin:xmax;
    %     ys=ymin:ymax;
    %     zs=zmin:zmax;
    % here we skipped every other node
    v=v(1:2:end,1:2:end,1:2:end);
    
    xs=xmin:2:xmax;
    ys=ymin:2:ymax;
    zs=zmin:2:zmax;
    [xnew,ynew,znew]=meshgrid(ys,xs,zs);
    vi=interp3(ys,xs,zs,v,xnew,ynew,znew);
    
    for n=1:size(vi,3)
        
        tmp(:,:)=vi(:,:,n);
        
        for i=1:size(tmp,1)
            for j=1:size(tmp,2)
                if tmp(i,j)==0, tmp(i,j)=NaN;
                end
            end
        end
        map1=jet;
        map1(1,:)=[0 0 0];
        figure
        imagesc(xs,ys,tmp',[min(E2) max(E2)])
        
        %transfer the original so xs is columns, ys is rows
        vinew(:,:,n)=tmp';
        
        axis equal
        axis tight
        axis xy
        colormap(map1)
        colorbar('SouthOutside')
        title(['layer at z = ',num2str(zs(n)),'  layer min = ',num2str(min(min(tmp)))])
        
    end
    
    
end



% For 21cha6M, between 130 and 135, the error space is well-defined,
% roughly in themiddle of the search volume.

%% try smoothing and plotting a subset
%
% [nx,ny,nz,vsub]=subvolume(vi,[3,12,3,19,21,27]);
% nx=nx+ys(1)-1;
% ny=ny+xs(1)-1;
% nz=nz+zs(1)-1;
if 0
    
    %         [nx,ny,nz,vsub]=subvolume(vinew,[3,19,3,12,21,27]);
    
    % using this for the subsets?
    %     vtmp=vi(1:2:end,1:2:end,1:2:end);
    [nx,ny,nz,vsub]=subvolume(vinew,[1,7,1,7,1,12]);
    
    nx=nx*2+xs(1)-2;
    ny=ny*2+ys(1)-2;
    nz=nz*2+zs(1)-2;
    
    %     vs=smooth3(vsub);
    vs=vsub;
    
    % xs(3:19),ys(3:12),zs(20:27),
    isolevel=0.3683*1.25;
    figure
    patch(isocaps(nx,ny,nz,vs,isolevel,'below'),'FaceColor','interp','EdgeColor','none');
    p1 = patch(isosurface(nx,ny,nz,vs,isolevel),'FaceColor','blue','EdgeColor','none');
    isonormals(vs,p1)
    axis equal
    view(3);
    axis vis3d
    grid on
    camlight left;
    lighting phong
    
end

%%
if usetilt && 0
    for nt=1:3:18
        figure
        for nc=0:2
            tmp=zeros(1,td);
            tmp(1:fd)=combdata(nt+nc,:);
            tmptd=ifft(tmp,'symmetric');
            tvec=(0:td-1)/50;
            figure(10)
            subplot(3,1,(1+nc))
            
            plot(tvec,tmptd,tvec,combdatatd(nt+nc,:))
            figure(11)
            subplot(3,1,(1+nc))
            
            plot(tvec,cumtrapz(tmptd)/50,tvec,cumtrapz(combdatatd(nt+nc,:))/50)
        end
        figure(10)
        suptitle2(sta{(nt-1)/3+1})
        figure(11)
        suptitle2(sta{(nt-1)/3+1})
        drawnow
        pause
    end
end

