% Invert VLP data from either 2009 or 2012 using fixed moment tensor. This
% is designed for either a full grid search, or for use in Monte-Carlo
% style inversion. It is only designed to evaluate sources at a single
% location. Use the code invVLPtiltdecon_vol.m to identify the best-fit
% location. You may also want to evaluate the stability of the moment
% tensor source time function to choose the source with both a stable
% solution and low error. I used the code evalSTFstats.m to compute the
% stability (Matoza's gamma metric) of all the solutions within 10% of the
% minimum E2 value.

%
% This version properly accounts for tilt
% It is written for GFs with 60 channels. The first 30 channels are for the
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
% tilt data following the method of Maeda et al. GJI, 2013. Note that all
% the green functions were computed using the TOPO code (topo_v17_F09_20.f)
% from Chouet and Dawson (see Ohminato and Chouet, BSSA 1998). They were
% then converted from time-domain to frequency domain using the code
% fdm2greenft.f These files have complex spectra for each station shown
% above in order following a single header line that gives sample rate,
% number of points, etc.
%%

clear
close all

% Add directory where functions to convert from lune to moment tensors and
% vice versa are stored. These are from Carl Tape
addpath('~/prgs/misc_matlab/compearth/momenttensor/matlab/')

% These are predifined filter bands used to explore the affect of tilt on
% the solutions. Most go from 10 seconds down (VLP band), but one is from
% 400 to 60 seconds, below the corner of all the stations used in 2009
%high corner of filter
HCall=[1/60, .1 .1 .1 .1 .1];
% decide with low corner filter to use. These are a good sample
LCall=1./[400 400 120  60  30 90];
tp=0.5; %disp(['using tp = ',num2str(tp)]);

% This works well for 2009 data
%chnwgtstr='17cha';
year_flag=2009;
%2009
stasub=[1:3,6:8];

%choose event 1, 2, or 3, or the phase-weighted stacks (event 4)
event=4;
eventfile=cell(5,1);
if event==1
    % this one is really good: Event 2009.1
    datastart=datenum(2009,1,19,16,9,25)-10/(24*60)-40.96/86400;
    for n=1:length(LCall)
        tmp=sprintf('event1_%03d.mat',1/LCall(n));
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
%%
dataend = datastart+20/(24*60);


plot_flag=0; % use plot_flag=2 for more plots, incl. the GFs. This is for testing only
usetilt = 1;

% if this is 1, STF.dat files will be written at each step
outputfiles_flag=0;

% choose 6 for moment components, 3 for single forces, or 9 for moment plus
% forces - ONLY TESTED with nummom=6;
nummom=6;

% This is the directory where the GFs are located. They are in individual
% directories labeled with node x-y-z (like 117-099-122) and labeled like
% greenf.1.ft and greenf.1.rot.ft for the translation and rotation GFs
% respectivley. Actually, the rotation GFs are simply the curl of the
% displacement field so the rotation is half of this.
% rootdir='/Users/gpwaite/Data/Guate/Fuego2012/FuegoVLP2012/fdmtilt';
%  if tp<1
%
%         allGFdirs=dir([rootdir,'/*-*-???']);
%     else
%         allGFdirs=dir([rootdir,'/*-*-*-tp',num2str(tp)]);
%
%  end
load('allGFdirectories.mat')
% rootdir=pwd;
rootdir='/Users/gpwaite/Data/Guate/Fuego2012/FuegoVLP2012/fdmtilt';

%make sure the data are consistent with padded GFs
td=2^15;
fd=td/2;

% this simply trims off the ends of the seismograms and synthetics for
% error calculations
timetrim=[10,645];

% 2012
%stasub=[11:15,17:20];

% determine weighting here in order to downweight or effectively eliminate
% some channels

% 27 weight
weight=ones(1,27);


% 21cha
if strcmp(chnwgtstr,'17cha')
    weight=ones(1,18);
    weight(13)=1e-20;
    weight(14)=1e-20;
    
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
    
    
    %%%%%%%%%%%%%%%%%
    % load the seismograms - raw w/o any deconvolution
    load(char(eventfile(locor)))
    sps=50;
    
    
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
    transtransmat=zeros(length(inst)*3,fd);
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
    % I put this in to ensure that all the directories had GFs in them, but
    %    % we can skip it for these 4
    %     gind=[];
    %     cnt=1;
    %     for nd=1:length(allGFdirs)
    %
    %         if exist([rootdir,'/',allGFdirs(nd).name,'/greenf.7.ft'],'file')
    %             gind(cnt)=nd;
    %             cnt=cnt+1;
    %         end
    %     end
    
    
    allGFs=allGFdirs;
    
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
    
    
    disp('results are from 22chaPWS_6moment_output_nofilt__flippedrotx2016_04_25_21_35_11_2009PWS.mat')
    disp('and refined with Matozas gamma metric');
    best_PWstack=[461,190,190,190,103,487];
    for ng=best_PWstack(locor)
        %%% otherwise, do all
        %     disp(['inverting for ',num2str(size(allGFs,1)),' nodes'])
        %     for ng=1:length(allGFs)
        Gtd=zeros(length(inst)*3,length(allmom),512); Grottd=Gtd;
        G=zeros(length(inst)*3,length(allmom),fd);Grotuncorr=G; Grot=G;
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
            gfile=sprintf('%s/%s/greenf.%d.ft',rootdir,allGFs(ng).name,nm);
            disp(['loading ',gfile]);
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
                %                 tmp=zeros(datahead{4},1,'gpuArray');
                tmp(1:datahead{1},1)=complex(datar,datai);
                tmp1g=ifft(tmp,'symmetric');
                tmp1=tmp1g*scalarmomentnormfac;
                %                 tmp1=gather(tmp1g)*scalarmomentnormfac;
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
                G(nc,nm,:)=gather(tmpfd(1:fd));
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
            disp(['loading ',gfile]);
            
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
                
                tmptd=zeros(td,1);
                tmptd(1:length(tmp1),1)=flipud(tmp1);
                
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
            fvecfd=linspace(0,sps,fd);
            constant=-9.8./(2*pi*fvecfd).^2;
            constant(1)=-1e-9;
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % invert and generate FD synthetics
        % it isn't necessary to invert over the whole range of frequencies
        % since the sample freq is 50 and the largest signals are well under 1
        % Hz. It would be quicker to restrict this to just the first 100
        % frequencies
        
        
        %%%%%%%%%%%%%%%%
        %A priori weighting:constrain the inversion at each point of the
        %lune space
        M00 = 1;
        
        %loop through all the possible lat and long on the lune
        %         delta=[0,15,30,45,60,75,90]; %latitude on the lune
        %         gamma=[-30,-22.5,-15,-7.5,0,7.5,15,22.5,30]; %longitude on the lune
        
        %         Use a  smarter gamma and delta; this requires just one
        %         loop for gamma-delta pair
        %         [gamma,delta]=getspheregrid([-30 30 0 90],5,5,1);
        [gamma,delta]=getspheregrid([-30 30 0 90],4,4,1);
        
        %                 gamma=gamma(2:25)
        %                 delta=delta(2:25)
        %loop through variations in angles
        % due to the symmetry of the symmetric tensor, we could go from
        %         % 0-180 in each, but I choose to use all dip angles and
        %         angx=0.1:30:360; %
        %         angy=0.1:30:360; %
        %         angz=0.1:30:360;
        % compute these otuside of the script in checksymmetry.m
        % load('45degangvec.mat'); %takes about  0.5 hour with 223 gamma-delta pairs
        % load('30degangvec.mat'); %takes about 1.5 hours with 223 gamma-delta pairs
        % load('25degangvec.mat');%takes about 3.5 hours with 223 gamma-delta pairs
        % load('20degangvec.mat');%takes about 6 hours with 223 gamma-delta pairs
        load('15.0degangvec.mat');%takes about 13 hours with 223 gamma-delta pairs
        % load('10degangvec.mat');%takes about 46 hours with 223 gamma-delta pairs
        %         disp(['Run will take ',num2str(totsec/86400),' days'])
        %fill these with large values, so if they get skipped, they won't
        %be have minima at zero
        a=cell2mat(ang(1));
        b=cell2mat(ang(2));
        disp(['using ',num2str(length(ang)),' rotation angles at ',num2str(b(3)-a(3)),' degrees'])
        E1=1e3*ones(length(ang),1);
        E2=E1;
        E1min=zeros(length(delta),1);
        E2min=E1min;
        % not sure if this actually works:
        storeE1=cell(length(delta),1);
        storeE2=storeE1;
        
        %move this out of the loop
        timetrimsamp=timetrim*sps;
        tts=timetrimsamp;
        
        tic
        for il=1:length(delta)
            disp([num2str( (100*(il-1))/length(delta)),'% done'])
            
            
            lam=lune2lam(gamma(il),delta(il),M00);
            
            m1=lam(1,1);
            m2=lam(2,1);
            m3=lam(3,1);
            m4=0; m5=m4; m6=m4;
            
            cost=0.00001;
            
            %             if m1 == 0 || m2 == 0 || m3 == 0 || m4 == 0 || m5 == 0 || m6 == 0
            m1=m1+cost;
            m2=m2+cost;
            m3=m3+cost;
            m4=m4+cost;
            m5=m5+cost;
            m6=m6+cost;
            %             end
            
            Min=[m1,m2,m3,m4,m5,m6]';
            M(1,1)=m1;
            M(2,2)=m2;
            M(3,3)=m3;
            M(1,2)=m4;
            M(1,3)=m5;
            M(2,3)=m6;
            M(2,1)=M(1,2);
            M(3,1)=M(1,3);
            M(3,2)=M(2,3);
            
            % loop over predetermined angle triples
            for zk=1:length(ang)
                angtmp=cell2mat(ang(zk));
                R=rotmat3(angtmp,4);
                Mout=R*M*R';
                
                
                
                % The order of components 5 & 6 is different in
                % topo
                % xx yy zz xy yz xz (11 22 33 12 23 13)
                newM=[Mout(1,1);Mout(2,2);Mout(3,3);Mout(1,2);Mout(2,3);Mout(1,3)];
                
                m1to2=(newM(1,1)*normwgt(1))/(newM(2,1)*normwgt(2));
                m1to3=(newM(1,1)*normwgt(1))/(newM(3,1)*normwgt(3));
                m1to4=(newM(1,1)*normwgt(1))/(newM(4,1)*normwgt(4));
                m1to5=(newM(1,1)*normwgt(1))/(newM(5,1)*normwgt(5));
                m1to6=(newM(1,1)*normwgt(1))/(newM(6,1)*normwgt(6));
                
                mi=[m1to2,m1to3,m1to4,m1to5,m1to6];
                
                %disp(['input moment ratios: ',num2str([Mout(1,1)/Mout(2,1), Mout(1,1)/Mout(3,1), Mout(1,1)/Mout(4,1), Mout(1,1)/Mout(5,1), Mout(1,1)/Mout(6,1)])])
                %disp(['input norm ratios: ',num2str(mi)])
                
                H = zeros(5,6);
                H(:,1)=ones(1,1);
                
                for n=1:length(mi)
                    H(n,n+1)=-mi(n);
                end
                
                h=zeros(5,1);
                
                F=zeros(11,11);
                
                % m=zeros(nummom,fd);
                m=zeros(11,fd);
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
                        F(1:6,1:6)=tmpG'*tmpG;
                        F(7:11,1:6)=H;
                        F(1:6,7:11)=H';
                        f=[tmpG'*(tmpdata*stf(nf));h];
                        m(:,nf)=F\f;
                        
                        %m(:,nf)=tmpG\(tmpdata*stf(nf)); %#ok<*SAGROW>
                    end
                    clear tmpG
                end
                %disp(['output norm ratios: ',num2str([real(m(1,2)/m(2,2)), real(m(1,2)/m(3,2)),real(m(1,2)/m(4,2)),real(m(1,2)/m(5,2)),real(m(1,2)/m(6,2))])])
                
                % un-norm the moments
                for nm=1:nummom
                    m(nm,:)=m(nm,:)/normwgt(nm);
                end
                
                synthfd=zeros(size(data,1),fd);
                if usetilt
                    synthfd=zeros(size(data,1),fd);
                end
                
                % skip the frequencies above 1 Hz
                % for nf=1:fd
                for nf=1:657
                    
                    %         tmpG(:,:)=Gorig(chasub,:,nf);
                    tmpG(:,:)=G(chasub,:,nf);
                    % we need to divide by the stf since the GF is convolved with it
                    % moment components are multiplied by the weighting vector used to
                    % invert
                    %         if ~usetilt
                    synthfd(:,nf)=tmpG*(m(1:6,nf).*normwgt(:))/stf(nf);
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
                synthtd=zeros(length(sta)*3,td);
                for nc=1:size(synthfd,1)
                    tmp=zeros(1,td);
                    tmp(1:fd)=synthfd(nc,:);
                    synthtd(nc,:)=ifft(tmp,'symmetric');
                end
                
                sysub=synthtd(:,tts(1):tts(2));
                    dasub= datatd(:,tts(1):tts(2));
                
                num1=zeros(1,21);
                num2=zeros(1,3);
                den1=zeros(1,21);
                den2=zeros(1,3);
                E2tmp=zeros(1,length(sta));
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
                
                E1(zk)=sum(num1)/sum(den1);
                E2(zk)=sum(E2tmp)/length(stasub);
                
            end %end of loop of ang (zk)
            elapsedtime=toc;
            disp(['elapsed time: ',num2str(elapsedtime/3600),' hours'])
            disp(['total estimated time: ',num2str((length(delta)/il)*(elapsedtime/3600)),' hours'])
            % change this to save everything
            storeE1(il)={E1};
            storeE2(il)={E2};
            E1min(il)=min(E1(:));
            E2min(il)=min(E2(:));
            
        end%end of loop for il (delta)
        
        %         E1all(locor,ng)=E1(ng);
        %         E2all(locor,ng)=E2(ng);
        
        if ~mod((100*ng/length(E1)),10)
            disp([num2str(100*ng/length(E1)),'% done with LC ',num2str(1/LC),' HC ',num2str(1/HC)])
        end
        
        if outputfiles_flag
            stfoutfile=sprintf('%s/%s/STF_%03d.dat',rootdir,allGFs(ng).name,round(1/LC));
            fid=fopen(stfoutfile,'w');
            fprintf(fid,'%d %d %f %d\n',td,size(m,1),1/sps,fd);
            mtd=zeros(length(allmom),td);
            for nm=1:size(m,1)
                tmpm=zeros(1,td);
                tmpm(1:fd)=m(nm,:);
                mtd(nm,:)=ifft(tmpm,'symmetric');
                fprintf(fid,'%f\n',mtd(nm,:));
            end
            fclose(fid);

        end
        
    end
    % save some output
    dte = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    
    savefile=sprintf('Full_lune_output_%d_%d_PWS_%s.mat',1/LC,1/HC,dte);
    
    save(savefile,'allGFs','E1min','E2min','storeE1','storeE2','LC','HC','gamma','delta','ang');
    
    
end
%if event==1
%elseif event==2
%    savefile=sprintf('%s_%dmoment_lune_%s_Event2009.2.mat',chnwgtstr,nummom,dte);
%elseif event==3
%    savefile=sprintf('%s_%dmoment_lune_%s_Event2009.3.mat',chnwgtstr,nummom,dte);
%end
if 1
end
if plot_flag
    figure;
    subplot(211)
    plot(E1min'); title('E1')
    legend(num2str(1./LCall'))
    
    hold on
    [a,b]=min(E1min,[],2);
    plot(b,a,'*')
    subplot(212)
    plot(E2min'); title('E2')
    legend(num2str(1./LCall'))
    
    hold on
    [a,b]=min(E2min,[],2);
    plot(b,a,'*')
    
end
stfdir=sprintf('%s/%s',rootdir,allGFs(ng).name);


%% plot residuals
if plot_flag
    tit=[datestr(datastart+10/(24*60)+40.96/86400),' ',allGFs(ng).name];
    %     for nt=1:size(combdata,1)
    %     tmp1=zeros(1,td);
    %     tmp1(1:fd)=combdata(nt,:);
    %     dataraw(nt,:)=ifft(tmp1,'symmetric');
    %     end
    if usetilt
        
        plotresiduals(datatd,synthtd(1:18,:),weight,sta,tit,270,380,1/sps,0);
        print('-dpsc',['res_vel',num2str(sum(weight)),'cha_',num2str(round(1/LC)),'s_',allGFs(ng).name,'.ps'])
        
        %integrate the data
        for nt=1:size(datatd,1)
            intdatatd(nt,:) =cumtrapz(datatd(nt,:))/50;
            intsynthtd(nt,:)=cumtrapz(synthtd(nt,:))/50;
        end
        plotresiduals(intdatatd,intsynthtd,weight,sta,tit,1,655,1/sps,0);
        %         print('-dpsc','residualsint.ps')
        print('-dpsc',['res_disp',num2str(sum(weight)),'cha_',num2str(round(1/LC)),'s_',allGFs(ng).name,'.ps']);
        
        tmp=pwd;
        cd (stfdir)
        
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



%%
a1=cell2mat(ang(1));
a2=cell2mat(ang(2));
dz=a2(3)-a1(3);
sorted=sortrows([gamma,delta,E2min]);
[x,y]=meshgrid(-30:1:30,0:1.5:90);
int=griddata(gamma,delta,E2min,x,y,'linear');
figure;imagesc(x(1,:),y(:,1),int);axis image
colorbar
axis xy
tittxt=sprintf('BP%d-%d PWS phi-theta-%d',1/LC,1/HC,dz);
lunefname=sprintf('LunePlotBP%d-%d_PWS_phi-theta-%d.ps',1/LC,1/HC,dz);

title(tittxt)
print('-dpsc',lunefname)
