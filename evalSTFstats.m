% evaluate gamma of low E2 models
clear, close all

% choose a cut-off E2 to evaluate
cutoff=.1;

% choose a normalization ratio
normrat=2;
% first load a file with E2 results from all the bands
% this one has the 6 event PWS
% load 22chaPWS_6moment_output_expandedZ_2016_03_19_19_19_11_2009PWS.mat
% load('22chaPWS_6moment_output_expandedZ_2016_03_28_13_55_50_2009PWS.mat')
% load('22chaPWS_6moment_output_expandedZ__flippedrotx2016_04_08_16_20_49_2009PWS.mat')
% this one only has the new 90-10 in band 6
% load('22chaPWS_6moment_output_expandedZ__flippedrotx2016_04_15_16_53_01_2009PWS.mat')
%

% this one has bands 1-6 (400 400 120 60 30 90; 60 10 10 10 10 10)
% load 22chaPWS_6moment_output_expandedZ__flippedrotx2016_04_18_18_28_14_2009PWS.mat

load 22chaPWS_6moment_output_nofilt__flippedrotx2016_04_25_21_35_11_2009PWS.mat


for band=1:6;
    disp(['BP ',num2str(1/LCall(band)),'-',num2str(1/HCall(band))])
    %
    rootdir='/Users/gpwaite/Data/Guate/Fuego2012/FuegoVLP2012/fdmtilt';
    
    % get the min
    [minE2,b]=min(E2all(band,:));
    disp(['min E2 solution is ',num2str(minE2),' at ',allGFs(b).name,' at index ',num2str(b)])
    
    
    [~,d]=find(E2all(band,:)<=minE2*(1+cutoff));
    
    % now loop over all these to calculate the stats
    
    for i=1:length(d)
        stfdir=sprintf('%s/%s',rootdir,allGFs(d(i)).name);
        stffile=sprintf('%s/STF_%03d-%03d.dat',stfdir,1/LCall(band),1/HCall(band));
        [gam(i),eigvals]=stfgammacalc(stffile,.5,normrat);
        disp(['  ',stffile(end-26:end-4),',  gam: ',num2str(gam(i)),',  E2: ',num2str(E2all(band,d(i)))]);
        if gam(i)<0.615
            stfmfilename(stffile,1,1,655,['BP ',num2str(1/LCall(band)),'-',num2str(1/HCall(band)),' ',allGFs(d(i)).name,' gamma: ',num2str(gam(i))])
            % plot ratio histograms
%             figure
            subplot(222);
            m=mean(eigvals(1,:)); bincent=m-.5:.1:m+.5;
            hist(eigvals(1,:),bincent); xlabel('minimum eigenvalue')
            subplot(224);
            m=mean(eigvals(2,:)); bincent=m-.5:.1:m+.5;
            hist(eigvals(2,:),bincent);xlabel('intermediate eigenvalue')
%             suptitle2([stffile(end-26:end-4),' normalized to ',num2str(normrat)]);
        end
    end
    [e,f]=min(gam);
    disp(['best solution is ',allGFs(d(f)).name,' at index ',num2str(d(f)),' with Matoza-gamma =',num2str(e)])
    disp(['with E2: ',num2str(E2all(band,d(f)))])
    disp('  ')
    
        pause
        close all
        clear gam b c d e f minE2
    
    
end