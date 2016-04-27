function [gamma,ff]=stfgammacalc(infile,peramp,ratio)
% stfgammacalc(infile,peramp,ratio,titl)
% infile is the name of the STF.dat file, peramp is the fraction of the
% maximum used as a cut-off for the point-by-point eigenvectors, 

if nargin==1
% set this to some fraction of the peak amplitude that you want to use to
% avoid including statistics from low amplitude, noisey points
 peramp=0.60;
 ratio=2;
end
% get information from the header

fi=fopen(infile);
a=fscanf(fi,'%d %d %f %d',[1 4]);
np=a(1);nt=a(2);dt=round(1000*a(3))/1000;
d=fscanf(fi,'%f',[nt*np 1]);
fclose(fi);
t=(1:np)*dt;
dat=zeros(nt,np);

rmaxx=dt*(np-1);
cut1=1;cut2=t(end)-1;
maxx=cut2;
st=cut1/dt+1;
en=cut2/dt+1;

% integrate so plotting moment instead of moment rate
for i=1:nt;for j=1:np;dat(i,j)=d((i-1)*np+j);end;end
for i=1:nt;dat(i,:)=cumtrapz(t,dat(i,:));end


Mij=dat;

% deg=180/pi;
[a,~] = max(abs(Mij),[],2);
[maxtraceamp,maxamptraceNo] = max((a));
ampD=maxtraceamp;

% use ampD to constrain plotting vectors with some meaningful amplitude
ampD=ampD*peramp;
cnt=0;
%
for i=cut1/dt:1:cut2/dt % 10 - 140 s
    % for i=130/dt:1:180/dt %:length(Mij)
    if abs(Mij(maxamptraceNo,i)) > ampD
        for j=1:6
            p(j)=Mij(j,i);
        end
        %  set the matrix
        for j=1:3;m(j,j)=p(j); end
        m(1,2)=p(4);m(2,1)=p(4);
        m(2,3)=p(5);m(3,2)=p(5);
        m(1,3)=p(6);m(3,1)=p(6);
        %  cal eigen values and vectors
        [~, D]=eig(m);
        cnt=cnt+1;
        ll(cnt)=i;
        ee=abs(D)/max(max(abs(D)))*ratio;
        maxD=ijkabs(D);
        jmax=maxD(1,1);
        if jmax == 3
%             disp ('jmax = 3')
            ff(1,cnt)=ee(1,1);
            ff(2,cnt)=ee(2,2);
            ff(3,cnt)=ee(3,3);
            
        else
%             disp ('jmax = 1')
            ff(1,cnt)=ee(3,3);
            ff(2,cnt)=ee(2,2);
            ff(3,cnt)=ee(1,1);
        end
        % check for eigenvalues with the opposite sign
        if ff(1,cnt)<0, error('negative miminum eigenvalue'); end
        if ff(2,cnt)<0, error('negative intermediate eigenvalue'); end

    end
    
end

% ff(1,:) is a vector of the ratios between the minimum and maximum
% (normalized to ratio). f(2,:) is a vector of the ratios between the
% intermediate and maximum. So following the equation (6) in Matoza et al.,
% JGR 2015, we should multiply these vectors by the ratio before computing
% the standard deviation
s11=std(ratio*ff(1,:));
s22=std(ratio*ff(2,:));
% s33=std(ff(3,:)); %%% this is always 0 since ff(3) is always the ratio value
% gamma=sqrt(s11^2+s22^2+s33^2); % this is what Robin told me he did
gamma=sqrt(s11^2+s22^2); %but this makes more sense