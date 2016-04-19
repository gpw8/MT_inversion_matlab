% gather results from grid searches
clear, close all

for band=4%[1:3,5]%:5;

if band==1
load('/Users/gpwaite/Data/Guate/Fuego2009/LuneGrid/rotxr/Full_lune_output_400_60_PWS_2016_04_17_02_23_02.mat')
    % now grid the data
    sorted=sortrows([gamma,delta,E2min]);
    [x,y]=meshgrid(-30:1:30,0:1.5:90);
    int=griddata(gamma,delta,E2min,x,y,'linear');
    figure;imagesc(x(1,:),y(:,1),int);axis image
    colorbar
    axis xy
    tittxt=sprintf('BP%d-%d PWS',1/LC,1/HC);
    lunefname=sprintf('LunePlotBP%d-%d_PWS_.ps',1/LC,1/HC);
    
    title(tittxt)
%     print('-dpsc',lunefname)

    % create grdfile
    outfile=sprintf('lune_%03d_%02d_int1_E2_grid.grd',1/LC,1/HC);
    grdwrite2(x(1,:),y(:,1),int,outfile)

elseif band==2
%     LC=1/400;
%     HC=1/10;
    % now grid the data
    sorted=sortrows([gamma,delta,E2min]);
    [x,y]=meshgrid(-30:1:30,0:1.5:90);
    int=griddata(gamma,delta,E2min,x,y,'linear');
    figure;imagesc(x(1,:),y(:,1),int);axis image
    colorbar
    axis xy
    tittxt=sprintf('BP%d-%d PWS',1/LC,1/HC);
    lunefname=sprintf('LunePlotBP%d-%d_PWS_.ps',1/LC,1/HC);
    
    title(tittxt)
    %print('-dpsc',lunefname)
    
    % create grdfile
    outfile=sprintf('lune_%03d_%02d_int1_E2_grid.grd',1/LC,1/HC);
    grdwrite2(x(1,:),y(:,1),int,outfile)

elseif band==3
    load('Full_lune_output_120_10_PWS_2016_04_17_02_00_22.mat')
% now grid the data
    sorted=sortrows([gamma,delta,E2min]);
    [x,y]=meshgrid(-30:1:30,0:1.5:90);
    int=griddata(gamma,delta,E2min,x,y,'linear');
    figure;imagesc(x(1,:),y(:,1),int);axis image
    colorbar
    axis xy
    tittxt=sprintf('BP%d-%d PWS',1/LC,1/HC);
    lunefname=sprintf('LunePlotBP%d-%d_PWS_.ps',1/LC,1/HC);
    
    title(tittxt)
    %print('-dpsc',lunefname)
    
    % create grdfile
    outfile=sprintf('lune_%03d_%02d_int1_E2_grid.grd',1/LC,1/HC);
    grdwrite2(x(1,:),y(:,1),int,outfile)


elseif band==4
load('Full_lune_output_60_10_PWS_2016_04_17_02_47_32.mat')
    % now grid the data
    sorted=sortrows([gamma,delta,E2min]);
    [x,y]=meshgrid(-30:1:30,0:1.5:90);
    int=griddata(gamma,delta,E2min,x,y,'linear');
    figure;imagesc(x(1,:),y(:,1),int);axis image
    colorbar
    axis xy
    tittxt=sprintf('BP%d-%d PWS',1/LC,1/HC);
    lunefname=sprintf('LunePlotBP%d-%d_PWS_.ps',1/LC,1/HC);
    
    title(tittxt)
    %print('-dpsc',lunefname)
    
    % create grdfile
    outfile=sprintf('lune_%03d_%02d_int1_E2_grid.grd',1/LC,1/HC);
    grdwrite2(x(1,:),y(:,1),int,outfile)


elseif band==5
load('Full_lune_output_30_10_PWS_2016_04_17_02_46_38.mat')
    % now grid the data
    sorted=sortrows([gamma,delta,E2min]);
    [x,y]=meshgrid(-30:1:30,0:1.5:90);
    int=griddata(gamma,delta,E2min,x,y,'linear');
    figure;imagesc(x(1,:),y(:,1),int);axis image
    colorbar
    axis xy
    tittxt=sprintf('BP%d-%d PWS',1/LC,1/HC);
    lunefname=sprintf('LunePlotBP%d-%d_PWS_.ps',1/LC,1/HC);
    
    title(tittxt)
    %print('-dpsc',lunefname)
    
    % create grdfile
    outfile=sprintf('lune_%03d_%02d_int1_E2_grid.grd',1/LC,1/HC);
    grdwrite2(x(1,:),y(:,1),int,outfile)


end
end