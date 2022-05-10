%% Code uses a digital elevation model to compute glacier thickness based 
%  on the Glabtop2 approach (reference below) presented in Linsbauer et al
% (2015). Code was initially written by S. Ragettli, and has been extensively 
% adapted by E. Miles, but precisely follows the Glabtop2 approach.
     % AUTHOR    : Evan Miles (evan.miles@wsl.ch)
     % $DATE     : 20-Sep-2018
     % DEVELOPED : R2016a
     % FILENAME  : Glabtop2_ESM.m
%
% Linsbauer, A., Frey, H., Haeberli, W., Machguth, H., Azam, M. F., & 
% Allen, S. (2015). Modelling glacier-bed overdeepenings and possible 
% future lakes for the glaciers in the Himalaya�Karakoram region. Annals 
% of Glaciology, 57(71), 119�130. https://doi.org/10.3189/2016AoG71A627

clear all
close all

%UPDATE/SET WORKING DIRECTORY
%workdir = 'E:\Teaching\GEOG3669 - The Cryosphere\201819\Practical 2 GlabTop';

%UPDATE/SET DEM USED FOR INPUT; must be geotiff, in working directory,
%include file extension; do not use geographic projection (use UTM or
%similar projection so that pixels are equal area; be sure to clip to the
%relevant area
DEM_file = 'N028E085_AVE_DSM_UTM45_crop.tif';

%UPDATE/SET DEM USED FOR INPUT; must be shapefile, in working directory, in
%same projection as DEM; be sure to only include outlines for which you
%want results, can include multiple outlines
OUTLINES_file = 'Langtang_glacier_outlines_23092015_DG2014Jan';

N = 5; %number of times to run the code - output will be the mean thickness

%UPDATE/SET OUTPUT NAME, with file extension (.tif)
output = ['GlabTop2Output_' date '.tif'];


%% load DEM
%cd(workdir)
[DEM,Dx,Dy,DI] = geoimread(fullfile(DEM_file));

[Xv,Yv] = meshgrid(Dx,Dy);
Xv=Xv(:);Yv=Yv(:);

%% create glacier mask from shapefile of outlines
    gmaskfile = fullfile(regexprep(DEM_file,'.tif','_mask.tif'));

    OUTLINES = shaperead(OUTLINES_file);
    AP = repmat(0*DEM,[1,1,length(OUTLINES)]);

    %convert glacier outlines to mask for entirety of DEM
    for iP = 1:12%length(glaSHP)
    %     curM = mask;
        curM = inpolygon(Xv,Yv,OUTLINES(iP).X,OUTLINES(iP).Y);
        curM = reshape(curM,size(DEM));
        AP(:,:,iP) = int8(curM).*iP;
    end

    GLAC = nansum(AP,3);
    geotiffwrite(gmaskfile,GLAC,DI.RefMatrix,'CoordRefSysCode',DI.GeoTIFFCodes.PCS)

%% smooth DEM by 100m Gaussian filter, calculate slope
DEM = double(DEM);
DEM(DEM==0)=NaN; %remove unrealistic values
dx = mode(diff(Dx)); %determine pixel spacing (assumes dx = dy)
pixdist100 = ceil(100/dx); %determine how many pixels in 100 m
DEM2 = imgaussfilt(DEM,pixdist100); %filter DEM to remove noise at the 100-m distance
SLO = imgradient(DEM2)./(dx.^2); %percent slope
SLO(SLO<tan(2*pi/180))=tan(2*pi/180); %limit lower-bound to 2-degrees; Glabtop2 produces unrealistic estimates for lower slope values
SLOa = atan(SLO); %slope in radians
SLOd = SLOa*180/pi; %slope in degrees
% SLO2 = imgaussfilt(SLO,3);
figure
imagesc(SLOd);colorbar
title('Surface slope in degrees')

%% calculate basal shear stress for each glacier: compute delta H
%this is directly from Linsbauer et al (2015)
glaid=unique(GLAC(GLAC>0));
for ig=1:length(unique(GLAC(GLAC>0)))
    deltaH(ig)=(max(DEM(GLAC==glaid(ig)))-min(DEM(GLAC==glaid(ig))))/1000;
    if deltaH(ig)>1.6
        tau(ig)=150000;
    else
        tau(ig)=min([0.005+1.598*deltaH(ig)-0.435*deltaH(ig)^2 1.5])*100000;
    end
end

%% determine buffer cells: 
%inner buffer (cells at inside margin of outlines)
slo_buffer=SLO;
slo_buffer(GLAC>0)=nan;
slo_buffer=ndnanfilter(slo_buffer,'rectwin', [1 1]);
slo_buffer(isnan(GLAC)==1)=nan;
slo_buffer_in=slo_buffer;
slo_buffer_in(slo_buffer_in>-9999)=2;
% figure
% imagesc(slo_buffer_in)
% impixelinfo

%outer buffer (cells at outside margin of outlines)
slo_buffer=SLO;
slo_buffer(isnan(GLAC)==1)=nan;
slo_buffer=ndnanfilter(slo_buffer,'rectwin', [1 1]);
slo_buffer(GLAC>0)=nan;
slo_buffer_out=slo_buffer;
slo_buffer_out(slo_buffer_out>-9999)=1;
slo_buffer=slo_buffer_in;
slo_buffer(slo_buffer_out==1)=1;
slo_buffer(isnan(slo_buffer)==1&GLAC>0)=3;

% figure
% imagesc(slo_buffer)
% caxis([0 3])
% impixelinfo

slo_gla=SLOa; %calculations are in radians
slo_gla(isnan(GLAC)==1)=nan;

%% calculate ice thickness for sample points
%varying window size slope calculation as in Linsbauer et al (2015)
% pixdist990 = ceil(990/dx); %determine how many pixels in 990 m
% igmax = ceil((pixdist990-1)/2);%calculates npixels for a 990 m max window size (as in Linsbauer et al (2015))
igmax = 35;
for ig=1:igmax%i=5 corresponds to a window of 2*5+1=>11x11 cells
    slo_b(ig).gla= ndnanfilter(slo_gla,'rectwin', [ig ig]);
    slo_b(ig).gla(isnan(GLAC)==1)=nan;
end

dem_gla=DEM2;
dem_gla(isnan(GLAC)==1)=nan;

%choose dh value that must be exceeded within the window (ie to resolve
%longitudual stress reasonably well)
hga=20;
% hga=15;
%hga=50;  %as in Frey 2014

%determine sampling rate
pctsamp = 0.3*(dx./90).^2;%Frey used 0.3 for 3-arc-sec DEM (~90m), reduce by ratio of pixel sizes squared

for n=1:N %iterate N times
    disp(n)
Vint=slo_buffer;
Vint(slo_buffer==1)=0;
Vint(slo_buffer~=1)=nan;
for ig=1:length(unique(GLAC(GLAC>0)))
    sample(ig).celid=find(GLAC==glaid(ig)&slo_buffer==3);
    ncells(ig)=length(GLAC(GLAC==glaid(ig)&slo_buffer==3));
    ncells_r(ig)=ceil(ncells(ig)*pctsamp);
    sample(ig).cells=sample(ig).celid(randsample(ncells(ig),ncells_r(ig)));
    for jc=1:length(sample(ig).cells)
        [xx,yy]=ind2sub(size(GLAC),sample(ig).cells(jc));
%         [xx yy]=find(cel==sample(ig).celid(sample(ig).cells(j)));
        b=0;
        dh=0;
        while dh<hga
            b=b+1;
            dhmin=nanmin(nanmin((dem_gla(xx-b:xx+b,yy-b:yy+b))));
            dhmax=nanmax(nanmin((dem_gla(xx-b:xx+b,yy-b:yy+b))));
            dh=dhmax-dhmin;
        end
        alpha=atan(slo_b(b).gla(xx,yy));%* pi / 180;
%         alpha=atan(slo_b(b).gla(cel==sample(ig).celid(sample(ig).cells(jc))));%* pi / 180;
%         alpha=slo_b(b).gla(cel==sample(i).celid(sample(i).cells(j)))* pi / 180;

        Vint(xx,yy)=tau(ig)./(0.8*914*9.81*sin(alpha));
%         glab(n).glah(cel==sample(ig).celid(sample(ig).cells(jc)))=tau(ig)./(0.8*890*9.81*sin(alpha));
        
    end
%     glab(:,:,n) = Vint;
end
    

    % [x y]=find(isnan(glab(n).glah)==1&slo_buffer>1);
    [x] = 1:size(DEM,2);
    [y] = [1:size(DEM,1)]';
    ic = find(Vint>0);
    [yc xc]= ind2sub(size(Vint),ic);
    vc=Vint(ic);
    Vint2=griddata(xc,yc,vc,x,y,'cubic');
%     Vint2(Vint2==hga)=NaN;
    Vint2(slo_buffer==1)=nan;
    glab(:,:,n)=Vint2;
    
    
% figure
% imagesc(glab(:,:,n));colorbar
clear Vint xc yc ic vc
end 
% impixelinfo


%% calculate mean and export
glahmean = nanmean(glab,3);
glahmedian = nanmedian(glab,3);


figure
imagesc(glahmean)
impixelinfo
title('Mean estimated glacier thickness')
% caxis([0 400])
colorbar

figure
imagesc(glahmedian)
impixelinfo
% caxis([0 400])
title('Median estimated glacier thickness')
colorbar

geotiffwrite(output,glahmean,DI.RefMatrix,'CoordRefSysCode',DI.GeoTIFFCodes.PCS)


% glahglab(glahglab==-9999)=nan;
% figure
% contourf(glahmean,20,'LineColor','none')
% % caxis([0 100])
% caxis([0 500])
% set(gca,'YDir','reverse');
% colorbar

% save('ice_thickness')

