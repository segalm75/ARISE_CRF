%% Details of the function:
%  NAME:
% readAMSR2data
%----------------------------------
% PURPOSE:
%  - read AMSR2 sea ice concentration hdf files
%
% CALLING SEQUENCE:
%  d = readAMSR2data
%
% INPUT:
% daystr = yyyymmdd - date to process
% example: {'20140901' '20140902'}
% 
% OUTPUT:
%  plots of sea ice concentration, MIZ edge line, sea ice conc matrix
%  sea ice mask (yes/no)
%
%
% DEPENDENCIES:
%  - startup_plotting.m
%  - save_fig.m
% 
%
% NEEDED FILES/INPUT:
% 
%  - AMSR2 daily .hdf files of sea ice concentation
%  - AMSR2 Longitude/Latitude .hdf file
%  - e.g. 'LongitudeLatitudeGrid-n3125-NorthWestPassage.hdf'
%  - AMSR data are from:http://www.iup.uni-bremen.de:8084/amsr2data/asi_daygrid_swath/n3125/2014/sep/
%
% EXAMPLE:
%  - out = readAMSR2data
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Feb-10-2015
% -------------------------------------------------------------------------
%% function routine
function out = readAMSR2data(daystr)

startup_plotting; 
movie = false;
cplots= false;
plotting = false;
%% load data
% load AMSR2
% this is for NorthWestPassage
adir     = 'F:\ARISE\GoogleEarth\AMSR\AMSR2\';
gridfile = strcat(adir,'LongitudeLatitudeGrid-n3125-NorthWestPassage.hdf');

% prepare variable for saving a movie
if movie
    mov(length(daystr)) = struct('cdata',[],'colormap',[]);
end
for i=1:length(daystr)
    datafile = strcat(adir,'asi-AMSR2-n3125-',daystr{:,i},'.hdf');

    d.iceconc = double(hdfread(datafile,'/ASI Ice Concentration')); 
    d.longrid = double(hdfread(gridfile,'/Longitudes')); 
    d.latgrid = double(hdfread(gridfile,'/Latitudes')); 
    dx = diff(d.longrid); dx = [dx;dx(end,:)];
    dy = diff(d.latgrid); dy = [dy;dy(end,:)];

    % prepare for plotting sea ice conc
    if plotting
        figure(1);set(gcf, 'Color','white');
        %POLARSTEREO_FWD transforms lat/lon data to map coordinates for a polar stereographic system
        %[X,Y]=POLARSTEREO_FWD(LAT,LONG,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY) 
        axis tight manual
        % getm('MapProjection') for the various projections
        ax = gca;
        %ax = axesm('MapProjection','mercator');
        ax.NextPlot = 'replaceChildren';
        hold all;
        colormap cool
        surf((d.longrid),(d.latgrid),(d.iceconc),'FaceColor','interp',...
             'EdgeColor','none');
        caxis([0 100]);
        axis([-180 -100 70 85 0 100]);grid off;
        set(gca,'Ytick',70:2.5:85);
        set(gca,'YTickLabel',{'70.0','72.5','75.0',...
            '77.5','80.0','82.5','85.0'});
        xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12);
        title(daystr{:,i});
        cb=colorbarlabeled('Sea Ice [%]');
        set(gca,'FontSize',12);
        view(2)
        
        % save figure
        
        fi=[strcat(adir, daystr{:,i}, 'seaIceConc')];
        save_fig(1,fi,false);
    end
    
    if movie
        mov(i) = getframe(gcf);
    end
    
    
    % create gradient
    %[px,py] = gradient(d.iceconc,d.longrid,d.longrid);
    [d.icegradx,d.icegrady] = gradient(d.iceconc);
    
    % create gradient line
    
    % contour plots
    if cplots
        figure;
        ax1=contour(d.longrid,d.latgrid,d.icegradx);title('gradient along x');
        axis([-180 -100 70 85]);
        caxis([-50 50]);colorbar;
        figure;
        ax2=contour(d.longrid,d.latgrid,d.icegrady);title('gradient along y');
        axis([-180 -100 70 85]);
        caxis([-50 50]);colorbar;
        figure;
        ax3=contour(d.longrid,d.latgrid,d.icegradx+d.icegrady);title('gradient along xy');
        axis([-180 -100 70 85]);
        caxis([-50 50]);colorbar;
        %contour(d.longrid,d.latgrid,d.iceconc), hold on, 
        %quiver(d.longrid,d.latgrid,px,py),      hold off
    end
    out.(strcat('amsr',daystr{:,i})) = d;
    clear d;
end
%% save processed struct
si = [adir,'AMSR2_',daystr{:,1},'_',daystr{:,end},'.mat'];
disp(['saving to ' si]);
save(si,'-struct','out');
return;