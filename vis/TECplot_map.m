addpath ../script_utils;

%SIMULATIONS LOCAITONS
simname='tohoku20113D_highres_long/';
basedir='~/zettergmdata/simulations/'
direc=[basedir,simname];
system(['mkdir ',direc,'/TECplots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/TECplots_eps']);    %store output plots with the simulation data


%LOAD THE COMPUTED MAGNETIC FIELD DATA
load([direc,'/vTEC_hemis.mat']);
lt=numel(t);
mlon=mlong;


%SIMULATION META-DATA
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs/config.dat']);


%TABULATE THE SOURCE LOCATION
mlatsrc=mloc(1);
mlonsrc=mloc(2);
thdist=pi/2-mlatsrc*pi/180;    %zenith angle of source location
phidist=mlonsrc*pi/180;


figure(1);
%set(gcf,'PaperPosition',[0 0 4 8]);


%MAKE THE PLOTS AND SAVE TO A FILE
for it=1:lt
    fprintf('Printing TEC plots...\n');
    %CREATE A MAP AXIS
    figure(1);
    clf;
%    FS=16;
    FS=10;
    
    datehere=simdate_series(it,:);
    ymd=datehere(1:3);
    UTsec=datehere(4)*3600+datehere(5)*60+datehere(6);
    filename=datelab(ymd,UTsec);
    filename=[filename,'.dat'];
    titlestring=datestr(datenum(datehere));
    
    mlatlimplot=[min(mlat)-0.5,max(mlat)+0.5];
    mlonlimplot=[min(mlon)-0.5,max(mlon)+0.5];
    axesm('MapProjection','Mercator','MapLatLimit',mlatlimplot,'MapLonLimit',mlonlimplot);
    param=dvTEC(:,:,it);
    %imagesc(mlon,mlat,param);
    mlatlim=[min(mlat),max(mlat)];
    mlonlim=[min(mlon),max(mlon)];
    [MLAT,MLON]=meshgrat(mlatlim,mlonlim,size(param));
    pcolorm(MLAT,MLON,param);
    colormap(parula(256));
    set(gca,'FontSize',FS);
    tightmap;
%    caxis([-3,3]);
    caxis([-0.25,0.25]);
    c=colorbar
    set(c,'FontSize',FS)
    xlabel(c,'\Delta vTEC (TECU)')
    xlabel(sprintf('magnetic long. (deg.)\n\n'))
    ylabel(sprintf('magnetic lat. (deg.)\n\n\n'))
    hold on;
    ax=axis;
    plotm(mlatsrc,mlonsrc,'r^','MarkerSize',10,'LineWidth',2);
    hold off;
    titlestring=datestr(datenum(simdate_series(it,:)));
    title(sprintf([titlestring,'\n\n']));
    %gridm;
    
    
    %ADD A MAP OF COASTLINES
    load coastlines;
    [thetacoast,phicoast]=geog2geomag(coastlat,coastlon);
    mlatcoast=90-thetacoast*180/pi;
    mloncoast=phicoast*180/pi;
    
    if (360-mlonsrc<20)
        inds=find(mloncoast>180);
        mloncoast(inds)=mloncoast(inds)-360;
    end
    
    hold on;
    ax=axis;
    plotm(mlatcoast,mloncoast,'k-','LineWidth',1);
    hold off;
    setm(gca,'MeridianLabel','on','ParallelLabel','on','MLineLocation',10,'PLineLocation',5,'MLabelLocation',10,'PLabelLocation',5);
%    setm(gca,'MeridianLabel','on','ParallelLabel','on','MLineLocation',1,'PLineLocation',1,'MLabelLocation',1,'PLabelLocation',1);
    gridm on;
    
    
    %PRINT THE THING
    print('-dpng',[direc,'/TECplots/',filename,'.png'],'-r300');
%    print('-depsc2',[direc,'/TECplots_eps/',filename,'.eps']);
end

rmpath ../script_utils;
