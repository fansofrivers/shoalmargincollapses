function ZZ2=ShoalmarginCollapses(directory,ID1,i,Ec,Ec2,FSavg,Lsm,GridSize,GridSize2)
%%
% This Matlab-script is created by WM van Dijk (26 October 2018)
% The code below includes Matlab scripts associated to Delft3D
% The following lines determine the location where collapses could occur
% based on analysis of Van Dijk et al (2018) in Earth Sruface Processes &
% Landforms.
% 
load('E:\MedianElevation.mat') %load median elevation of the Western Scheldt
%% Open trim-file
% Open the data from the last timestep
NFS = vs_use(strcat(directory, 'work\','trim-', ID1, '_', num2str(i-1),'.dat'),'quiet'); % run trim file that changes over time

X = vs_get(NFS,'map-const','XCOR','quiet'); % get X values of the grid
Y = vs_get(NFS,'map-const','YCOR','quiet'); % get Y values of the grid
depth = vs_get(NFS,'map-sed-series','DPS','quiet'); % read all depth files
depth0 = vs_get(NFS,'map-sed-series',{1},'DPS','quiet'); % get depth file at start
depth1 = vs_get(NFS,'map-sed-series',{length(vs_get(NFS,'map-series','S1','quiet'))},'DPS','quiet'); % get the final depth

dZ=-depth1--depth0; % calculate erosion-deposition during simulation
depth = depth{end}; %   get the depth of the last timestep
dX = abs(X(3,3) - X(2,2)); % Grid spacing X-direction
dY = abs(Y(3,3) - Y(2,2)); % Grid spacing Y-direction
z = -depth;
z(z==0) = NaN;

% Distribution of the Shoal margin collapse sizes and volumes based on analysis for the Western Scheldt by Van Dijk et al (2017)
xa = (10000:1000:300000)'; % size of the collapes
ya = lognpdf(xa,10.38,0.88); % log-normal distribution (10.38 and 0.88 are site specific)
xv = (10000:10000:3000000)'; % volume of the collapses
yv = lognpdf(xv,11.59,1.21); % log-normal distribution (11.59 and 1.21 are site specific)
ya = ya./sum(ya);
yv = yv./sum(yv);
yac = cumsum(ya); % cumulative size
yav = cumsum(yv); % cumalative volume
%% Calculate the relative slope height (H), gradient (Gr), Aspect (A) and Probability (P)
Hr=NaN(size(z));
Gr=NaN(size(z));
A=NaN(size(z));
P=NaN(size(z));
Z=z;
X(X==0)=NaN;
Y(Y==0)=NaN;
Z(Y==0)=NaN;
warning('off')
[Xh,Yh]=find(isnan(Z)==0);
% Calculate the relative slope height and gradient and the probability for
% cells with data
for j=1:length(Xh)
    % identify cells within a window of the GridSize (Grid-size can be
    % chosen seperately
    [c11,c22]=find(X>X(Xh(j),Yh(j))-GridSize & X<X(Xh(j),Yh(j))+GridSize & Y>Y(Xh(j),Yh(j))-GridSize & Y<Y(Xh(j),Yh(j))+GridSize);
    [IND]=find(X>X(Xh(j),Yh(j))-GridSize & X<X(Xh(j),Yh(j))+GridSize & Y>Y(Xh(j),Yh(j))-GridSize & Y<Y(Xh(j),Yh(j))+GridSize);
    F=zeros(size(Z));
    F(c11,c22)=1;
    % Calculate difference between center cell and adjacent cells within
    % the window associated to the grid size
    Dif=Z(Xh(j),Yh(j))-Z(IND);
    Dif(Dif<0)=0;
    % calculate distance between center cell and adjacent cells
    dxy1=sqrt(abs(X(IND)-X(Xh(j),Yh(j))).^2+abs(Y(IND)-Y(Xh(j),Yh(j))).^2);
    S5=atand(Dif./dxy1); %calculate the slope gradient
    Freql=(double(Dif)/11).^2.5.*(9.5./cotd(double(S5))).^5.*(FSavg./Lsm);
    Freqb=(double(Dif)/11).^5.*(9.5./cotd(double(S5))).^5.*(FSavg./Lsm);
    Freq2=0.5.*Freql+0.5.*Freqb; %calculate the frequency of flow slides
    Probl=1-exp(-Freql);
    Probb=1-exp(-Freqb);
    Prob2=1-exp(-Freq2); % transfer frequence to probability
    if nansum(Prob2)==0 % if probability is zero, everything is zero
        A(Xh(j),Yh(j))=0;
        P(Xh(j),Yh(j))=0;
        Gr(Xh(j),Yh(j))=0;
        Hr(Xh(j),Yh(j))=0;
    else
        AA=find(Prob2==nanmax(Prob2)); % select cell with highest probability
        cx=c11(AA(1));
        cy=c22(AA(1));
        A(Xh(j),Yh(j))=90+atan2d(Y(cx,cy)-Y(Xh(j),Yh(j)),X(cx,cy)-X(Xh(j),Yh(j))); %determine the aspect and transpose for correct direction in Matlab
        P(Xh(j),Yh(j))=Prob2(AA(1)); % give the probability
        Gr(Xh(j),Yh(j))=S5(AA(1)); % give the associated slope gradient
        Hr(Xh(j),Yh(j))=Dif(AA(1)); % give the associated slope height
    end
end

%% Determine shoal margin of the DEM
M=interp2(XX,YY,zM,X,Y); % interpolate median elevation on X,Y model
ROI= shaperead('D:\Vaklodingen\ROI_Prob.shp'); % Load edge of the Western Scheldt
[ROIp,ROIn]=inpolygon(X,Y,ROI.X,ROI.Y); % interpolate Region of Interest on X,Y model
% Determine the shoals and the shoal margins
Platen=Z-M; % Platen is the Dutch translation for Shoal, didn't change this before uploading the code
Platen(Platen<0)=0;
Platen(Platen>0)=1;
Platen=double(Platen.*ROIp);
Plaat=zeros(size(Platen));
Plaat(find(Platen==1))=1;

se = strel('sphere',1);
Pd=imdilate(Plaat,se); % Increase shoal size
Pe=imerode(Plaat,se); % Recuded shoal size
Randen=Pd-Pe; % Difference is the shoal margin
% remove the edges of the model, to limit the change of instable model at
% the edge
Randen(1:end,1:10)=0;
Randen(1:end,end-9:end)=0;
Randen(1:10,1:end)=0;
Randen(end-9:end,1:end)=0;

[B1,L1,N1,A1]= bwboundaries(Randen);
L1(L1>N1)=0; % Map that identifies seperate shoals in the WS

for jj=1:N1
    if length(B1{jj})<100
        L1(L1==jj)=0;
    end
end
[B1,L1,N1,A1]= bwboundaries(L1);
L1(L1>N1)=0; % Map that identifies seperate shoals in the WS

L3=zeros(size(L1));
L3(L1>0)=1;
L3=L3.*dZ;
L4=zeros(size(L1));
L4(L3<-0.25)=1;
se = strel('sphere',2);
L5=imdilate(L4,se); % Increase shoal size
L1(L5==1)=0;

%% Go through the regions of interest and add the shoal margin collapse
ZZ2=Z; % make changes on a new Z-matrix
for ip=1:N1 % go to all tidal flats
    while 1
        % select the first shoal margin and determine locations with a
        % probability of at least 10^-6
        Shoal=zeros(size(L1));
        Shoal(L1==ip)=1;
        ROI1=P.*Shoal;
        
        ROI1(ROI1<10.^-4)=0;
        ROI1(isnan(ROI1)==1)=0;
        if nansum(ROI1(:))==0
            break
        end
        [Xa,Ya]=find(ROI1==nanmax(ROI1(:))); %ONLY THE HIGHEST VALUE
        Shoal=zeros(size(L1));
        Shoal(Xa(1),Ya(1))=1;
        %% Determine the shape of the shoal margin collapse
        random=rand(1);
        PlaatA=xa(min(find(random<yac)));
        PlaatV=xv(min(find(random<yav)));
        Major=(sqrt(PlaatA))./(sqrt(pi()).*sqrt(sqrt(1-Ec.^2))); %calculte major axis
        Minor=sqrt(Major.^2-Ec.^2*Major.^2); %calculta minor axis
        F=sqrt((Major./2).^2-(Minor./2).^2); %focu length
        Al=sqrt(F.^2+(Minor./2).^2); %length to the focu
        ec=F./Al; %check eccentricity (NOT necessary)
        PlaatB=(Major).*(sqrt(Major.^2-Ec.^2*Major.^2)).*pi(); %check area of collapse
        Thickness=(PlaatV.*3)./Major./Minor./(4/3)./pi(); %determine thickness of the collapse
        [x,y,z]=ellipsoid(X(Xa(1),Ya(1)),Y(Xa(1),Ya(1)),0,Major,Minor,Thickness); %make an ellipsoid at point with the highest probability for a shoal margin collapse
        %remove values on one side of the ellipse
        x(y<(Y(Xa(1),Ya(1))-Minor/2))=[];
        z(y<(Y(Xa(1),Ya(1))-Minor/2))=[];
        y(y<(Y(Xa(1),Ya(1))-Minor/2))=[];
        % remove values that are negative, as we only have 1/4 of an
        % ellipsoid
        x(z<0)=[];
        y(z<0)=[];
        z(z<0)=[];
        % get X and Y from the center of the location of the collapse
        X1=X(Xa(1),Ya(1))-Major:1:X(Xa(1),Ya(1))+Major;
        Y1=Y(Xa(1),Ya(1))-Major:1:Y(Xa(1),Ya(1))+Major;
        [XP,YP]=meshgrid(X1,Y1); % create new grid with equal lengths
        ZP=griddata(x,y,z,XP,YP); % determine the Z of the shoal margin collapse
        % Get one side reduced by the slope
        curve=max(ZP');
        cs=curve(min(find(curve>0)))./(min(find(curve>0))-floor(length(curve)./2));
        cx=1:(floor(length(curve)./2)-min(find(curve>0)));
        cy=[cs.*cx+curve(min(find(curve>0)))];
        ZPc=zeros(size(ZP));
        for doe=1:length(ZPc);
            ZPc(min(find(curve>0)):floor(length(curve)./2)-1,doe)=cy;
        end
        ZP2=ZP-ZPc;
        ZP2(isnan(ZP)==1)=NaN;
        ZP2(ZP2<0)=0;
        % Get collapse into DEM, by also rotating to slope direction
        ZPr=imrotate(ZP2,A(Xa(1),Ya(1))+180); % rotate image towards the actual slope direction (aspect, A)
        [AB,BA]=find(ZPr==nanmax(ZPr(:)));%(round(ZP(ceil(length(ZP)./2+1),ceil(length(ZP)./2+1)).*10000)./10000)); % find center
        %         ZPr2=ZPr(AB(1)-length(ZP2)./2+1:AB(1)+length(ZP2)./2,BA(1)-length(ZP2)./2:BA(1)+length(ZP2)./2-1); % resize image
        ZPrB=ZPr(length(ZPr)./2-length(XP)./2:length(ZPr)./2+length(XP)./2,length(ZPr)./2-length(XP)./2:length(ZPr)./2+length(XP)./2); % resize image
        ZPr2=ZPrB(1:length(XP),1:length(XP));
        ZPr2(isnan(ZPr2)==1)=0;
        ZZp=interp2(XP,flipud(YP),ZPr2,X,Y); %interpolate on actual grid the collapse
        ZZp(isnan(ZZp)==1)=0;
        
        ZZ2=ZZ2-ZZp; % Create new DEM from the original minus the collapse
        %% Determine the location of the shoal margin collapse deposit
        PlaatA2=PlaatA.*1.5;
        Major=(sqrt(PlaatA2))./(sqrt(pi()).*sqrt(sqrt(1-Ec2.^2))); %calculte major axis
        Minor=sqrt(Major.^2-Ec2.^2*Major.^2); %calculta minor axis
        F=sqrt((Major./2).^2-(Minor./2).^2); %focu length
        Al=sqrt(F.^2+(Minor./2).^2); %length to the focu
        ec=F./Al; %check eccentricity
        PlaatB2=(Major).*(sqrt(Major.^2-Ec2.^2*Major.^2)).*pi(); %check area of collapse
        Thickness=(PlaatV.*2)./Major./Minor./(4/3)./pi(); %determine thickness of the collapse
        % depending on slope direction determine which cells are affected,
        % where collapsed sediment should be transported
        if A(Xa(1),Ya(1))<0 && A(Xa(1),Ya(1))>=-90
            xt=X(Xa(1),Ya(1))+sind(A(Xa(1),Ya(1))).*GridSize2;
            yt=Y(Xa(1),Ya(1))-cosd(A(Xa(1),Ya(1))).*GridSize2;
        elseif A(Xa(1),Ya(1))<-90 && A(Xa(1),Ya(1))>=-180
            xt=X(Xa(1),Ya(1))+sind(A(Xa(1),Ya(1))).*GridSize2;
            yt=Y(Xa(1),Ya(1))-cosd(A(Xa(1),Ya(1))).*GridSize2;
        elseif A(Xa(1),Ya(1))<-180 && A(Xa(1),Ya(1))>=-270
            xt=X(Xa(1),Ya(1))+sind(A(Xa(1),Ya(1))).*GridSize2;
            yt=Y(Xa(1),Ya(1))-cosd(A(Xa(1),Ya(1))).*GridSize2;
        else
            xt=X(Xa(1),Ya(1))+sind(A(Xa(1),Ya(1))).*GridSize2;
            yt=Y(Xa(1),Ya(1))-cosd(A(Xa(1),Ya(1))).*GridSize2;
        end
        X = vs_get(NFS,'map-const','XCOR','quiet');
        Y = vs_get(NFS,'map-const','YCOR','quiet');
        
        zt=griddata(X,Y,Z,xt,yt);
        xc=xt(min(find(diff(zt)>-0.5)));
        yc=yt(min(find(diff(zt)>-0.5)));
        if isempty(xc)
            xc=xt(min(find(diff(zt)>-1.5)));
            yc=yt(min(find(diff(zt)>-1.5)));
        end
        [x,y,z]=ellipsoid(xc,yc,0,Major,Minor,Thickness); %make an ellipsoid at point with the highest probability for a shoal margin collapse
        x(z<0)=[];
        y(z<0)=[];
        z(z<0)=[];
        X1=xc-Major:1:xc+Major;
        Y1=yc-Major:1:yc+Major;
        [XP,YP]=meshgrid(X1,Y1); % create new grid with equal lengths
        ZP=griddata(x,y,z,XP,YP); % determine the Z of the shoal margin collapse
        X(X==0)=NaN;
        Y(Y==0)=NaN;
        ZZd=interp2(XP,YP,ZP,X,Y);
        ZZd(isnan(ZZd)==1)=0;
        
        ZZ2=ZZ2+ZZd; %create new DEM from the collapsed shoal and add the deposited sediment
        break
    end
end

warning('on')
