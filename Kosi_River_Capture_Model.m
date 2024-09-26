
clc
clear
close all

addpath('data')
addpath('data/kml_shapefile')
addpath('data/topotoolbox-2.4')

% read in DEM

load KosiDEM.mat

%genearate flow direction and upstream area
FD = FLOWobj(DEM,'preprocess','carve');
DEM = imposemin(FD,DEM,0.0001);
A = flowacc(FD);

%% read in kml file
fault = kml_shapefile('data/MFT.kml');
 
lon = fault.X;
lat = fault.Y;
 
%extract UTM zone from lat and lon
p1 = [lat,lon];
z1 = utmzone(p1);
 
%get geoid and construct projection structure
[ellipsoid,estr] = utmgeoid(z1);
utmstruct = defaultm('utm');
utmstruct.zone = z1;
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct);

%transform coordinates to easting and northing
[fault.X, fault.Y] = mfwdtran(utmstruct,lat,lon);

fault_line = line2GRIDobj(DEM, fault);
 
%make the fault line twice as wide - prevents streams that flow diagonally
%being missed
fault_line = dilate(fault_line,ones(2));
 
% imagesc(fault_line);

%define caputure event location
capture_lat=28.049;
capture_lon=87.379;
[capture.X,capture.Y]=mfwdtran(utmstruct,capture_lat,capture_lon);



%define a point for cut the trunk and its tributaries
cut_lat2=26.936;
cut_lon2=87.155;
[cut2.X,cut2.Y]=mfwdtran(utmstruct,cut_lat2,cut_lon2);


%% extract stream network for model
S = STREAMobj(FD,'minarea',1e6,'unit','map');
S = modify(S, 'upstreamto', fault_line);
S = klargestconncomps(S,1);

% remove part tributaries to reduce the amount of computation
S1 = modify(S,'tributaryto',trunk(S));
S_3tri = klargestconncomps(S1,3);
S_tri = STREAMobj2cell(S_3tri);
S_Left = S_tri{1};
S_Right = S_tri{2};

S_R1 = modify(S_Right,'tributaryto',trunk(S_Right));
S_R2 = modify(S_R1,'tributaryto',trunk(S_R1));

S_Mid = modify(S,'rmnodes',S_Left);
S_Mid = modify(S_Mid,'rmnodes',S_Right);
S_M1 = modify(S_Mid,'tributaryto',trunk(S_Mid));
S_M2 = modify(S_M1,'tributaryto',trunk(S_M1));

[temp.X,temp.Y]=mfwdtran(utmstruct,27.59,85.32);
[xn,yn] = snap2stream(S,temp.X,temp.Y);
[IX] = coord2ind(DEM,xn,yn);
S_Left_t = modify(S_Left,'downstreamto',IX);

S_L1 = modify(S_Left,'tributaryto',S_Left_t);
S_L2 = modify(S_L1,'tributaryto',trunk(S_L1));
S_L3 = modify(S_L2,'tributaryto',trunk(S_L2));

S_L4 = klargestconncomps(S_L1,6);
S_L5 = modify(S_Left,'rmnodes',S_L4);
S_L5 = modify(S_L5,'rmnodes',S_Left_t);


%remove tributaries of tributaries
S3 = modify(S,'rmnodes',S_R2);
S3 = modify(S3,'rmnodes',S_M2);
S3 = modify(S3,'rmnodes',S_L3);
S3 = modify(S3,'rmnodes',S_L5);
S3 = removeshortstreams(S3,10000);
S = S3;

%plot the stream network
figure
imagesc(DEM);
hold on;
plot(S,'r');

clear S1 S2 S3 S22 S23 S24 S25 S_L1 S_L2 S_L3 S_L4 S_L5 S_Left S_M1 S_M2 S_Mid S_R1 S_R2 S_Right S_3tri S_tri S_Left_t;

%% prepare for model


% upstream area.
a = getnal(S,A)*DEM.cellsize^2;

% elevation of river nodes
z=DEM.Z(S.IXgrid);

%distance to fault vector
dist_fault = zeros(length(S.ix), 1);
fault_points = find(fault_line.Z == 1); % or "find(fault_line);"
[x_fault, y_fault] = getcoordinates(fault_line, 'matrix');
x = S.x;
y = S.y;
 
for i = 1:length(S.x)
    
    for j = 1:length(fault_points)
        ij = fault_points(j);
      
        distance_to_fault(j) = sqrt((((x(i) - x_fault(ij))^2)+((y(i)-y_fault(ij))^2)));
        
    end
    min_dist = min(distance_to_fault);
    dist_fault(i) = min_dist; %distance to fault vector
end

% find capture location on trunk
min_dist = (1e9);
for i = 1:length(S.x)
    
        distance_to_capture = sqrt((((x(i) - capture.X)^2)+((y(i)-capture.Y))^2));
        
    if distance_to_capture < min_dist
        min_dist = distance_to_capture;
        capture_loc = i;
    end
end
 

% build the river network downstream of the capture point for the steady state model

loc_ix=find(S.ix==capture_loc);
upstream_capture=modify(S,'upstreamto',S.IXgrid(S.ixc(loc_ix)));
S_modify = modify(S,'rmnodes',upstream_capture);

% list trunk points downstream the capture point
list_points=zeros(size(S.ix));
n_points=0;
ij=capture_loc;
for i=1:length(S.ix)
    list_points(i)=ij;
    loc_ix=find(S.ix==ij);
    n_points=n_points+1;
     
     if isempty(loc_ix)==1
         break        
     end
    
    ij=S.ixc(loc_ix);
    
end
list_points(list_points==0)=[];

%calculate the upstream area before capture
area_orig=a;
for i=1:n_points
    area_orig(list_points(i))=area_orig(list_points(i))-a(capture_loc);
end
area_orig(area_orig<0)=0.1;

 for i = numel(S.ix):-1:1
        if i==capture_loc
            area_orig(i)=-9999;
        end
 end

for i = numel(S.ix):-1:1
        if area_orig(S.ixc(i))==-9999
           area_orig(S.ix(i))=-9999;  
        end
end

area_orig(area_orig<0)=0.1;

% build a array 'record' to mark river nodes downstream(1) and upstream(0) the knickpoint
record_modify=zeros(length(area_orig),1);
for i=1:length(area_orig)
    
    if area_orig(i) == 0.1
    record_modify(i)=0;
    else
        record_modify(i)=1;
    end
end

% delete upstream point to fit the modifed S
area_orig_dc=area_orig.*record_modify;
area_orig_dc(area_orig_dc==0)=[];

a_modify=a.*record_modify;
a_modify(a_modify==0)=[];



% build an array for calculating misfit
record_misfit=ones(length(S_modify.x),1);

% Select the river nodes for the two tributaries on the left
[cut.X,cut.Y]=mfwdtran(utmstruct,27.252,86.194);
[cut_x,cut_y,cut_IX]=snap2stream(S_modify,cut.X,cut.Y);
loc_cut_ix=find(S_modify.IXgrid==cut_IX);

list_tri=zeros(size(S_modify.ix));
ij=loc_cut_ix;
loc_temp=[];
for i=1:length(S_modify.ix)
    list_tri(i)=ij;
    loc_ix=find(S_modify.ixc==ij);
    
    if isempty(loc_ix)==1
        if isempty(loc_temp)==0
          loc_ix=loc_temp(1);
          loc_temp(1)=[];
        end
    end
    
    if length(loc_ix)==1
        ij=S_modify.ix(loc_ix);
    else
        loc_temp=[loc_temp;loc_ix(2)];
        loc_ix=loc_ix(1);
        ij=S_modify.ix(loc_ix);
    end
    
    loc_ix2=find(S_modify.ixc==ij);
    if isempty(loc_ix2)==1
        if isempty(loc_temp)==1
            break
        end
    end

end


list_tri(list_tri==0)=[];


for i=1:length(list_tri)
    record_misfit(list_tri(i))=0;
       
end

% Select the river nodes with elevation >4000 m
for i=1:length(record_misfit)
    if DEM.Z(S_modify.IXgrid(i))>4000
        record_misfit(i)=0;
    end
end

% Select the river nodes with upstream area < 5e6
for i=1:length(record_misfit)
    if a_modify(i)<5e6
       record_misfit(i)=0;
    end
end


%% genearate initial river profile
%discretise time steps
dt=1000;

%assign values to parameters
capture_time=89550;
k = 1.159704386891183e-05;
m = 0.42637643;
n_fluv= 0.96752578;

tot_time=capture_time;
t=[0:dt:tot_time];

%initialise elevation values
z_mod = zeros(length(S.x), length(t));
k_mat=z_mod;
z_mod = z_mod+min(z);

%calculate distance between each node
dx = distance(S, 'node_to_node');


u_mat=k_mat;
for j=1:length(t)
k_mat(:,j)=k; 
end

%define uplift rates for every nodes
u_mat(:)=0;
for i = 1:length(S.x)
     u_mat(i,:)=0.5*1000/1e6;
     if dist_fault(i)<75000        
         u_mat(i,:)=1*1000/1e6;
     elseif dist_fault(i)>75000 && dist_fault(i)<120000
        u_mat(i,:)= 4*1000/1e6;
     end
end

figure
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false,'ticklabels','nice');
 plotc(S,u_mat(:,2))
 colorbar
title('uplift rate map')
 

% steady state model

z_modify = DEM.Z(S_modify.IXgrid);
z_ss=z_modify;
z_ss(:)=0;
z_ss = z_ss+min(z_modify);

u_mat2=u_mat(:,1).*record_modify;
u_mat2(u_mat2==0)=[];

dx_modify = distance(S_modify, 'node_to_node');

    for i = numel(S_modify.ix):-1:1
        z_ss(S_modify.ix(i)) =  z_ss(S_modify.ixc(i)) + ((u_mat2(S_modify.ix(i))/(k*(area_orig_dc(S_modify.ix(i))^m)))^(1/n_fluv))*dx_modify(S_modify.ix(i));       
    end



figure
plot(S_modify.distance,z_ss, '.')
hold on 
plotdz(S_modify,DEM)
hold off




%% river evolution after capture event

z_mod_modify=zeros(length(S_modify.x), length(t));

for i=1:length(t)
    z_mod_modify(:,i)=z_ss;
end


nk=10

figure
for j = 2:length(t)-1
    for i = numel(S_modify.ix):-1:1
        if n_fluv==1
        c = k_mat(S_modify.ix(i),j)*(a_modify(S_modify.ix(i))^m)*dt/dx_modify(S_modify.ix(i));       
        z_mod_modify(S_modify.ix(i), j) =  (z_mod_modify(S_modify.ix(i), j-1) + z_mod_modify(S_modify.ixc(i), j) * c)/(1+c);       
        else
              c=dt*k_mat(S_modify.ix(i),j)*(a_modify(S_modify.ix(i))^m)/(dx_modify(S_modify.ix(i))^n_fluv);
               hnn= z_mod_modify(S_modify.ix(i), j-1);
                 for ik=1:nk
                   hn=max(hnn,z_mod_modify(S_modify.ixc(i), j)+1.e-10);
                   hnn=hn-(hn-z_mod_modify(S_modify.ix(i), j-1)+c*(hn-z_mod_modify(S_modify.ixc(i), j))^n_fluv)/(1.d0+n_fluv*c*(hn-(z_mod_modify(S_modify.ixc(i), j)))^(n_fluv-1.));
                 end
               dh=hn-z_mod_modify(S_modify.ix(i), j-1);
               z_mod_modify(S_modify.ix(i), j)= z_mod_modify(S_modify.ix(i), j-1)+dh;
         end
    end
    
    z_mod_modify(:, j) = z_mod_modify(:, j)+ (u_mat2(:)*dt);
    z_mod_modify(1,j)=0; 
    
    if mod(j,10)==0
         
    plot(S_modify.distance,z_mod_modify(:, j), '.')

    %hold on
    %time = (1-j/(length(t)-1))*(tot_time/1e3);
    %txt = [num2str(time) 'ka'];
    %text(1e5,10000,txt)
    %hold off  
    pause(0.0001)
    
    end
   
    
end
hold on 

plotdz(S_modify,DEM);
hold off

%%
% calculate misfit
diff_elev=z_mod_modify(:,length(t)-1)-z_modify(:);
diff_elev=diff_elev.*record_misfit;

misfit_model=0;

for i=1:length(diff_elev)
    misfit_model = misfit_model+abs(diff_elev(i));
end

misfit_model = misfit_model/length(diff_elev)



% make a figure show trunk(Arun River) model profile vs. real profile 
[temp.X,temp.Y]=mfwdtran(utmstruct,27.25,86.18);
[xn,yn] = snap2stream(S,temp.X,temp.Y);
[IX] = coord2ind(DEM,xn,yn);
S_modelpart = modify(S_modify,'upstreamto',IX);
S_modelpart = modify(S_modify,'rmnodes',S_modelpart);


% Rivers to the east and west of the Kosi
[temp.X,temp.Y]=mfwdtran(utmstruct,26.95,87.16);
[xn,yn] = snap2stream(S,temp.X,temp.Y);
[IX] = coord2ind(DEM,xn,yn);

S_temp1 = modify(S_modelpart,'upstreamto',IX);
S_temp1 = modify(S_modelpart,'rmnodes',S_temp1);

S_trunk = modify(S_modelpart,'rmnodes',S_temp1);


DEM_mod=DEM;
for i=1:length(z_mod_modify(:, length(t)-1))
    DEM_mod.Z(S_modify.IXgrid(i))=z_mod_modify(i,length(t)-1);
    
end

S_trunk_modify = modify(S_trunk,'rmnodes',upstream_capture);

figure 
hold on
plotdz(S_trunk_modify,DEM_mod)
plotdz(S_trunk_modify,DEM)
axis([0 2.5e5 0 8000]) 
hold off


