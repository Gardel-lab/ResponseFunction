%% Inputs
%Optical flow field files
Ind_files=1;  %if optical flow files for each frame is stored in different files. If stacked, change to zero
%file name format
OF_fmt='/Users/mehdi/Downloads/POackes/Analyzed/OpticalFlow/cell01/dt1/flowField_%03d.mat';

%stack file name if optical flow is saved in one mat file;
% OF_name='/Volumes/BackUp_Chicago_Nov20/Daniel/T31_p_c/FlowField_sep/flowField_name.mat';

%lambd file
filename_lambda='/Users/mehdi/Downloads/POackes/Analyzed/OpticalFlow/cell01/dt1/lambda_cell1.mat';

%frames
frames=1:2:6;

% spacing for lambda calculation
de=100;     %ditance from edge
ds=20;      %grid spacing
meshsize=8;  % mesh size in pixel to measure the displacement gradient


%% running the code

if Ind_files==1
    load(sprintf(OF_fmt,frames(1)));
    u=flowField(:,:,1);
    ly=length(u(:,1));
    lx=length(u(1,:));
else
    load(OF_name);
    u=flowFiled(frames(1)).vu;
    xs=Sz(1); ys=Sz(2);
    ly=length(u(:,1));
    lx=length(u(1,:));
end


MT2=[];

[Xm,Ym]=meshgrid((de:ds:lx-de-1),(ds:ds:ly-ds-1));
[X,Y]=meshgrid((1:1:lx),(1:1:ly));
xr=round(Xm(:));
yr=round(Ym(:));
f = waitbar(0,'calculating lambdas');
for i=1:length(frames);
    
    if Ind_files==1
        load(sprintf(OF_fmt,frames(i)));
        u=flowField(:,:,1);
        v=flowField(:,:,2);
    else
        kk=frames(i);
        u=flowField(kk).u;
        v=flowField(kk).v;
    end  
    [G,uf,vf]=grad2D_higher(X,Y,u,v,xr,yr,meshsize);
    a=G.u_x-G.v_y;
    b=G.u_y+G.v_x;
    l=sqrt(a.^2+b.^2);
    cx1=-b./(a-l);cy1=1;n1=sqrt(cx1.^2+cy1.^2);
    cx1=cx1./n1;cy1=cy1./n1;
    cx2=-b./(a+l);cy2=1;n2=sqrt(cx2.^2+cy2.^2);
    cx2=cx2./n2;cy2=cy2./n2;
    mt2=[xr yr uf vf cx1.*l cy1.*l cx2.*l cy2.*l xr*0+frames(i)];
    MT2=[MT2;mt2];
    waitbar(i/length(frames),f)
end
MT=MT2;
save(filename_lambda,'MT')
