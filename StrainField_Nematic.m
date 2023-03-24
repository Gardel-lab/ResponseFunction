%% Inputs
%Optical flow field files
Ind_files=1;  %if optical flow files for each frame is stored in different files. If stacked, change to zero
%file name format
OF_fmt='/Users/mehdi/Downloads/POackes/Analyzed/OpticalFlow/cell01/dt1/flowField_%03d.mat';

%stack file name if optical flow is saved in one mat file;
% OF_name='/Volumes/BackUp_Chicago_Nov20/Daniel/T31_p_c/FlowField_sep/flowField_name.mat';

%lambda file
filename_lambda='/Users/mehdi/Downloads/POackes/Analyzed/OpticalFlow/cell01/dt1/lambda_cell1.mat';
%outputname
Response_name='/Users/mehdi/Downloads/POackes/Analyzed/OpticalFlow/cell01/dt1/FF_name.mat';

%lagtime: the difference between measurement frame of pertubation field and observation
%filed 
tau=[0,1];

%frames:   frames that pertubation fields are measured. 
frames=1:2:6;

%size of the box to measure response function
lmax= 500; %set to a value smaller then the optical flow field

%% running the code
load(filename_lambda);
xytuv=MT(:,[1,2,9,5,6]);

if Ind_files==1
    load(sprintf(OF_fmt,frames(1)+tau(1)));
    v=flowField;
    Sz=size(v);
    xs=Sz(2); ys=Sz(1);
else
    load(OF_name);
    v=flowFiled(frames(1)).v;
    Sz=size(v);
    xs=Sz(2); ys=Sz(1);
end

[Xi,Yi]=meshgrid(1:xs,1:ys);
[Xmesh,Ymesh]=meshgrid(1:2*lmax+1,1:2*lmax+1);
Response=struct('tau',[],'ux',[],'uy',[]);
for j=1:length(tau)
    f = waitbar(0,['Running correlation function  ',num2str(tau(j))]);
    FF=zeros(length(Xmesh(:)),11);
    for i=1:length(frames)-tau(j);
        if Ind_files==1
            load(sprintf(OF_fmt,frames(i)+tau(j)));
            v=flowField;
        else
            kk=frames(i)+tau(j);
            v=flowField(kk).u;
            v(:,:,2)=flowField(kk).v;
        end
        idx=find(xytuv(:,3)==frames(i));
        xyuv=xytuv(idx,[1,2,4,5]);
        FF=corrFlow_fast(Xi,Yi,v,xyuv,FF,lmax);
        waitbar(i/length(frames),f)
    end
    uxs2=FF(:,10)./FF(:,1);
    uys2=FF(:,9)./FF(:,1);
    uxs2=real(reshape(uxs2,2*lmax+1,2*lmax+1));
    uys2=real(reshape(uys2,2*lmax+1,2*lmax+1));
    uxs=uxs2/2-fliplr(uxs2)/2;
    uys=uys2/2+fliplr(uys2)/2;
    ux=uxs/2+flipud(uxs)/2;
    uy=uys/2-flipud(uys)/2;
    Response(j).tau=tau(j);
    Response(j).ux=ux;
    Response(j).uy=uy;
    
    save(Response_name,'Xmesh','Ymesh','Response','frames');
    close(f)
    
end