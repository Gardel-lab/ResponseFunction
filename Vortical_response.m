%% Inputs
%Optical flow field files
Ind_files=0;  %if optical flow files for each frame is stored in different files. If stacked, change to zero
%file name format
OF_fmt='Data/FlowField_%03d.mat';

%stack file name if optical flow is saved in one mat file;
% OF_name=''/Users/mehdi/Dropbox/GardelLab/Manuscripts/Method/GitHubPackage/Data/Nematic_flowField.mat';

%outputname
Response_name='/Users/mehdi/Dropbox/GardelLab/Manuscripts/Method/GitHubPackage/Data/VorticalResponse.csv';

%lagtime
tau=0;% 2 3 6 10 20 30];

%frames
frames=1:2;

%maximum radial distance to calculate the response
maxR=500;

%grids spacing to calulate reponse function
dsi=10;


% spacing for lambda calculation
de=110;     %ditance from edge
ds=20;      %grid spacing
meshsize=8;  % mesh size in pixel to measure the displacement gradient



[Xmesh,Ymesh]=meshgrid(-maxR:ds:maxR,-maxR:ds:maxR);
Xmesh=Xmesh(:);
Ymesh=Ymesh(:);



if Ind_files==1
    load(sprintf(OF_fmt,frames(1)+tau(1)));
    u=flowField(:,:,1);
    v=flowField(:,:,2);
else
    load(OF_name);
    u=flowField(frames(1)).v;
    v=flowField(frames(1)).u;
end



ly=length(u(:,1)); lx=length(u(1,:));

[Xm,Ym]=meshgrid((de:ds:lx-de-1),(ds:ds:ly-ds-1));
Xm=Xm(:);
Ym=Ym(:);
[X,Y]=meshgrid((1:1:lx),(1:1:ly));


[Xi,Yi]=meshgrid((1:dsi:lx),(1:dsi:ly));
Xi=Xi(:); Yi=Yi(:);

[amat,bmat]=meshgrid(1:length(Xm),1:length(Xi(:)));
amat=amat(:);
bmat=bmat(:);
Ex=Xi(bmat)-Xm(amat);
Ey=Yi(bmat)-Ym(amat);
r=(Ex.^2+Ey.^2).^.5;
ex=Ex./r;
ey=Ey./r;

Ur=zeros(length(amat),length(tau));
Ut=Ur; Oc=Ur;


R2=round(r/2);
Ri=unique(R2);
Ri=Ri(Ri<maxR);
idx=unique(round(logspace(log10(1),log10(length(Ri)),150)));
Ri=Ri(idx);
Utc=zeros(length(Ri),length(tau));
Urc=Utc;
f = waitbar(0,['Running correlation function  ',num2str(frames(1))]);
header{1}='R';
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
    
    
    [G,uf,vf]=grad2D_higher(X,Y,u,v,Xm,Ym,meshsize);
    oc=G.u_y-G.v_x;
    for j=1:length(tau);
        if i+tau(j)<=max(frames)
            
            if Ind_files==1
                load(sprintf(OF_fmt,frames(i)+tau(j)));
                u=flowField(:,:,1);
                v=flowField(:,:,2);
            else
                kk=frames(i)+tau(j);
                u=flowField(kk).u;
                v=flowField(kk).v;
            end            
            u=u(1:dsi:ly,1:dsi:lx);
            v=v(1:dsi:ly,1:dsi:lx);
            u=u(:);
            v=v(:);
            urdot=ex.*u(bmat)+ey.*v(bmat);
            utdot=ey.*u(bmat)-ex.*v(bmat);
            Ur(:,j)=Ur(:,j)+urdot.*oc(amat);
            Ut(:,j)=Ut(:,j)+utdot.*oc(amat);
            Oc(:,j)=Oc(:,j)+abs(oc(amat));
        end
    end
    waitbar(i/length(frames),f)
    
end

tic
for i=1:length(Ri)
    %     idx=find(R2==Ri(i));
    %     Utc(i,:)=sum(Ut(idx,:))./sum(Oc(idx,:));
    Utc(i,:)=sum(Ut(R2==Ri(i),:))./sum(Oc(R2==Ri(i),:));
    %     Urc(i,:)=sum(Ur(idx,:))./sum(Oc(idx,:));
end
for i=1:length(tau)
    header{i+1}=sprintf('tau_%d',tau(i));
end

T=array2table([Ri Utc],'VariableNames',header);
writetable(T,Response_name);
close(f);
