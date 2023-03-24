function [G,uf,vf]=grad2D_higher(X,Y,u,v,X0,Y0,l)
%G=struct(u_x,u_y,v_x,v_y);
uf=zeros(length(X0),1);
vf=zeros(length(Y0),1);
uf_x=uf;uf_y=uf;vf_x=uf;vf_y=uf;
for i=1:length(X0(:))
    xs=X0(i)-l-min(X(:))+1:X0(i)+l-min(X(:))+1;
    ys=Y0(i)-l-min(Y(:))+1:Y0(i)+l-min(Y(:))+1;
%     [xs',ys']
    Xs=X(ys,xs);Xs=Xs(:);
    Ys=Y(ys,xs);Ys=Ys(:);
    us=u(ys,xs);us=us(:);
    vs=v(ys,xs);vs=vs(:);
    cu = fit2dPolySVD(Xs,Ys,us, 2 );
    cv = fit2dPolySVD(Xs,Ys,vs, 2 );
    uf(i)=eval2dPoly( X0(i), Y0(i), cu );
    vf(i)=eval2dPoly( X0(i), Y0(i), cv );
    uf_x(i)=cu(4)+cu(5)*Y0(i)+2*cu(6)*X0(i);
    uf_y(i)=cu(2)+2*cu(3)*Y0(i)+cu(5)*X0(i);
    vf_x(i)=cv(4)+cv(5)*Y0(i)+2*cv(6)*X0(i);
    vf_y(i)=cv(2)+2*cv(3)*Y0(i)+cv(5)*X0(i);
end
uf=reshape(uf,size(X0));
vf=reshape(vf,size(X0));
uf_x=reshape(uf_x,size(X0));
uf_y=reshape(uf_y,size(X0));
vf_x=reshape(vf_x,size(X0));
vf_y=reshape(vf_y,size(X0));
G.u_x=uf_x;
G.u_y=uf_y;
G.v_x=vf_x;
G.v_y=vf_y;


