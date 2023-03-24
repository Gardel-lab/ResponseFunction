function FF=corrFlow_fast(Xi,Yi,v,xyuv,FF,lmax)

ta=atan2(xyuv(:,4),xyuv(:,3));
ta(ta<0)=2*pi+ta(ta<0);

for i=1:length(xyuv(:,1))
   dxa=xyuv(i,3);
   dya=xyuv(i,4);
   dxb=v(:,:,1);
   dyb=v(:,:,2);
   va=(dxa.^2+dya.^2).^.5;
   vb=(dxb.^2+dyb.^2).^.5;
   rx=Xi-xyuv(i,1);
   ry=Yi-xyuv(i,2);
   r=(rx.^2+ry.^2).^.5;
   tr=atan2(ry,rx);
   tr(tr<0)=2*pi+tr(tr<0);
   theta=ta(i)-tr;
%    theta(theta<0)=2*pi+theta(theta<0);
%    theta(theta>pi)=theta(theta>pi)-2*pi;
   udv=(dxa*dxb+dya*dyb)./va./vb;
   sg=(rx.*dxb+ry.*dyb)-(rx.*dxa+ry.*dya).*(dxb.*dxa+dyb.*dya)./va.^2;
   signt=sign(sin(theta));signt(signt==0)=1;
   upv=sin(acos(udv)).*sign(sg).*signt;
   ddx1=udv;
   ddy1=upv;
   ddx2=udv.*vb;
   ddy2=upv.*vb;
   ddx3=udv.*vb./va;
   ddy3=upv.*vb./va;
   ddx4=udv.*vb.*va;
   ddy4=upv.*vb.*va;
   va=repmat(va,size(vb));
   C=[va(:),vb(:),ddx1(:),ddy1(:),ddx2(:),ddy2(:),ddx3(:),ddy3(:), ddx4(:),ddy4(:)];
   rxn=round(r.*sin(theta))+lmax+1;
   ryn=round(r.*cos(theta))+lmax+1;
   idx=find(rxn>0 & ryn>0 & rxn<(2*lmax+1) & ryn<(2*lmax+1) & ~isnan(vb));
   idx=unique(idx);
   rxn=rxn(idx);
   ryn=ryn(idx);
   C=C(idx,:);
   id=(rxn-1)*(2*lmax+1)+ryn;
   FF(id,1:end-1)=FF(id,1:end-1)+C;
   FF(id,end)=FF(id,end)+1;
end

