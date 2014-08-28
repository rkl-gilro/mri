function [fima,B,mm]=biascorrector3Dv2(ima,umbral,factor,res)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Usage: [fima,B]=BiasCorrector3D(ima,umbral,factor,res)
%
%   Input:
%     ima: volume to be filtered
%     umbral: Background elimination (0:Dont erase (already done by the user), 1: Erase)
%     factor: Subresolution factor (default value=4)
%     res: voxel resolution (e.g. [1,1,3])
%
%   Output:
%     fima: filtered volume
%        B: estimated bias field
%
%  Authors: Jose Vicente Manjon Herrera
%           Marco Andre Ferreira Rodrigues
%    Email: jmanjon@fis.upv.es
%     Date: 6-02-2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Initiating bias correction ...');

ima=double(ima);
imao=ima;
so=size(ima);

% erase background
if(umbral>0)
    
    if(umbral==1)
        
        nima=log(ima+1);
        nima=nima/max(nima(:));
        ind=find(nima);
        i = graythresh(nima(ind));
        mm=nima>i;
        clear nima;
        
        if(length(so)==2) so(3)=1;end;
        for i=1:so(3)
            mm(:,:,i)=imerode(mm(:,:,i),ones(3,3));
            mm(:,:,i)=imdilate(mm(:,:,i),ones(7,7));
            mm(:,:,i)=imfill(mm(:,:,i),'holes');
            mm(:,:,i)=imerode(mm(:,:,i),ones(3,3));
        end
        
        ima=ima.*mm;
        
    end
    
else  mm=ima>0; end;

rima=ima;
s=size(ima);
if(length(s)==2) s(3)=1;end;

disp('Subsampling ...');

% filtering

if(s(3)==1)
    fh=fspecial('average');
    ima=imfilter(ima,fh);
else
    ima=cvolfilter(ima,1);
end

%subsampling
% r=factor;
% nr=round(s./((s.*res)/r));
% ima=volresize(ima,s./nr,'spline');
s=size(ima);
if(length(s)==2) s(3)=1;end;

m1 = max(max(max(ima(:))));
size(m1)
ima=(ima*256)./m1;

media=mean(ima(:));

disp('processing ...');

lv(1)=centropiaconjunta(ima)
ii=2;
ind=1;
gg=7;%[3,5,7,9,11];

% start process
for k=1%1:5,
    
    g=gg(k);
    gx=g;
    t1=(s(1))/gx;
    gy=round((s(2))/t1);
    gz=round((s(3))/t1);
    if(gz<3) gz=3;end;
    if(s(3)==1) gz=1;end;
    
    matrix=[gx,gy,gz]
    if(k==1)
        inc=0.1;
        coe=ones(gx,gy,gz);
    else
        inc=0.01;
        coe=volresize(coe,[gx,gy,gz],'spline');
    end
    t=size(coe);
    if(length(t)==2) t(3)=1;end;
    
    % initial estimate
    
    B=volresize(coe,s,'spline');
    
    oldB=B;
    F=ima./B;
    F=(F*media)/mean(F(:));
    
    if(ind==1) v=centropiaconjunta(F);
    else v=lv(ii);
    end;
    
    while(1)
        
        for c1=1:t(1)
            for c2=1:t(2)
                for c3=1:t(3)
                    ss=-1;
                    while(ss<2)
                        oldcoe=coe;
                        coe(c1,c2,c3)=coe(c1,c2,c3)+inc*ss;
                        if(coe(c1,c2,c3)<0.001)
                            coe=oldcoe;
                            ss=ss+2;
                            continue;
                        end;
                        
                        % create spline surface
                        
                        i=(c3-1)*t(2)*t(1)+(c2-1)*t(1)+c1;
                        
                        if(t(3)>1)
                            d=zeros(t);
                            [x,y,z] = ndgrid(1:(gx-1)/(s(1)-1):gx,1:(gy-1)/(s(2)-1):gy,1:(gz-1)/(s(3)-1):gz);
                            d(i)=1;
                            coeff  = spm_bsplinc(d,[3 3 3 0 0 0]);
                            base=spm_bsplins(coeff,x,y,z,[3 3 3 0 0 0]);
                        else
                            d=zeros(t);
                            [x,y,z] = ndgrid(1:(gx-1)/(s(1)-1):gx,1:(gy-1)/(s(2)-1):gy,1:1);
                            d(i)=1;
                            coeff  = spm_bsplinc(d,[3 3 1 0 0 0]);
                            base=spm_bsplins(coeff,x,y,z,[3 3 1 0 0 0]);
                        end
                        
                        B = B + base*(coe(i)-oldcoe(i));
                        
                        clear base;
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        if(min(B(:))<=0)
                            coe=oldcoe;
                            B=oldB;
                            ss=ss+2;
                            continue;
                        end
                        
                        % correct data
                        F=ima./B;
                        F=(F*media)/mean(F(:));
                        
                        % cost function
                        va=centropiaconjunta(F);
                        
                        % update
                        if(va<v)
                            v=va;
                            oldB=B;
                            ss=ss+3;
                        else
                            coe=oldcoe;
                            B=oldB;
                        end
                        ss=ss+2;
                    end
                    
                end
            end
        end
        
        lv(ii)=v;
        disp(v);
        
        if(ii>1)
            if((lv(ii-1)-lv(ii))<0.001)
                inc=inc/2;
            end;
        end
        if(inc<0.001) break; end;
        
        ii=ii+1;
        
    end
    
    lvd(ind)=lv(ii-1);
    ind=ind+1;
    
    if(ind>2)
        if((lvd(ind-2)-lvd(ind-1))<0.001) break;end;
    end
    
end

disp('finishing ...');

% final solution
clear ima;
clear B;
clear F;
%clear bases;

s=size(rima);
B=volresize(coe,s,'spline');
ind=find(mm==0);
B(ind)=1;
media=mean(rima(:));
fima=imao./B;
fima=(fima*media)/mean(fima(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%            Auxuliary functions
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol]=volresize(data,g,method)

s=size(data);
if(length(s)==2)
    [x,y,z] = ndgrid(1:(s(1)-1)/(g(1)-1):s(1),1:(s(2)-1)/(g(2)-1):s(2),1:1);
    coeff  = spm_bsplinc(data,[3 3 1 0 0 0]);
    vol=spm_bsplins(coeff,x,y,z,[3 3 1 0 0 0]);
else
    [x,y,z] = ndgrid(1:(s(1)-1)/(g(1)-1):s(1),1:(s(2)-1)/(g(2)-1):s(2),1:(s(3)-1)/(g(3)-1):s(3));
    coeff   = spm_bsplinc(data,[3 3 3 0 0 0]);
    vol     = spm_bsplins(coeff,x,y,z,[3 3 3 0 0 0]);
end


