%% function to combine Z T P RH and Modes of 4STAR modes 
% the altitude(Z), temperature(T), pressure(P), relative humidity(RH), 
% the operation mode (Str), the mode (Md) from the different modes 
% of 4star star.mat; s is yyyymmddstar.mat file, e.g.
% 'F:\ARISE\starmat_ARISE\20140904star.mat'
% utc_range=[19.3,22.9]; [0 28] is default
function so=combine_star_params(s,utc_range)

if nargin < 2; 
    utc_range=[0,28];
end
uu=fieldnames(s);
so.Tst=[0.];
so.Pst=[0.];
so.Z=[0.];
so.RH=[0.];
so.t=[0.];
so.Str=[0.];
so.Md=[0.];
for u=1:length(uu);
nn=findstr(uu{u},'vis');
    if nn; 
        kk=size(s.(uu{u}));
        if kk(2)==1;
          tt=s.(uu{u}).Tst;  
          so.Tst=[so.Tst,tt'];
          pp=s.(uu{u}).Pst;  
          so.Pst=[so.Pst,pp'];
          rr=s.(uu{u}).RH;   
          so.RH=[so.RH,rr'];
          zz=s.(uu{u}).Alt;  
          so.Z=[so.Z,zz'];
          ti=s.(uu{u}).t;    
          so.t=[so.t,ti'];
          st=s.(uu{u}).Str;    
          so.Str=[so.Str,st'];
          md=s.(uu{u}).Md;    
          so.Md=[so.Md,md'];
        else
          for i=1:kk(2);
            tt=s.(uu{u})(i).Tst;  
            so.Tst=[so.Tst,tt'];
            pp=s.(uu{u})(i).Pst;  
            so.Pst=[so.Pst,pp'];
            rr=s.(uu{u})(i).RH;   
            so.RH=[so.RH,rr'];
            zz=s.(uu{u})(i).Alt;  
            so.Z=[so.Z,zz'];
            ti=s.(uu{u})(i).t;    
            so.t=[so.t,ti'];
            st=s.(uu{u})(i).Str;    
            so.Str=[so.Str,st'];
            md=s.(uu{u})(i).Md;    
            so.Md=[so.Md,md'];
          end;
        end;
    end;
end;
[so.t,is]=sort(so.t(2:end));
so.utc=t2utch(so.t);
if nargin >= 2; 
    it=find(so.utc >= utc_range(1) & so.utc <=utc_range(2));
    is=is(it);
    so.utc=so.utc(it);
    so.t=so.t(it);
end;
so.Tst=so.Tst(is)';
so.Pst=so.Pst(is)';
so.Z=so.Z(is)';
so.RH=so.RH(is)';
so.Str=so.Str(is)';
so.Md=so.Md(is)';
return;
