clear all;
clc;

no_wires=2;
max_harmonic=1;
no_harmonic=2*max_harmonic+1;
matrix_rank=no_wires*no_harmonic*4*2;

x1=-2e-3;
y1=0;
x2=2e-3;
y2=0;
dx=x1-x2;
dy=y1-y2;
gammapq=[atan2(x1-x2,y1-y2), atan2(x2-x1,y2-y1)];
sigmapq=sqrt(dx*dx+dy*dy);
radius_out=1e-3;
radius_in=0.8e-3;
freq=1e6;

muo=4e-7*pi;
eo=8.8541878176e-12;
co=1/sqrt(eo*muo);
ko=2*pi*freq/co;
conductivity=1e3;
omega=2*pi*freq;
omega_muo=omega*muo;
omega_eo=omega*eo;

erd=2.3-i*2.3*2.5e-4;
erw=1-i*conductivity/(2*pi*1e6*eo);
omega_ed=omega_eo*erd;
omega_ew=omega_eo*erw;
kod=ko*sqrt(erd);
kow=ko*sqrt(erw);

% beta=0.99*ko;
% kt=sqrt(ko*ko-beta*beta);
% kt2=ko*ko-beta*beta;
% ktd=sqrt(kod*kod-beta*beta);
% ktd2=kod*kod-beta*beta;
% ktw=sqrt(kow*kow-beta*beta);
% ktw2=kow*kow-beta*beta;
k=1;
for beta=0.9*ko:(1.09*ko-0.9*ko)/100:1.1*ko
    kt=sqrt(ko*ko-beta*beta);
    kt2=ko*ko-beta*beta;
    ktd=sqrt(kod*kod-beta*beta);
    ktd2=kod*kod-beta*beta;
    ktw=sqrt(kow*kow-beta*beta);
    ktw2=kow*kow-beta*beta;
for ii=1:no_wires;
    
    for ja=-max_harmonic:max_harmonic
    
    % bessel terms in  atmosphere------------------------------------------------------------------------
        hankad=kt*(besselh(ja-1,2,kt*radius_out)-ja/(kt*radius_out)*besselh(ja,2,kt*radius_out));
        hanka=besselh(ja,2,kt*radius_out);
        besa=besselj(ja,kt*radius_out);
        besad=kt*(besselj(ja-1,kt*radius_out)-ja/(kt*radius_out)*besselj(ja,kt*radius_out));
    % bessel terms in dielectric------------------------------------------------------------------     
        besj1=besselj(ja,ktd*radius_in);
        besj1d=ktd*(besselj(ja-1,ktd*radius_in)-ja/(ktd*radius_in)*besselj(ja,ktd*radius_in));
        besy1=bessely(ja,ktd*radius_in);
        besy1d=ktd*(bessely(ja-1,ktd*radius_in)-ja/(ktd*radius_in)*bessely(ja,ktd*radius_in));
        besj2=besselj(ja,ktd*radius_out);
        besj2d=ktd*(besselj(ja-1,ktd*radius_out)-ja/(ktd*radius_out)*besselj(ja,ktd*radius_out));
        besy2=bessely(ja,ktd*radius_out);
        besy2d=ktd*(bessely(ja-1,ktd*radius_out)-ja/(ktd*radius_out)*bessely(ja,ktd*radius_out));
   
    % bessel terms in conductor------------------------------------------------------------------------------------------------- 
        besb=besselj(ja,ktw*radius_in);
        besbd=ktw*(besselj(ja-1,ktw*radius_in)-ja/(ktw*radius_in)*besselj(ja,ktw*radius_in));
    
    % coefficient derivation
        
%{          
        d/dz B<n>(z) = 0.5* (B<n-1> - B<n+1>)
        B<n+1> = 2n/z B<n> - B<n-1>
        hence,
        d/dz B<n>(z) = 0.5* (B<n-1> - 2n/z B<n> + B<n-1>) = B<n-1> - n/z B<n>
        d/dr B(ktr) = kt d/dz B(z)
%}
    
    % matrix element calculatiohn for self terms------------------------------------------------
        matrix_index=(ii-1)*no_harmonic+ja+max_harmonic;    %%ii
        
        cmatrix(8*matrix_index+1,8*matrix_index+1)=kt2;
        cmatrix(8*matrix_index+1,8*matrix_index+3)=-besj1*kt2;
        cmatrix(8*matrix_index+1,8*matrix_index+4)=-besy1*kt2;
        
        cmatrix(8*matrix_index+2,8*matrix_index+2)=kt2;
        cmatrix(8*matrix_index+2,8*matrix_index+5)=-besj1*kt2;
        cmatrix(8*matrix_index+2,8*matrix_index+6)=-besy1*kt2;
        
        cmatrix(8*matrix_index+3,8*matrix_index+1)=i*(-beta/radius_in)*(-i*ja)*kt2/ktw2;
        cmatrix(8*matrix_index+3,8*matrix_index+2)=i*omega_muo*kt2/ktw2*besbd/besb;
        cmatrix(8*matrix_index+3,8*matrix_index+3)=-i*(kt2/ktd2)*(-beta/radius_in)*(-i*ja)*besj1;
        cmatrix(8*matrix_index+3,8*matrix_index+4)=-i*(kt2/ktd2)*(-beta/radius_in)*(-i*ja)*besy1;
        cmatrix(8*matrix_index+3,8*matrix_index+5)=-i*(kt2/ktd2)*omega_muo*besj1d;
        cmatrix(8*matrix_index+3,8*matrix_index+6)=-i*(kt2/ktd2)*omega_muo*besy1d;
        
        cmatrix(8*matrix_index+4,8*matrix_index+1)=i*(-omega_ew)*kt2/ktw2*besbd/besb;
        cmatrix(8*matrix_index+4,8*matrix_index+2)=i*(-beta/radius_in)*(-i*ja)*kt2/ktw2;
        cmatrix(8*matrix_index+4,8*matrix_index+3)=-i*(-omega_ed)*kt2/ktd2*besj1d;
        cmatrix(8*matrix_index+4,8*matrix_index+4)=-i*(-omega_ed)*kt2/ktd2*besy1d;
        cmatrix(8*matrix_index+4,8*matrix_index+5)=-i*(kt2/ktd2)*(-beta/radius_in)*(-i*ja)*besj1;
        cmatrix(8*matrix_index+4,8*matrix_index+6)=-i*(kt2/ktd2)*(-beta/radius_in)*(-i*ja)*besy1;
        
        cmatrix(8*matrix_index+5,8*matrix_index+3)=besj2*kt2;
        cmatrix(8*matrix_index+5,8*matrix_index+4)=besy2*kt2;
        cmatrix(8*matrix_index+5,8*matrix_index+7)=-kt2;
        
        cmatrix(8*matrix_index+6,8*matrix_index+5)=besj2*kt2;
        cmatrix(8*matrix_index+6,8*matrix_index+6)=besy2*kt2;
        cmatrix(8*matrix_index+6,8*matrix_index+8)=-kt2;
        
        cmatrix(8*matrix_index+7,8*matrix_index+3)=i*(kt2/ktd2)*(-beta/radius_out)*(-i*ja)*besj2;
        cmatrix(8*matrix_index+7,8*matrix_index+4)=i*(kt2/ktd2)*(-beta/radius_out)*(-i*ja)*besy2;
        cmatrix(8*matrix_index+7,8*matrix_index+5)=i*(kt2/ktd2)*omega_muo*besj2d;
        cmatrix(8*matrix_index+7,8*matrix_index+6)=i*(kt2/ktd2)*omega_muo*besy2d;
        cmatrix(8*matrix_index+7,8*matrix_index+7)=-i*(-beta/radius_out)*(-i*ja);
        cmatrix(8*matrix_index+7,8*matrix_index+8)=-i*omega_muo*hankad/hanka;
        
        cmatrix(8*matrix_index+8,8*matrix_index+3)=i*(-omega_ed)*kt2/ktd2*besj2d;
        cmatrix(8*matrix_index+8,8*matrix_index+4)=i*(-omega_ed)*kt2/ktd2*besy2d;
        cmatrix(8*matrix_index+8,8*matrix_index+5)=i*(kt2/ktd2)*(-beta/radius_out)*(-i*ja)*besj2;
        cmatrix(8*matrix_index+8,8*matrix_index+6)=i*(kt2/ktd2)*(-beta/radius_out)*(-i*ja)*besy2;
        cmatrix(8*matrix_index+8,8*matrix_index+7)=-i*(-omega_eo)*hankad/hanka;
        cmatrix(8*matrix_index+8,8*matrix_index+8)=-i*(-beta/radius_out)*(-i*ja);
        
    %matrix elements for coupling terms---------------------------------
        
        for ib=1:no_wires   
           
            if ib~=ii
            
               for jb=-max_harmonic:max_harmonic
                    hanknm=besselh(jb-ja,2,kt*sigmapq);
                    hankm=besselh(jb,2,kt*radius_out);
                    besjn=besselj(ja,kt*radius_out);
            
                    Tnm=exp(-i*(jb-ja)*gammapq(ii))*hanknm*besjn/hankm;
                      
                    iia=(ii-1)*no_harmonic+ja+max_harmonic;
                    iib=(ib-1)*no_harmonic+jb+max_harmonic;
                   
                    cmatrix(8*iia+5,8*iib+7)=-1*Tnm*kt2;
                    cmatrix(8*iia+6,4*iib+8)=-1*Tnm*kt2;
                  
                    cmatrix(8*iia+7,8*iib+7)=-1*Tnm*i*(-beta/radius_out)*(-i*ja);
                    cmatrix(8*iia+7,8*iib+8)=-1*Tnm*i*omega_muo*besad/besa;
                    
                    cmatrix(8*iia+8,8*iib+7)=-1*Tnm*i*(-omega_eo)*besad/besa;
                    cmatrix(8*iia+8,8*iib+8)=-1*Tnm*i*(-beta/radius_out)*(-i*ja);
               end
            
            end 
        end  
            
    end
    
end



[U,S,V]=svd(cmatrix);
r=rank(cmatrix)
k=k+1;
% sv(k)=S(35,35);
end
