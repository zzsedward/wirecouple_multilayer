#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <time.h>
#include <stdarg.h>
#include <functional>
#include <vector>

#include "cable.h"

//g++ -o cabler cable.cpp lapack_LINUX.a blas_LINUX.a NAGSX.lib -lm -lgfortran
int main(int argc, char* argv[]){
    cout<<"Hello World!"<<endl;

    const int max_harmonic(1); // ie will use from -this to +this

    int number_layers(2);

    double* radii= new double[number_layers];

    radii[0]=1.e-3;
    radii[1]=1.2e-3;
    //radii[2]=1.4e-3;
//
    const double freq_MHz(1);

    const double ko(2.*constants::get_pi()*freq_MHz*1e6/constants::get_co());
        cout<<"ko: "<<2.*constants::get_pi()*freq_MHz*1e6/constants::get_co()<<endl;
       // cout<<"co: "<<constants::get_co()<<endl;
    const double kappa (1e9); //conductivity of conductor

    const complex epsilon_r(1.-complex(0.,1.)* kappa/(2.*constants::get_pi()*freq_MHz*1e6*constants::get_eo()));
        cout<<"epsilon_r: "<<epsilon_r<<endl;

    const complex kow(ko*sqrt(epsilon_r));
        cout<<"kow: "<<kow<<endl;

    const double er_plastic(1);
    const double loss_tan(0);
    const complex epsilon_rd(er_plastic, -er_plastic*loss_tan);
      const complex epsilon_rd2(1,0);
        cout<<"epsilon dielectric: "<<epsilon_rd<<endl;

    const complex kod(ko*sqrt(epsilon_rd));
        cout<<"kod: "<<kod<<endl;

    complex* relative_epsilon=new complex[number_layers];

    relative_epsilon[0]=epsilon_r;
    relative_epsilon[1]=epsilon_rd;
    //relative_epsilon[2]=epsilon_rd2;

    vector<wire> wires;

    wires.push_back(wire(-3e-3,0,radii,relative_epsilon,number_layers));

    wires.push_back(wire(3e-3,0,radii,relative_epsilon,number_layers));

    //wires.push_back(wire(1e-3,3e-3,radii,relative_epsilon,number_layers));

     const double min_beta(0.95*ko);
        cout<<"min beta: "<<min_beta<<endl;

     const double max_beta(1.1*ko);
        cout<<"max beta: "<<max_beta<<endl;

     const int no_beta_steps(100);

     vector<complex> amps;

     double min_sv(1e50);

     complex min_sv_beta(0.,0.);

     vector<complex> best_amps;

     vector<double> sv;

     ofstream fout("sv.csv");

     int which_mode(0); // 0 means lowest sv, 1 next lowest etc

/*     for(int i=0;i<no_beta_steps;++i) {fout<<","<<(min_beta+double(i)/double(no_beta_steps-1)*(max_beta-min_beta))/ko;}

     for(int i=0;i<no_beta_steps;++i) {

            double beta_step_imag(-0.01+double(i)/double(no_beta_steps-0.01)*(0.01-(-0.01)));

            fout<<"\n"<<beta_step_imag;

        for(int ii=0;ii<no_beta_steps;++ii){

            double beta_step_real(min_beta+double(ii)/double(no_beta_steps-1)*(max_beta-min_beta));

            const complex beta(beta_step_real,beta_step_imag);

             //cout<<"beta: "<<imag(beta)<<"\r"<<flush;
            cout<<"step: "<<i*no_beta_steps+ii<<"\r"<<flush;

            int no_solutions=get_determinant(wires,ko,beta,max_harmonic,amps,which_mode,sv);

             //cout<<"no_solution: "<<no_solutions<<endl;

            if(sv[no_solutions-1-which_mode]<min_sv){

                 min_sv=sv[no_solutions-1-which_mode];

                 min_sv_beta=beta;

                 best_amps=amps;}

                 fout<<","<<sv[no_solutions-1];
        }
     }

     cout<<"\nBest="<<min_sv_beta<<" "<<min_sv;cout.flush();*/

    //------------beta variation with real part only-----------------------------------
       // double frequency_min(1e6),frequency_max(9.99e8);

     /*   for(int i=0;i<no_beta_steps;++i) {
            //const double frequency(frequency_min+double(i)/double(no_beta_steps-1)*(frequency_max-frequency_min));
            //const complex beta(2.*constants::get_pi()*frequency/constants::get_co(),0);
         const complex beta(min_beta+double(i)/double(no_beta_steps-1)*(max_beta-min_beta),0.);

         cout<<"beta: "<<real(beta)<<"\r"<<flush;

         int no_solutions=get_determinant(wires,ko,beta,max_harmonic,amps,which_mode,sv);

         if(sv[no_solutions-1-which_mode]<min_sv){

             min_sv=sv[no_solutions-1-which_mode];

             min_sv_beta=beta;

	     best_amps=amps;}

        fout<<"\n"<<real(beta)/ko<<","<<sv[no_solutions-1];
        //fout<<"\n"<<frequency<<","<<sv[no_solutions-1];

        if(no_solutions>1){fout<<","<<sv[no_solutions-2];}
        if(no_solutions>2){fout<<","<<sv[no_solutions-3];}
        if(no_solutions>3){fout<<","<<sv[no_solutions-4];}
        if(no_solutions>4){fout<<","<<sv[no_solutions-5];}
        if(no_solutions>5){fout<<","<<sv[no_solutions-6];}

        //for(int s=0;s<no_solutions;s++){fout<<","<<sv[no_solutions-s-1];}
	  }


         fout.close();
        cout<<"\nBest="<<real(min_sv_beta)/ko<<" "<<min_sv;cout.flush();*/
         min_sv_beta=complex(0.020958,1e-6);
         cout<<"\nbeta: "<<min_sv_beta;
         int no_solutions=get_determinant(wires,ko,min_sv_beta,max_harmonic,best_amps,which_mode,sv);
          //cout<<"\ncondition: "<<sv[0]/sv[no_solutions-1]<<endl;
          //cout<<"\nlargest sv: "<<sv[0]<<endl;
         //cout<<"\nsmallest sv: "<<sv[no_solutions-1]<<endl;*/
         // plot_field2(wires,ko,epsilon_rd,epsilon_r,complex(0.99*ko,0.0),max_harmonic,amps);
         complex ex,ey,ez,hx,hy,hz;
         get_fields(wires,ko,min_sv_beta,max_harmonic,best_amps,-1.80001e-3,0.,ex,ey,ez,hx,hy,hz);
         cout<<"\nEx1: "<<ex;
         get_fields(wires,ko,min_sv_beta,max_harmonic,best_amps,-1.8e-3,0.,ex,ey,ez,hx,hy,hz);
         cout<<"\nEx2: "<<ex;

//         plot_field(wires,ko,min_sv_beta,max_harmonic,best_amps);
//         plot_field2(wires,ko,min_sv_beta,max_harmonic,best_amps);
}


int get_determinant(vector<wire>& wires,
                    const double ko,
                    const complex beta,
                    const int max_harmonic,
                    vector<complex>& amps,
                    int which_amps,
                    vector<double>& sv){

     const complex kt2(ko*ko-beta*beta);
     const complex kt(sqrt(kt2));
     //cout<<"\nkt: "<<kt<<endl;
//-----------------matrix initialization--------------------------------------
     const int no_wires(wires.size());     //number of wires

     const int no_harmonics(2*max_harmonic+1);    //number of harmonics

     int matrix_rank(0);   //number of boundary conditions * number of wires * number of harmonics
     for(int i=0;i<no_wires;i++){
        matrix_rank+=4*wires[i].no_layers*no_harmonics;
     }

     //cout<<"\nmatrix rank: "<<matrix_rank<<endl;

     static vector<complex> cmatrix;

     cmatrix.reserve(matrix_rank*matrix_rank); //the overall matrix

     memset(&cmatrix[0],0,sizeof(complex)*matrix_rank*matrix_rank);

     static vector<double> dmatrix;   //the overall matrix with double type

     dmatrix.reserve(4*matrix_rank*matrix_rank);

     memset(&dmatrix[0],0,sizeof(double)*4*matrix_rank*matrix_rank);

//----------------Bessel function parameters-------------------------------------

     int nz,ifail,scale_len=1;

     static complex jj(0.,1.);

     char scale='u';            //return unscaled results

     double zc[2], fnu=0., gampq;   //zc: argument value; fnu: starting order; gampq: gammapq(angle between pth and qth wire, with respect to y axis)

     int no_bes_needed(max(2,max_harmonic+1)), hank_kind(2), hank_its;

     static vector<complex> besj1,besy1,besj1d,besy1d,besj2,besy2,besj2d,besy2d,besjm1,besjm1d,besym1,besym1d,besjm2,besym2,besjm2d,besym2d;
     static vector<complex> hanka,hankad,besa,besad,besj3,besj3d,besy3,besy3d,besja,besya,besjb,besyb;
     static vector<double> temp,cwrk;

     besj1.reserve(no_bes_needed);   //bessel term in mth layer
     besj1d.reserve(no_bes_needed);

     besj2.reserve(no_bes_needed);   //bessel term in m+1th layer
     besj2d.reserve(no_bes_needed);

     besy1.reserve(no_bes_needed);  //Bessel(ktd*r)   term in mth layer
     besy1d.reserve(no_bes_needed);

     besy2.reserve(no_bes_needed);  //Bessel(ktd*r)   term in m+1th layer
     besy2d.reserve(no_bes_needed);

     besj3.reserve(no_bes_needed);   //besselj term in last layer
     besj3d.reserve(no_bes_needed);

     besy3.reserve(no_bes_needed);  //Bessely term in last layer
     besy3d.reserve(no_bes_needed);

     besjm1.reserve(no_bes_needed);   //besselj term in last layer
     besjm1d.reserve(no_bes_needed);

     besym1.reserve(no_bes_needed);  //Bessely term in last layer
     besym1d.reserve(no_bes_needed);

     besjm2.reserve(no_bes_needed);   //besselj term in last layer
     besjm2d.reserve(no_bes_needed);

     besym2.reserve(no_bes_needed);  //Bessely term in last layer
     besym2d.reserve(no_bes_needed);

     hanka.reserve(no_bes_needed);
     hankad.reserve(no_bes_needed);

     besa.reserve(no_bes_needed);
     besad.reserve(no_bes_needed);

     besja.reserve(no_bes_needed);   //besselj term in last layer

     besya.reserve(no_bes_needed);  //Bessely term in last layer

     besjb.reserve(no_bes_needed);   //besselj term in last layer

     besyb.reserve(no_bes_needed);  //Bessely term in last layer

     temp.reserve(2*no_bes_needed); //temporary store of bessel and hankel terms
     cwrk.reserve(2*no_bes_needed); // work space for hankel

//--------------------------omega_mu and omega_epsilon in the air--------------
     const double w_ua(ko*constants::get_zo());
     const double w_ea(ko*constants::get_yo());

//--------------------------matrix element calculation------------------------
    int block_index=0;  //boundary block

    for(int ip=0;ip<no_wires;++ip){
        //cout<<"\nip: "<<ip;
        const int number_layers(wires[ip].no_layers);

            const double radius_m0(wires[ip].radius[0]);
            //cout<<"\nfirst layer radius: "<<radius_m0<<endl;

            const complex epsilonr_m0(wires[ip].epsilonr[0]);

            const complex kom0(ko*sqrt(epsilonr_m0));

            const complex kt2m0(kom0*kom0-beta*beta);

            const complex ktm0(sqrt(kt2m0));

            const complex w_em0(w_ea*epsilonr_m0);

            zc[0]=real(ktm0*radius_m0);
            zc[1]=imag(ktm0*radius_m0);
            //cout<<"\nko0: "<<kom0<<endl;
            S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

                for(int k=0;k<no_bes_needed;++k){besj1[k]=complex(temp[2*k],temp[2*k+1]);}

                besj1d[0]=-ktm0*besj1[1];

                for(int k=1;k<no_bes_needed;++k){besj1d[k]=ktm0*(besj1[k-1]-double(k)/complex(zc[0],zc[1])*besj1[k]);}

            const complex epsilonr_m1(wires[ip].epsilonr[1]);
                //cout<<"\nepsilonr for second layer: "<<epsilonr_m1<<endl;

            const complex kom1(ko*sqrt(epsilonr_m1));

            const complex kt2m1(kom1*kom1-beta*beta);

            const complex ktm1(sqrt(kt2m1));

            const complex w_em1(w_ea*epsilonr_m1);

            zc[0]=real(ktm1*radius_m0);
            zc[1]=imag(ktm1*radius_m0);

            if(zc[0]==0&&zc[1]==0){cout<<"\nktm1: "<<ktm1<<" beta: "<<beta;}

            S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

                for(int k=0;k<no_bes_needed;++k){besj2[k]=complex(temp[2*k],temp[2*k+1]);}//cout<<"\nbesj2 ["<<k<<"]: "<<besj2[k]<<endl;}

                besj2d[0]=-ktm1*besj2[1];

                for(int k=1;k<no_bes_needed;++k){besj2d[k]=ktm1*(besj2[k-1]-double(k)/complex(zc[0],zc[1])*besj2[k]);}

            S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);
                for(int k=0;k<no_bes_needed;++k){besy2[k]=complex(temp[2*k],temp[2*k+1]);}

                besy2d[0]=-ktm1*besy2[1];

                for(int k=1;k<no_bes_needed;++k){besy2d[k]=ktm1*(besy2[k-1]-double(k)/complex(zc[0],zc[1])*besy2[k]);}

    //-------------------last layer field------------------------------------
            const double radius_lastlayer(wires[ip].radius[number_layers-1]);
                //cout<<"\nradius of last layer: "<<radius_lastlayer<<endl;

            const complex epsilonr_lastlayer(wires[ip].epsilonr[number_layers-1]);
                //cout<<"\nepsilonr for last layer: "<<epsilonr_lastlayer<<endl;

            const complex kom_last(ko*sqrt(epsilonr_lastlayer));

            const complex kt2mlast(kom_last*kom_last-beta*beta);

            const complex ktmlast(sqrt(kt2mlast));

            const complex w_emlast(w_ea*epsilonr_lastlayer);

            zc[0]=real(ktmlast*radius_lastlayer);
            zc[1]=imag(ktmlast*radius_lastlayer);

            if(zc[0]==0&&zc[1]==0){cout<<"\nktmlast: "<<ktmlast<<" beta: "<<beta;}

            S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

                for(int k=0;k<no_bes_needed;++k){besj3[k]=complex(temp[2*k],temp[2*k+1]);}

                besj3d[0]=-ktmlast*besj3[1];

                for(int k=1;k<no_bes_needed;++k){besj3d[k]=ktmlast*(besj3[k-1]-double(k)/complex(zc[0],zc[1])*besj3[k]);}

            S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);
                for(int k=0;k<no_bes_needed;++k){besy3[k]=complex(temp[2*k],temp[2*k+1]);}

                besy3d[0]=-ktmlast*besy3[1];

                for(int k=1;k<no_bes_needed;++k){besy3d[k]=ktmlast*(besy3[k-1]-double(k)/complex(zc[0],zc[1])*besy3[k]);}

            //----normalisation------------
            const double radius_norm(wires[ip].radius[number_layers-2]);

            zc[0]=real(ktmlast*radius_norm);
            zc[1]=imag(ktmlast*radius_norm);
            if(zc[0]==0&&zc[1]==0){cout<<"\nktmlast: "<<ktmlast<<" beta: "<<beta;}
                S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

                    for(int k=0;k<no_bes_needed;++k){besjb[k]=complex(temp[2*k],temp[2*k+1]);}

                S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);

                    for(int k=0;k<no_bes_needed;++k){besyb[k]=complex(temp[2*k],temp[2*k+1]);}

//----------------------outside field----------------------------------

            zc[0]=real(kt*radius_lastlayer);
            zc[1]=imag(kt*radius_lastlayer);

            S17DLF(&hank_kind,&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&hank_its,&ifail);

            for(int k=0;k<no_bes_needed;k++) {
                hanka[k]=complex(temp[2*k],temp[2*k+1]);}

            hankad[0]=-kt*hanka[1];
                                            //derivative of hankel for n=0;
            for(int k=1;k<no_bes_needed;k++) {
                hankad[k]=kt*(hanka[k-1]-double(k)/complex(zc[0],zc[1])*hanka[k]);}

//-----------------------------matrix fill-----------------------------------

            for(int ih=-max_harmonic;ih<=max_harmonic;++ih){

//-------------boundary conditions for first layer--------------
                const int ii(block_index+(ih+max_harmonic)*number_layers*4);

            //------------ez from ez inside and ez outside--------------------
                cmatrix[ii+0+(ii+0)*matrix_rank]=-1.*kt2;
                cmatrix[ii+0+(ii+2)*matrix_rank]=1.*kt2;
                cmatrix[ii+0+(ii+3)*matrix_rank]=1.*kt2;

                //cout<<"\nrow: "<<4*ii+0<<"  column: "<<(4*ii+3)<<": "<<besy2[abs(ih)]<<endl;
            //------------hz from hz inside and hz outside---------------------
                cmatrix[ii+1+(ii+1)*matrix_rank]=-1.*kt2;
                cmatrix[ii+1+(ii+4)*matrix_rank]=1.*kt2;
                cmatrix[ii+1+(ii+5)*matrix_rank]=1.*kt2;

            //------------ephi from ez inside and ez outside-------------------
                cmatrix[ii+2+(ii+0)*matrix_rank]=-jj*kt2/kt2m0*(-beta/radius_m0*(-jj*double(ih)));
                cmatrix[ii+2+(ii+2)*matrix_rank]=jj*kt2/kt2m1*(-beta/radius_m0*(-jj*double(ih)));
                cmatrix[ii+2+(ii+3)*matrix_rank]=jj*kt2/kt2m1*(-beta/radius_m0*(-jj*double(ih)));

            //-----------ephi from hz inside and hz outside--------------------
                cmatrix[ii+2+(ii+1)*matrix_rank]=-jj*kt2/kt2m0*w_ua*(besj1d[abs(ih)]/besj1[abs(ih)]);
                cmatrix[ii+2+(ii+4)*matrix_rank]=jj*kt2/kt2m1*w_ua*besj2d[abs(ih)]/besj2[abs(ih)];
                cmatrix[ii+2+(ii+5)*matrix_rank]=jj*kt2/kt2m1*w_ua*besy2d[abs(ih)]/besy2[abs(ih)];

            //-----------hphi from ez inside and ez outside-----------------------
                cmatrix[ii+3+(ii+0)*matrix_rank]=-jj*kt2/kt2m0*(-w_em0)*besj1d[abs(ih)]/besj1[abs(ih)];
                cmatrix[ii+3+(ii+2)*matrix_rank]=jj*kt2/kt2m1*(-w_em1)*besj2d[abs(ih)]/besj2[abs(ih)];
                cmatrix[ii+3+(ii+3)*matrix_rank]=jj*kt2/kt2m1*(-w_em1)*besy2d[abs(ih)]/besy2[abs(ih)];

            //----------hphi from hz inside and hz outside------------------------
                cmatrix[ii+3+(ii+1)*matrix_rank]=-jj*kt2/kt2m0*(-beta/radius_m0*(-jj*double(ih)));
                cmatrix[ii+3+(ii+4)*matrix_rank]=jj*kt2/kt2m1*(-beta/radius_m0*(-jj*double(ih)));
                cmatrix[ii+3+(ii+5)*matrix_rank]=jj*kt2/kt2m1*(-beta/radius_m0*(-jj*double(ih)));


    //----------boundary conditions for middle layers--------------------------
                for(int im=1;im<number_layers-1;im++){
                        //cout<<"\ndid I go in?"<<endl;

                //------------inside field----------------------------------------
                    const double radius_m(wires[ip].radius[im]);

                    const double radius_mminus1(wires[ip].radius[im-1]);

                    const complex epsilon_m(wires[ip].epsilonr[im]);

                    const complex ko_m(ko*sqrt(epsilon_m));

                    const complex kt2_m(ko_m*ko_m-beta*beta);

                    const complex kt_m(sqrt(kt2_m));

                    const complex w_em(w_ea*epsilon_m);

                    zc[0]=real(kt_m*radius_m);
                    zc[1]=imag(kt_m*radius_m);
                    if(zc[0]==0&&zc[1]==0){cout<<"\nktm: "<<kt_m<<" beta: "<<beta;}
                    S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

                        for(int k=0;k<no_bes_needed;++k){besjm1[k]=complex(temp[2*k],temp[2*k+1]);}

                        besjm1d[0]=-kt_m*besjm1[1];

                        for(int k=1;k<no_bes_needed;++k){besjm1d[k]=kt_m*(besjm1[k-1]-double(k)/complex(zc[0],zc[1])*besjm1[k]);}

                    S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);

                        for(int k=0;k<no_bes_needed;++k){besym1[k]=complex(temp[2*k],temp[2*k+1]);}

                        besym1d[0]=-kt_m*besym1[1];

                        for(int k=1;k<no_bes_needed;++k){besym1d[k]=kt_m*(besym1[k-1]-double(k)/complex(zc[0],zc[1])*besym1[k]);}

                 //----------------normalization-------------------------------
                    zc[0]=real(kt_m*radius_mminus1);
                    zc[1]=imag(kt_m*radius_mminus1);

                    S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

                        for(int k=0;k<no_bes_needed;++k){besja[k]=complex(temp[2*k],temp[2*k+1]);}

                    S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);

                        for(int k=0;k<no_bes_needed;++k){besya[k]=complex(temp[2*k],temp[2*k+1]);}

                //--------------outside field----------------------------------
                    const complex epsilon_mplus1(wires[ip].epsilonr[im+1]);

                    const complex ko_mplus1(ko*sqrt(epsilon_mplus1));

                    const complex kt2_mplus1(ko_mplus1*ko_mplus1-beta*beta);

                    const complex kt_mplus1(sqrt(kt2_mplus1));

                    const complex w_emplus1(w_ea*epsilon_mplus1);

                    zc[0]=real(kt_mplus1*radius_m);
                    zc[1]=imag(kt_mplus1*radius_m);
                    if(zc[0]==0&&zc[1]==0){cout<<"\nktmp1: "<<kt_mplus1<<" beta: "<<beta;}
                    S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

                        for(int k=0;k<no_bes_needed;++k){besjm2[k]=complex(temp[2*k],temp[2*k+1]);}

                        besjm2d[0]=-kt_mplus1*besjm2[1];

                        for(int k=1;k<no_bes_needed;++k){besjm2d[k]=kt_mplus1*(besjm2[k-1]-double(k)/complex(zc[0],zc[1])*besjm2[k]);}

                    S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);
                        for(int k=0;k<no_bes_needed;++k){besym2[k]=complex(temp[2*k],temp[2*k+1]);}

                        besym2d[0]=-kt_mplus1*besym2[1];

                        for(int k=1;k<no_bes_needed;++k){besym2d[k]=kt_mplus1*(besym2[k-1]-double(k)/complex(zc[0],zc[1])*besym2[k]);}

                //---------------
                    const int iia(block_index+(ih+max_harmonic)*number_layers*4+im*4);
                    const int iib(block_index+(ih+max_harmonic)*number_layers*4+(im+1)*4);

                //------------------Ez from ez------------------------------
                    cmatrix[iia+0+(iia-2)*matrix_rank]=-besjm1[abs(ih)]/besja[abs(ih)]*kt2;
                    cmatrix[iia+0+(iia-1)*matrix_rank]=-besym1[abs(ih)]/besya[abs(ih)]*kt2;
                    cmatrix[iia+0+(iib-2)*matrix_rank]=1.*kt2;
                    cmatrix[iia+0+(iib-1)*matrix_rank]=1.*kt2;

    //
    //                cout<<"\nez from ez inside and outside: "<<4*iia+0+(4*iib-1)*matrix_rank<<endl;

                //------------------Hz from hz------------------------------
                    cmatrix[iia+1+(iia+0)*matrix_rank]=-besjm1[abs(ih)]/besja[abs(ih)]*kt2;
                    cmatrix[iia+1+(iia+1)*matrix_rank]=-besym1[abs(ih)]/besya[abs(ih)]*kt2;
                    cmatrix[iia+1+(iib+0)*matrix_rank]=1.*kt2;
                    cmatrix[iia+1+(iib+1)*matrix_rank]=1.*kt2;

                //------------------Ephi from ez----------------------------
                    cmatrix[iia+2+(iia-2)*matrix_rank]=-jj*kt2/kt2_m*(-beta/radius_m*(-jj*double(ih)))*besjm1[abs(ih)]/besja[abs(ih)];
                    cmatrix[iia+2+(iia-1)*matrix_rank]=-jj*kt2/kt2_m*(-beta/radius_m*(-jj*double(ih)))*besym1[abs(ih)]/besya[abs(ih)];
                    cmatrix[iia+2+(iib-2)*matrix_rank]=jj*kt2/kt2_mplus1*(-beta/radius_m*(-jj*double(ih)));
                    cmatrix[iia+2+(iib-1)*matrix_rank]=jj*kt2/kt2_mplus1*(-beta/radius_m*(-jj*double(ih)));

                //------------------Ephi from hz----------------------------
                    cmatrix[iia+2+(iia+0)*matrix_rank]=-jj*kt2/kt2_m*w_ua*besjm1d[abs(ih)]/besja[abs(ih)];
                    cmatrix[iia+2+(iia+1)*matrix_rank]=-jj*kt2/kt2_m*w_ua*besym1d[abs(ih)]/besya[abs(ih)];
                    cmatrix[iia+2+(iib+0)*matrix_rank]=jj*kt2/kt2_mplus1*w_ua*besjm2d[abs(ih)]/besjm2[abs(ih)];
                    cmatrix[iia+2+(iib+1)*matrix_rank]=jj*kt2/kt2_mplus1*w_ua*besym2d[abs(ih)]/besym2[abs(ih)];

                //------------------Hphi from ez------------------------------
                    cmatrix[iia+3+(iia-2)*matrix_rank]=-jj*kt2/kt2_m*(-w_em)*besjm1d[abs(ih)]/besja[abs(ih)];
                    cmatrix[iia+3+(iia-1)*matrix_rank]=-jj*kt2/kt2_m*(-w_em)*besym1d[abs(ih)]/besya[abs(ih)];
                    cmatrix[iia+3+(iib-2)*matrix_rank]=jj*kt2/kt2_mplus1*(-w_emplus1)*besjm2d[abs(ih)]/besjm2[abs(ih)];
                    cmatrix[iia+3+(iib-1)*matrix_rank]=jj*kt2/kt2_mplus1*(-w_emplus1)*besym2d[abs(ih)]/besym2[abs(ih)];

                //------------------Hphi from hz-------------------------------

                    cmatrix[iia+3+(iia+0)*matrix_rank]=-jj*kt2/kt2_m*(-beta/radius_m*(-jj*double(ih)))*besjm1[abs(ih)]/besja[abs(ih)];
                    cmatrix[iia+3+(iia+1)*matrix_rank]=-jj*kt2/kt2_m*(-beta/radius_m*(-jj*double(ih)))*besym1[abs(ih)]/besya[abs(ih)];
                    cmatrix[iia+3+(iib+0)*matrix_rank]=jj*kt2/kt2_mplus1*(-beta/radius_m*(-jj*double(ih)));
                    cmatrix[iia+3+(iib+1)*matrix_rank]=jj*kt2/kt2_mplus1*(-beta/radius_m*(-jj*double(ih)));

                }  //bracket for middle layer iteration


    //---------boundary condition for the very outside layer---------------------
                const int iii(block_index+(ih+max_harmonic)*number_layers*4+(number_layers-1)*4);

            //---------------------ez from ez inside and outside---------------
                cmatrix[iii+0+(iii-2)*matrix_rank]=-besj3[abs(ih)] / besjb[abs(ih)]*kt2;
                cmatrix[iii+0+(iii-1)*matrix_rank]=-besy3[abs(ih)] / besyb[abs(ih)]*kt2;
                cmatrix[iii+0+(iii+2)*matrix_rank]=1.*kt2;

                    //cout<<"\nrow: "<<4*iii+0<<"  column: "<<(4*iii-1)<<": "<<cmatrix[4*iii+0+(4*iii-1)*matrix_rank]<<endl;

            //---------------------hz from hz inside and outside------------------
                cmatrix[iii+1+(iii+0)*matrix_rank]=-besj3[abs(ih)] / besjb[abs(ih)]*kt2;
                cmatrix[iii+1+(iii+1)*matrix_rank]=-besy3[abs(ih)] / besyb[abs(ih)]*kt2;
                cmatrix[iii+1+(iii+3)*matrix_rank]=1.*kt2;

            //---------------------ephi from ez inside and outside----------------
                cmatrix[iii+2+(iii-2)*matrix_rank]=-jj*kt2/kt2mlast*(-beta/radius_lastlayer*(-jj*double(ih)))*besj3[abs(ih)] / besjb[abs(ih)];
                cmatrix[iii+2+(iii-1)*matrix_rank]=-jj*kt2/kt2mlast*(-beta/radius_lastlayer*(-jj*double(ih)))*besy3[abs(ih)] / besyb[abs(ih)];
                cmatrix[iii+2+(iii+2)*matrix_rank]=jj*kt2/kt2*(-beta/radius_lastlayer*(-jj*double(ih)));

            //---------------------ephi from hz inside and outside-----------------
                cmatrix[iii+2+(iii+0)*matrix_rank]=-jj*kt2/kt2mlast*w_ua*(besj3d[abs(ih)]/besjb[abs(ih)]);
                cmatrix[iii+2+(iii+1)*matrix_rank]=-jj*kt2/kt2mlast*w_ua*(besy3d[abs(ih)]/besyb[abs(ih)]);
                cmatrix[iii+2+(iii+3)*matrix_rank]=jj*kt2/kt2*w_ua*(hankad[abs(ih)]/hanka[abs(ih)]);

            //---------------------hphi from ez inside and outside-----------------
                cmatrix[iii+3+(iii-2)*matrix_rank]=-jj*kt2/kt2mlast*(-w_emlast)*(besj3d[abs(ih)]/besjb[abs(ih)]);
                cmatrix[iii+3+(iii-1)*matrix_rank]=-jj*kt2/kt2mlast*(-w_emlast)*(besy3d[abs(ih)]/besyb[abs(ih)]);
                cmatrix[iii+3+(iii+2)*matrix_rank]=jj*kt2/kt2*(-w_ea)*hankad[abs(ih)]/hanka[abs(ih)];

            //---------------------hphi from hz inside and outside-----------------
                cmatrix[iii+3+(iii+0)*matrix_rank]=-jj*kt2/kt2mlast*(-beta/radius_lastlayer*(-jj*double(ih)))*besj3[abs(ih)] / besjb[abs(ih)];
                cmatrix[iii+3+(iii+1)*matrix_rank]=-jj*kt2/kt2mlast*(-beta/radius_lastlayer*(-jj*double(ih)))*besy3[abs(ih)] / besyb[abs(ih)];
                cmatrix[iii+3+(iii+3)*matrix_rank]=jj*kt2/kt2*(-beta/radius_lastlayer*(-jj*double(ih)));


            }   //bracket for ih iteration

        block_index+=number_layers*no_harmonics*4;     //iterate for wire p

    }    //bracket for wire p iteration


//-----------------------coupling terms---------------------------------------
    int boundary_index=0;

    for(int ia=0;ia<no_wires;++ia){

        const int ia_no_layers(wires[ia].no_layers);

        const double radius_ia(wires[ia].radius[ia_no_layers-1]);

        zc[0]=real(kt*radius_ia);
        zc[1]=imag(kt*radius_ia);

        S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

            for(int k=0;k<no_bes_needed;k++) {besa[k]=complex(temp[2*k],temp[2*k+1]);}   //coupling bessel term

            besad[0]=-kt*besa[1];

            for(int k=1;k<no_bes_needed;k++) besad[k]=kt*(besa[k-1]-double(k)/complex(zc[0],zc[1])*besa[k]);

        int coeff_index=0;

        for(int ib=0; ib<no_wires;++ib){

            const int ib_no_layers(wires[ib].no_layers);

            if(ia!=ib){

                wire_coupling_termer T(wires[ia], wires[ib],max_harmonic,kt);

                for(int ja=-max_harmonic;ja<=max_harmonic;++ja){

                    const int iaih(boundary_index+(ja+max_harmonic)*ia_no_layers*4+(ia_no_layers-1)*4);

                    for (int jb=-max_harmonic;jb<=max_harmonic;++jb){

                        const int ibih(coeff_index+(jb+max_harmonic)*ib_no_layers*4+(ib_no_layers-1)*4);

                        cmatrix[iaih+0+(ibih+2)*matrix_rank]=T(ja,jb)*kt2;

                        cmatrix[iaih+1+(ibih+3)*matrix_rank]=T(ja,jb)*kt2;

                        cmatrix[iaih+2+(ibih+2)*matrix_rank]=T(ja,jb)*jj*kt2/kt2*(-beta/radius_ia*(-jj*double(ja)));
                        cmatrix[iaih+2+(ibih+3)*matrix_rank]=T(ja,jb)*jj*kt2/kt2*w_ua*(besad[abs(ja)]/besa[abs(ja)]);

                        cmatrix[iaih+3+(ibih+2)*matrix_rank]=T(ja,jb)*jj*kt2/kt2*(-w_ea)*(besad[abs(ja)]/besa[abs(ja)]);
                        cmatrix[iaih+3+(ibih+3)*matrix_rank]=T(ja,jb)*jj*kt2/kt2*(-beta/radius_ia*(-jj*double(ja)));

                    }
                }
            }

            coeff_index+=no_harmonics*ib_no_layers*4;
        }

        boundary_index+=no_harmonics*ia_no_layers*4;
    }

//-----------------------output the matrix------------------------------
     ofstream mout("matrix.csv");
     mout<<"\nbeta: "<<beta<<endl;

     for(int i=0;i<matrix_rank;++i){

         mout<<"\n";
         for(int j=0;j<matrix_rank;++j){

                mout<<real(cmatrix[i+j*matrix_rank])<<"+"<<imag(cmatrix[i+j*matrix_rank])<<"i, ";

         }

     }

//------------convert to a real matrix--------------------------------

    for(int i=0;i<matrix_rank;++i){

        for(int j=0;j<matrix_rank;++j){

             dmatrix[(2*i+0)+(2*j+0)*(2*matrix_rank)]=real(cmatrix[i+j*matrix_rank]);
             dmatrix[(2*i+0)+(2*j+1)*(2*matrix_rank)]=-imag(cmatrix[i+j*matrix_rank]);
             dmatrix[(2*i+1)+(2*j+0)*(2*matrix_rank)]=imag(cmatrix[i+j*matrix_rank]);
             dmatrix[(2*i+1)+(2*j+1)*(2*matrix_rank)]=real(cmatrix[i+j*matrix_rank]);
        }
    }


   //------------------calculate singular value----------------------------------------------

   int sz(2*matrix_rank);

   int info(0);

   const char jobz('A');

   int lwork(2*(3*sz*sz+4*sz*sz+4*sz));

     static vector<int> iwork;

     iwork.reserve(8*sz);

     static vector<double> Q,PT,wk;

     sv.reserve(sz);

     Q.reserve(sz*sz);

     PT.reserve(sz*sz);

     wk.reserve(lwork);

             dgesdd_(&jobz,
                 &sz,
                 &sz,
                 &dmatrix[0],
                 &sz,
                 &sv[0],
                 &Q[0],
                 &sz,
                 &PT[0],
                 &sz,
                 &wk[0],
                 &lwork,
                 &iwork[0],
                 &info);


 // dmatrix = Q sv PT

 // For an sv~0, the solution amplitudes are in The PT
    ofstream ampsout("amps.csv");

     amps.resize(matrix_rank);

     for(int i=0;i<matrix_rank;++i) {
            amps[i]=complex(PT[sz-1-which_amps+sz*(2*i)],PT[sz-1-which_amps+sz*(2*i+1)]);
            //cout<<"\namps"<<i<<": "<<amps[i];
            ampsout<<amps[i]<<endl;
     }
     //cout<<"sz: "<<sz<<endl;

     for(int i=0;i<matrix_rank;++i){
        complex resulta(0.);

        for(int j=0; j<matrix_rank;++j){
                //cout<<"\ncmatrix "<<j<<": "<<cmatrix[2+j*matrix_rank]*amps[j]<<endl;
            resulta+=cmatrix[i+j*matrix_rank]*amps[j];
        }

        cout<<"\nresult: "<<resulta<<endl;
     }

     return(sz);
}


void get_fields(vector<wire>& wires,
                const double ko,
                const complex beta,
                const int max_harmonic,
                vector<complex>& amps,
                double x,
                double y,
                complex& ex,
                complex& ey,
                complex& ez,
                complex& hx,
                complex& hy,
                complex& hz){

    ex=ey=ez=hx=hy=hz=0.;
    complex er,ephi,hphi,hr;

    const int no_wires(wires.size());
    const int no_harmonics(2*max_harmonic+1);

    const double w_ua(ko*constants::get_zo());
    const double w_ea(ko*constants::get_yo());

    complex kt2=ko*ko-beta*beta;
    complex kt=sqrt(kt2);

//--------------------bessel function declaration----------------------------
    int nz,ifail,scale_len=1;

         static complex jj(0.,1.);

         char scale='u';

         double r,zc[2],fnu=0.;

     int no_bes_needed(std::max(2,max_harmonic+1)),hank_kind(2),hank_its;

     static vector<complex> hanka,hankad,hankb,hankbd,besj,besjd,besy,besyd,besa,besad,besb,besbd,besjm1,besym1;

     static vector<double> temp;
     static vector<double> cwrk;

         hanka.reserve(no_bes_needed);  //Hank(kt*r)
         hankad.reserve(no_bes_needed);

         hankb.reserve(no_bes_needed);  //Hank(kt*r)
         hankbd.reserve(no_bes_needed);

         besj.reserve(no_bes_needed);  //Besselktd*r)   term in dielectric
         besjd.reserve(no_bes_needed);

         besy.reserve(no_bes_needed);  //Besselktd*r)   term in dielectric
         besyd.reserve(no_bes_needed);

         besjm1.reserve(no_bes_needed);

         besym1.reserve(no_bes_needed);

         besa.reserve(no_bes_needed);  //Besselktd*r)   term in first layer
         besad.reserve(no_bes_needed);

         besb.reserve(no_bes_needed);  //Besselktd*r)   term in first layer boundary
         besbd.reserve(no_bes_needed);

     temp.reserve(2*no_bes_needed);
     cwrk.reserve(2*no_bes_needed);
//----------------------------------------------------------------------------------------------

    int wire_number(21),layer_number(21);

    complex epsilon_relative(1.,0.);

    for(int iw=0;iw<no_wires;++iw){

        const double dx(x-wires[iw].centre[0]);
        const double dy(y-wires[iw].centre[1]);

        const double dist(sqrt(dx*dx+dy*dy));

        const int number_layer(wires[iw].no_layers);

        for(int il=number_layer-1;il>=0;--il){
                //cout<<"\nlayer index: "<<il;

                if(dist<wires[iw].radius[il]) {

                        wire_number=iw;
                        layer_number=il;
                        epsilon_relative=wires[iw].epsilonr[il];
                }
        }
    }
    //cout<<"\n\nFind wire number and layer number: "<<endl;
    //cout<<"wire number: " <<wire_number+1<<"\nlayer number: "<<layer_number+1<<endl;
    //cout<<"\nrelative epsilon: "<<epsilon_relative<<endl;

//----------------------when points in the outside region(air)--------------------------------------
    if(wire_number>20) {

        //cout<<"\nThe point is in the air."<<endl;
        int index_coeff(0);

        for(int index_wire=0;index_wire<no_wires;++index_wire){

            const double ddx(x-wires[index_wire].centre[0]);
            const double ddy(y-wires[index_wire].centre[1]);

            const double OPdist(sqrt(ddx*ddx+ddy*ddy));
            const double phi(atan2(ddx,ddy));

            const int number_layer(wires[index_wire].no_layers);
            //cout<<"\nwire "<<iw<<" number of layers: "<<number_layer;
            const double radius_norm(wires[index_wire].radius[number_layer-1]);

            zc[0]=real(kt*radius_norm);
            zc[1]=imag(kt*radius_norm);

                S17DLF(&hank_kind,&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&hank_its,&ifail);

                for(int k=0;k<no_bes_needed;k++) hankb[k]=complex(temp[2*k],temp[2*k+1]);

            zc[0]=real(kt*OPdist);
            zc[1]=imag(kt*OPdist);
                    //cout<<"kt*r air: "<<zc[0]+jj*zc[1]<<endl;
                S17DLF(&hank_kind,&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&hank_its,&ifail);

                for(int k=0;k<no_bes_needed;k++) hanka[k]=complex(temp[2*k],temp[2*k+1]);

                hankad[0]=-kt*hanka[1];

                for(int k=1;k<no_bes_needed;k++) hankad[k]=kt*(hanka[k-1]-double(k)/complex(zc[0],zc[1])*hanka[k]);

            ez=ephi=hz=0.;
            er=hphi=hr=0.;

            for(int ih=-max_harmonic;ih<=max_harmonic;++ih){

                complex Xpoe(amps[index_coeff+(ih+max_harmonic)*number_layer*4+ number_layer*4-2]*exp(-jj*double(ih)*phi));
                complex Xpoh(amps[index_coeff+(ih+max_harmonic)*number_layer*4+ number_layer*4-1]*exp(-jj*double(ih)*phi));

                //cout<<"\nXpoh: "<<iw*number_layer*no_harmonics*4+(ih+max_harmonic)*number_layer*4+ number_layer*4-1<<endl;

                ez+=hanka[abs(ih)]/hankb[abs(ih)]*Xpoe*kt2;
                hz+=hanka[abs(ih)]/hankb[abs(ih)]*Xpoh*kt2;

                ephi+=jj*kt2/kt2*(-beta/OPdist*(-jj*double(ih)))*hanka[abs(ih)]/hankb[abs(ih)]*Xpoe;

                ephi+=jj*kt2/kt2*w_ua*hankad[abs(ih)]/hankb[abs(ih)]*Xpoh;

                er+=jj*kt2/kt2*(-beta)*hankad[abs(ih)]/hankb[abs(ih)]*Xpoe;

                er+=jj*kt2/kt2*(-w_ua/OPdist)*(-jj*double(ih))*hanka[abs(ih)]/hankb[abs(ih)]*Xpoh;

                hphi+=jj*kt2/kt2*(-w_ea)*hankad[abs(ih)]/hankb[abs(ih)]*Xpoe;

                hphi+=jj*kt2/kt2*(-beta/OPdist*(-jj*double(ih)))*hanka[abs(ih)]/hankb[abs(ih)]*Xpoh;

                hr+=jj*kt2/kt2*(w_ea/OPdist)*(-jj*double(ih))*hanka[abs(ih)]/hankb[abs(ih)]*Xpoe;

                hr+=jj*kt2/kt2*(-beta)*hankad[abs(ih)]/hankb[abs(ih)]*Xpoh;
            }

            ex+=ephi*cos(phi)+er*sin(phi);
            ey+=-ephi*sin(phi)+er*cos(phi);
            hx+=hphi*cos(phi)+hr*sin(phi);
            hy+=-hphi*sin(phi)+hr*sin(phi);

            //ex=ephi;
            //hx=hphi;

            index_coeff+=number_layer*no_harmonics*4;
        }

        return;

    } else{
//----------------------when point is in the wire-------------------------------------------------
        //cout<<"\nPoint in the wire "<<wire_number<<endl;

        const double ddx(x-wires[wire_number].centre[0]);
        const double ddy(y-wires[wire_number].centre[1]);

        const double OPdist(sqrt(ddx*ddx+ddy*ddy));
        const double phi(atan2(ddx,ddy));

        //cout<<"distance from origin: "<<r<<"\nangle: "<<phi<<endl;
    //-----------------point in the first layer---------------------------------------
        int index(0);
        for(int iw=0;iw<wire_number;++iw){

            index+=wires[iw].no_layers*no_harmonics*4;
        }
        //cout<<"number of elements before wire "<<wire_number+1<<": "<< index<<endl;

        if(layer_number==0){

            //cout<<"\nPoint in the first layer."<<endl;
            //cout<<"layer number: "<<layer_number+1<<endl;

            double radius(wires[wire_number].radius[0]);
            //cout<<"layer radius: "<<radius<<endl;
            cout<<"relative epsilon: "<<epsilon_relative<<endl;

            complex kt2_firstlayer(ko*ko*epsilon_relative-beta*beta);
            complex kt_firstlayer(sqrt(kt2_firstlayer));
            //cout<<"kt at first layer: "<<kt_firstlayer<<endl;
            complex w_e0(w_ea*epsilon_relative);

            zc[0]=real(kt_firstlayer*radius);
            zc[1]=imag(kt_firstlayer*radius);

            S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

            for(int k=0;k<no_bes_needed;k++) besb[k]=complex(temp[2*k],temp[2*k+1]);

            zc[0]=real(kt_firstlayer*OPdist);
            zc[1]=imag(kt_firstlayer*OPdist);

            S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

            for(int k=0;k<no_bes_needed;k++) besa[k]=complex(temp[2*k],temp[2*k+1]);

             besad[0]=-kt_firstlayer*besa[1];

            for(int k=1;k<no_bes_needed;k++) besad[k]=kt_firstlayer*(besa[k-1]-double(k)/complex(zc[0],zc[1])*besa[k]);

            ex=ey=ez=ephi=hz=0.;
            hx=hy=er=hphi=hr=0.;

            for(int ih=-max_harmonic;ih<=max_harmonic;++ih){

                complex Xpie(amps[index+(ih+max_harmonic)*(wires[wire_number].no_layers)*4+0]*exp(-jj*double(ih)*phi));
                complex Xpih(amps[index+(ih+max_harmonic)*(wires[wire_number].no_layers)*4+1]*exp(-jj*double(ih)*phi));

                //cout<<"\nXpie: "<<index+(ih+max_harmonic)*wires[wire_number].no_layers*4+0<<endl;

                ez+=besa[abs(ih)]/besb[abs(ih)]*Xpie*kt2;
                hz+=besa[abs(ih)]/besb[abs(ih)]*Xpih*kt2;

                ephi+=jj*kt2/kt2_firstlayer*(-beta/OPdist*(-jj*double(ih)))*besa[abs(ih)]/besb[abs(ih)]*Xpie;
                ephi+=jj*kt2/kt2_firstlayer*w_ua*besad[abs(ih)]/besb[abs(ih)]*Xpih;

                er+=jj*kt2/kt2_firstlayer*(-beta)*besad[abs(ih)]/besb[abs(ih)]*Xpie;
                er+=jj*kt2/kt2_firstlayer*(-w_ua/OPdist)*(-jj*double(ih))*besa[abs(ih)]/besb[abs(ih)]*Xpih;

                hphi+=jj*kt2/kt2_firstlayer*(-w_e0)*besad[abs(ih)]/besb[abs(ih)]*Xpie;
                hphi+=jj*kt2/kt2_firstlayer*(-beta/OPdist*(-jj*double(ih)))*besa[abs(ih)]/besb[abs(ih)]*Xpih;

                hr+=jj*kt2/kt2_firstlayer*(w_e0/OPdist)*(-jj*double(ih))*besa[abs(ih)]/besb[abs(ih)]*Xpie;
                hr+=jj*kt2/kt2_firstlayer*(-beta)*besad[abs(ih)]/besb[abs(ih)]*Xpih;
            }

            ex+=ephi*cos(phi)+er*sin(phi);
            ey+=-ephi*sin(phi)+er*cos(phi);
            hx+=hphi*cos(phi)+hr*sin(phi);
            hy+=-hphi*sin(phi)+hr*sin(phi);
//            ex=ephi;
//            hx=hphi;
            //ex*=epsilon_relative;
            return;

        }else{
            //cout<<"\nPoint in the mid layer."<<endl;
            //cout<<"layer number: "<<layer_number+1<<endl;

            const double radius(wires[wire_number].radius[layer_number]);

            const double radius_mminus1(wires[wire_number].radius[layer_number-1]);
            //cout<<"layer radius: "<<radius<<endl;

            complex kt2m(ko*ko*epsilon_relative-beta*beta);
            complex ktm(sqrt(kt2m));
           // cout<<"kt at layer "<<layer_number+1<<endl;

            complex w_em(w_ea*epsilon_relative);

            zc[0]=real(ktm*OPdist);
            zc[1]=imag(ktm*OPdist);

            S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);
                for(int k=0;k<no_bes_needed;++k){besj[k]=complex(temp[2*k],temp[2*k+1]);}

                besjd[0]=-ktm*besj[1];

                for(int k=1;k<no_bes_needed;++k){besjd[k]=ktm*(besj[k-1]-double(k)/complex(zc[0],zc[1])*besj[k]);}

            S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);
                for(int k=0;k<no_bes_needed;++k){besy[k]=complex(temp[2*k],temp[2*k+1]);}

                besyd[0]=-ktm*besy[1];

                for(int k=1;k<no_bes_needed;++k){besyd[k]=ktm*(besy[k-1]-double(k)/complex(zc[0],zc[1])*besy[k]);}

            zc[0]=real(ktm*radius_mminus1);
            zc[1]=imag(ktm*radius_mminus1);

            S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);
                for(int k=0;k<no_bes_needed;++k){besjm1[k]=complex(temp[2*k],temp[2*k+1]);}

            S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);
                for(int k=0;k<no_bes_needed;++k){besym1[k]=complex(temp[2*k],temp[2*k+1]);}


            ex=ey=ez=ephi=hz=0.;
            hx=hy=er=hphi=hr=0.;

            for(int ih=-max_harmonic;ih<=max_harmonic;++ih){
                complex Xpej(amps[index+(ih+max_harmonic)*wires[wire_number].no_layers*4+layer_number*4-2]*exp(-jj*double(ih)*phi));
                complex Xpey(amps[index+(ih+max_harmonic)*wires[wire_number].no_layers*4+layer_number*4-1]*exp(-jj*double(ih)*phi));
                complex Xphj(amps[index+(ih+max_harmonic)*wires[wire_number].no_layers*4+layer_number*4-0]*exp(-jj*double(ih)*phi));
                complex Xphy(amps[index+(ih+max_harmonic)*wires[wire_number].no_layers*4+layer_number*4+1]*exp(-jj*double(ih)*phi));

                ez+=besj[abs(ih)]/besjm1[abs(ih)]*Xpej*kt2+besy[abs(ih)]/besym1[abs(ih)]*Xpey*kt2;
                hz+=besj[abs(ih)]/besjm1[abs(ih)]*Xphj*kt2+besy[abs(ih)]/besym1[abs(ih)]*Xphy*kt2;

                ephi+=jj*kt2/kt2m*(-beta/OPdist*(-jj*double(ih)))*(besj[abs(ih)]/besjm1[abs(ih)]*Xpej+besy[abs(ih)]/besym1[abs(ih)]*Xpey);
                ephi+=jj*kt2/kt2m*w_ua*(besjd[abs(ih)]/besjm1[abs(ih)]*Xphj+besyd[abs(ih)]/besym1[abs(ih)]*Xphy);

                er+=jj*kt2/kt2m*(-beta)*(besjd[abs(ih)]/besjm1[abs(ih)]*Xpej+besyd[abs(ih)]/besym1[abs(ih)]*Xpey);
                er+=jj*kt2/kt2m*(-w_ua/OPdist)*(-jj*double(ih))*(besj[abs(ih)]/besjm1[abs(ih)]*Xphj+besy[abs(ih)]/besym1[abs(ih)]*Xphy);

                hphi+=jj*kt2/kt2m*(-w_em)*(besjd[abs(ih)]/besjm1[abs(ih)]*Xpej+besyd[abs(ih)]/besym1[abs(ih)]*Xpey);
                hphi+=jj*kt2/kt2m*(-beta/OPdist*(-jj*double(ih)))*(besj[abs(ih)]/besjm1[abs(ih)]*Xphj+besy[abs(ih)]/besym1[abs(ih)]*Xphy);

                hr+=jj*kt2/kt2m*(w_em/OPdist)*(-jj*double(ih))*(besj[abs(ih)]/besjm1[abs(ih)]*Xpej+besy[abs(ih)]/besym1[abs(ih)]*Xpey);
                hr+=jj*kt2/kt2m*(-beta)*(besjd[abs(ih)]/besjm1[abs(ih)]*Xphj+besyd[abs(ih)]/besym1[abs(ih)]*Xphy);
            }

            ex+=ephi*cos(phi)+er*sin(phi);
            ey+=-ephi*sin(phi)+er*cos(phi);
            hx+=hphi*cos(phi)+hr*sin(phi);
            hy+=-hphi*sin(phi)+hr*sin(phi);
            //ex=ephi;
            //hx=hphi;
            return;
        }
    }

    //cout<<"\nephi: "<<ephi<<endl;
}

void plot_field(vector<wire>& wires,
                const double ko,
                const complex beta,
                const int max_harmonic,
                vector <complex>& amps){
    ofstream Exout("Ex.csv");
    ofstream Eyout("Ey.csv");
    ofstream Ezout("Ez.csv");
    ofstream Etout("Et.csv");
    ofstream Hxout("Hx.csv");
    ofstream Hyout("Hy.csv");
    ofstream Hzout("Hz.csv");
    ofstream Htout("Ht.csv");

    int npts(500);

    double xmin(-5e-3),xmax(5e-3);
    double ymin(-5e-3),ymax(5e-3);

    const double dx((xmax-xmin)/double(npts));
    const double dy((ymax-ymin)/double(npts));

    for(int j=0;j<npts;j++) Exout<<","<<ymin+double(j)*dy;
    for(int j=0;j<npts;j++) Eyout<<","<<ymin+double(j)*dy;
    for(int j=0;j<npts;j++) Ezout<<","<<ymin+double(j)*dy;
    for(int j=0;j<npts;j++) Etout<<","<<ymin+double(j)*dy;

    for(int j=0;j<npts;j++) Hxout<<","<<ymin+double(j)*dy;
    for(int j=0;j<npts;j++) Hyout<<","<<ymin+double(j)*dy;
    for(int j=0;j<npts;j++) Hzout<<","<<ymin+double(j)*dy;
    for(int j=0;j<npts;j++) Htout<<","<<ymin+double(j)*dy;

    for(int i=0;i<npts;++i) {
         Exout<<"\n"<<xmin+double(i)*dx;
         Eyout<<"\n"<<xmin+double(i)*dx;
         Ezout<<"\n"<<xmin+double(i)*dx;
         Etout<<"\n"<<xmin+double(i)*dx;

         Hxout<<"\n"<<xmin+double(i)*dx;
         Hyout<<"\n"<<xmin+double(i)*dx;
         Hzout<<"\n"<<xmin+double(i)*dx;
         Htout<<"\n"<<xmin+double(i)*dx;

         complex ex,ey,ez,hx,hy,hz,ephi;

         const double x(xmin+double(i)*dx);

         for(int j=0;j<npts;++j){
            const double y(ymin+double(j)*dy);

            get_fields(wires,ko,beta,max_harmonic,amps,x,y,ex,ey,ez,hx,hy,hz);

            Exout<<","<<abs(ex);
            Eyout<<","<<abs(ey);
            Ezout<<","<<abs(ez);
            Hxout<<","<<abs(hx);
            Hyout<<","<<abs(hy);
            Hzout<<","<<abs(hz);
            Etout<<","<<sqrt(abs(ex)*abs(ex)+abs(ey)*abs(ey));
            Htout<<","<<sqrt(abs(hx)*abs(hx)+abs(hy)*abs(hy));
         }
    }

    Exout.close();
    Eyout.close();
    Ezout.close();
    Etout.close();

    Hxout.close();
    Hyout.close();
    Hzout.close();
    Htout.close();

}
//------------------plot field with respect to angles phi at boundary--------------------------
void plot_field2(vector<wire>& wires,
                const double ko,
                const complex beta,
                const int max_harmonic,
                vector <complex>& amps)
{
    ofstream Ephiout("Ephi.csv");

    int npts(500);

    double xmin(-5e-3),xmax(5e-3);

    const double y(0);

    const double dx((xmax-xmin)/double(npts));

    for(int i=0;i<npts;++i) {

        Ephiout<<"\n"<<xmin+double(i)*dx;

        complex ex,ey,ez,hx,hy,hz,ephi;

        const double x(xmin+double(i)*dx);

        get_fields(wires,ko,beta,max_harmonic,amps,x,y,ex,ey,ez,hx,hy,hz);

        Ephiout<<","<<abs(ex);
    }

    Ephiout.close();
}

/*void plot_field2(vector<wire>& wires,
          const double ko,
          const complex erd,
          const complex erw,
          const complex beta,
          const int max_harmonic,
          vector<complex>& amps){

    ofstream ephiout_phi("ephi_phi.txt");
    ofstream ephiout_r("ephi_r.txt");
    ofstream ezout_phi("ez_phi.txt");
    ofstream hphiout_phi("hphi_phi.txt");

    const complex kod(ko*sqrt(erd));
    const complex kow(ko*sqrt(erw));

    complex kt(sqrt(ko*ko-beta*beta));
    complex kt2(ko*ko-beta*beta);

    complex ktd(sqrt(kod*kod-beta*beta));
    complex ktd2(kod*kod-beta*beta);

    complex ktw(sqrt(kow*kow-beta*beta));
    complex ktw2(kow*kow-beta*beta);

    const double w_uo(ko*constants::get_zo());      //omega*muo for angular field calculation
    const double w_eo(ko*constants::get_yo());      //omega*epsilon_o for angualr field calculation
    const complex w_ed(w_eo*erd);                   //omega*epsilon_rd=omega*epsilon_o*epsilon_rd
    const complex w_ew(w_eo*erw);                   //omega*epsilon_rw=omega*epsilon_o*epsilon_rw


    const int no_harmonics(2*max_harmonic+1);

    int nz,ifail,scale_len=1;

         static complex jj(0.,1.);

         char scale='u';

         double zc[2],fnu=0.;

    int no_bes_needed(std::max(2,max_harmonic+1)),hank_kind(2),hank_its;

    static vector<complex> hanka,hankad,hankb,hankbd,besa,besad,besb,besbd;
    static vector<complex> besy1,besy1d,besj1,besj1d,besj2,besj2d,besy2,besy2d,besj,besjd,besy,besyd;

    static vector<double> temp;
    static vector<double> cwrk;

         hanka.reserve(no_bes_needed);  //Hank(kt*r)
         hankad.reserve(no_bes_needed);

         hankb.reserve(no_bes_needed);  //Hank(kt*r)
         hankbd.reserve(no_bes_needed);

         besa.reserve(no_bes_needed);   //coupling term
         besad.reserve(no_bes_needed);

         besb.reserve(no_bes_needed);   //Bessel (ktw*r1), term inside wires
         besbd.reserve(no_bes_needed);

         besj1.reserve(no_bes_needed);   //term in dielectric
         besj1d.reserve(no_bes_needed);

         besj2.reserve(no_bes_needed);   //term in dielectric
         besj2d.reserve(no_bes_needed);

         besy1.reserve(no_bes_needed);  //Besselktd*r)   term in dielectric
         besy1d.reserve(no_bes_needed);

         besy2.reserve(no_bes_needed);  //Besselktd*r)   term in dielectric
         besy2d.reserve(no_bes_needed);

         besj.reserve(no_bes_needed);  //Besselktd*r)   term in dielectric
         besjd.reserve(no_bes_needed);

         besy.reserve(no_bes_needed);  //Besselktd*r)   term in dielectric
         besyd.reserve(no_bes_needed);


     temp.reserve(2*no_bes_needed);
     cwrk.reserve(2*no_bes_needed);

    const double cradius(wires[0].in_radius);  //conductor radius
    const double radius_out(wires[0].radius);

    int npts(100);

    double phimin(0.),phimax(2*constants::get_pi()),rmin(-2e-3),rmax(2e-3);

    const double dphi((phimax-phimin)/double(npts));
    const double dr((rmax-rmin)/double(npts));
    ephiout_phi<<","<<"ephi1"<<","<<"ephi21"<<","<<"ephi22"<<","<<"ephi3";
    ephiout_r<<","<<"ephi";

    ezout_phi<<","<<"real(ez1)"<<","<<"imag(ez1)"<<","<<"ez21_real"<<","<<"ephiz_imag";

    hphiout_phi<<","<<"hphi1"<<","<<"hphi21"<<","<<"hphi22"<<","<<"hphi3";
    for(int i=0;i<npts;i++){

        ephiout_phi<<"\n"<<phimin+double(i)*dphi;
        ezout_phi<<"\n"<<phimin+double(i)*dphi;
        hphiout_phi<<"\n"<<phimin+double(i)*dphi;
        const double phi(phimin+double(i)*dphi);

        zc[0]=real(kt*radius_out);
        zc[1]=imag(kt*radius_out);
                                            //argument for bessel and hankel at out boundary
        S17DLF(&hank_kind,&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&hank_its,&ifail);

            for(int k=0;k<no_bes_needed;k++) {
                hanka[k]=complex(temp[2*k],temp[2*k+1]);}

            hankad[0]=-kt*hanka[1];
                                            //derivative of hankel for n=0;
            for(int k=1;k<no_bes_needed;k++) {
                hankad[k]=kt*(hanka[k-1]-double(k)/complex(zc[0],zc[1])*hanka[k]);}

        zc[0]=real(kt*radius_out);
        zc[1]=imag(kt*radius_out);

        S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

            for(int k=0;k<no_bes_needed;k++) besa[k]=complex(temp[2*k],temp[2*k+1]);

            besad[0]=-ktw*besa[1];
            for(int k=1;k<no_bes_needed;k++) besad[k]=kt*(besa[k-1]-double(k)/complex(zc[0],zc[1])*besa[k]);

        zc[0]=real(ktw*cradius);
        zc[1]=imag(ktw*cradius);

        S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

            for(int k=0;k<no_bes_needed;k++) besb[k]=complex(temp[2*k],temp[2*k+1]);

            besbd[0]=-ktw*besb[1];
            for(int k=1;k<no_bes_needed;k++) besbd[k]=ktw*(besb[k-1]-double(k)/complex(zc[0],zc[1])*besb[k]);

        zc[0]=real(ktd*cradius);
        zc[1]=imag(ktd*cradius);

        S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

            for(int k=0;k<no_bes_needed;k++) {besj1[k]=complex(temp[2*k],temp[2*k+1]);}

            besj1d[0]=-ktd*besj1[1];
            for(int k=1;k<no_bes_needed;k++) besj1d[k]=ktd*(besj1[k-1]-double(k)/complex(zc[0],zc[1])*besj1[k]);

        S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);

            for(int k=0;k<no_bes_needed;k++) {besy1[k]=complex(temp[2*k],temp[2*k+1]);}

            besy1d[0]=-ktd*besy1[1];

            for(int k=1;k<no_bes_needed;k++) {besy1d[k]=ktd*(besy1[k-1]-double(k)/complex(zc[0],zc[1])*besy1[k]);}

        zc[0]=real(ktd*wires[0].radius);
        zc[1]=imag(ktd*wires[0].radius);
               //cout<<"dielectric zc: "<<zc[0]+jj*zc[1]<<endl;
        S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

            for(int k=0;k<no_bes_needed;k++) {besj2[k]=complex(temp[2*k],temp[2*k+1]);}

            besj2d[0]=-ktd*besj2[1];
            for(int k=1;k<no_bes_needed;k++) besj2d[k]=ktd*(besj2[k-1]-double(k)/complex(zc[0],zc[1])*besj2[k]);

        S17DCF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&cwrk[0],&ifail);

            for(int k=0;k<no_bes_needed;k++) {besy2[k]=complex(temp[2*k],temp[2*k+1]);}

            besy2d[0]=-ktd*besy2[1];
            for(int k=1;k<no_bes_needed;k++) besy2d[k]=ktd*(besy2[k-1]-double(k)/complex(zc[0],zc[1])*besy2[k]);

         complex ephi1(0.),ephi21(0.),hphi1(0.),hphi21(0.),ez1(0.),ez21(0.),hz1(0.),hz21(0.);
         complex ephi_out(0.),hphi_out(0.),ephi_dielectric_out(0.),hphi_dielectric_out(0.);
         for (int ia=-max_harmonic;ia<=max_harmonic;ia++){

            complex Xpie_over_kt2(amps[0*4*no_harmonics+4*(ia+max_harmonic)+2]*exp(-jj*double(ia)*phi));
            complex Xpih_over_kt2(amps[0*4*no_harmonics+4*(ia+max_harmonic)+3]*exp(-jj*double(ia)*phi));

            complex Xpe_over_kt2(amps[0*4*no_harmonics+4*(ia+max_harmonic)+0]*exp(-jj*double(ia)*phi));
            complex Xph_over_kt2(amps[0*4*no_harmonics+4*(ia+max_harmonic)+1]*exp(-jj*double(ia)*phi));
            //cout<<"\nXpe term"<<0*4*no_harmonics+4*(ia+max_harmonic)+0<<": "<<Xpe_over_kt2;

            complex TXqe(0.,0.), TXqh(0.,0.);

            wire_coupling_termer T(wires[0],wires[1],max_harmonic,kt);

            complex alpha(besj1[abs(ia)]*besy2[abs(ia)]-besj2[abs(ia)]*besy1[abs(ia)]);

            for(int ja=-max_harmonic;ja<=max_harmonic;ja++){

                    complex Xqe_over_kt2(amps[1*4*no_harmonics+4*(ja+max_harmonic)+0]*exp(-jj*double(ia)*phi));
                    complex Xqh_over_kt2(amps[1*4*no_harmonics+4*(ja+max_harmonic)+1]*exp(-jj*double(ia)*phi));
                    //cout<<"\nXqe term"<<1*4*no_harmonics+4*(ja+max_harmonic)+0<<": "<<Xqe_over_kt2;
                TXqe+=Xqe_over_kt2*T(ia,ja);
                TXqh+=Xqh_over_kt2*T(ia,ja);

            }

            complex coeffAe((Xpie_over_kt2*besy2[abs(ia)]-(Xpe_over_kt2+TXqe)*besy1[abs(ia)])/alpha);
            complex coeffBe(((Xpe_over_kt2+TXqe)*besj1[abs(ia)]-Xpie_over_kt2*besj2[abs(ia)])/alpha);

            complex coeffAh((Xpih_over_kt2*besy2[abs(ia)]-(Xph_over_kt2+TXqh)*besy1[abs(ia)])/alpha);
            complex coeffBh(((Xph_over_kt2+TXqh)*besj1[abs(ia)]-Xpih_over_kt2*besj2[abs(ia)])/alpha);

             ez1+=Xpie_over_kt2;

             ez21+=besj1[abs(ia)]*coeffAe+besy1[abs(ia)]*coeffBe;

             hz1+=Xpih_over_kt2;

             hz21+=besj1[abs(ia)]*coeffAh+besy1[abs(ia)]*coeffBh;

             ephi1+=jj*(kt2/ktw2)*(-beta/cradius*(-jj*double(ia)))*Xpie_over_kt2;

             ephi1+=jj*(kt2/ktw2)*w_uo*besbd[abs(ia)]/besb[abs(ia)]*Xpih_over_kt2;

             ephi21+=jj*(kt2/ktd2)*(-beta/cradius*(-jj*double(ia)))*(coeffAe*besj1[abs(ia)]+coeffBe*besy1[abs(ia)]);

             ephi21+=jj*(kt2/ktd2)*w_uo*(coeffAh*besj1d[abs(ia)]+coeffBh*besy1d[abs(ia)]);

             ephi_dielectric_out+=jj*(kt2/ktd2)*(-beta/radius_out*(-jj*double(ia)))*(coeffAe*besj2[abs(ia)]+coeffBe*besy2[abs(ia)]);

             ephi_dielectric_out+=jj*(kt2/ktd2)*w_uo*(coeffAh*besj2d[abs(ia)]+coeffBh*besy2d[abs(ia)]);

             ephi_out+=jj*(-beta/radius_out*(-jj*double(ia)))*(Xpe_over_kt2+TXqe);

             ephi_out+=jj*w_uo*(hankad[abs(ia)]/hanka[abs(ia)]*Xph_over_kt2+besad[abs(ia)]/besa[abs(ia)]*TXqh);


             hphi1+=jj*(kt2/ktw2)*(-w_ew)*besbd[abs(ia)]/besb[abs(ia)]*Xpie_over_kt2;

             hphi1+=jj*(kt2/ktw2)*(-beta/cradius*(-jj*double(ia)))*Xpih_over_kt2;

             hphi21+=jj*(kt2/ktd2)*(-w_ed)*(coeffAe*besj1d[abs(ia)]+coeffBe*besy1d[abs(ia)]);

             hphi21+=jj*(kt2/ktd2)*(-beta/cradius*(-jj*double(ia)))*(coeffAh*besj1[abs(ia)]+coeffBh*besy1[abs(ia)]);

             hphi_dielectric_out+=jj*(kt2/ktd2)*(-w_ed)*(coeffAe*besj2d[abs(ia)]+coeffBe*besy2d[abs(ia)]);

             hphi_dielectric_out+=jj*(kt2/ktd2)*(-beta/cradius*(-jj*double(ia)))*(coeffAh*besj2[abs(ia)]+coeffBh*besy2[abs(ia)]);

             hphi_out+=jj*(-w_eo)*(hankad[abs(ia)]/hanka[abs(ia)]*Xpe_over_kt2+besad[abs(ia)]/besa[abs(ia)]*TXqe);

             hphi_out+=jj*(-beta/radius_out*(-jj*double(ia)))*(Xph_over_kt2+TXqh);
         }

         ephiout_phi<<","<<abs(ephi1)<<","<<abs(ephi21)<<","<<abs(ephi_dielectric_out)<<","<<abs(ephi_out);

         //ezout_phi<<","<<real(ez1)<<","<<imag(ez1)<<","<<real(ez21)<<","<<imag(ez21);

         hphiout_phi<<","<<abs(hphi1)<<","<<abs(hphi21)<<","<<abs(hphi_dielectric_out)<<","<<abs(hphi_out);
//-----------------ephi plot with respect to phi end------------------------------------------

        const double x(rmin+double(i)*dr);
        const double y(0.);
        ephiout_r<<"\n"<<x;
        complex ex,ey,ez,hx,hy,hz,ephi;
        get_fields(wires,ko,erd,erw,beta,max_harmonic,amps,x,y,ex,ey,ez,hx,hy,hz,ephi);
        ephiout_r<<","<<abs(ephi);
    }
}*/

