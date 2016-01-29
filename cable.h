#ifndef CABLE_H_INCLUDED
#define CABLE_H_INCLUDED

//------------include--------------------------------------------------

#include <string.h>
#include <complex>
#include <vector>
#include <math.h>

#include "nagmk18.h"
//-------------------------------------------------------------------

#define complex complex<double>

using namespace std;

extern "C" void  dgesdd_(const char* jobz,
                         const int * m,
                         const int * n,
                         double * a,
                         const int *lda,
                         double * s,
                         double * u,
                         const int *ldu,
                         double * vt,
                         const int * ldvt,
                         double * work,
                         const int * lwork,
             int* iwork,
                         int* info);

struct constants {

    static const double get_pi(){return(M_PI);}
    static const double get_uo(){return(4.e-7*get_pi());}
    static const double get_eo(){return(8.8541878176e-12);}
    static const double get_co(){return(1./sqrt(get_eo()*get_uo()));}
    static const double get_yo(){return(sqrt(get_eo()/get_uo()));}
    static const double get_zo(){return(sqrt(get_uo()/get_eo()));}

 };

//-----------------------------------------------------------------------------------------
struct wire {
     double centre[2];

     double *radius;

     int no_layers;

     complex *epsilonr;

     wire(const double* const _centre,const double *_radius, const complex *_epsilonr, const int _no_layers)
             :no_layers(_no_layers){

         memcpy(centre,_centre,2*sizeof(double));

         radius=new double[no_layers];
         for(int i=0;i<no_layers;i++){

            radius[i]=_radius[i];
         }

         epsilonr=new complex[no_layers];
         for(int i=0;i<no_layers;i++){

            epsilonr[i]=_epsilonr[i];
         }
     }

     wire(const double x, const double y,const double* _radius, const complex* _epsilonr, const int _no_layers)
            :no_layers(_no_layers){

         centre[0]=x;centre[1]=y;

         radius=new double[no_layers];
         for(int i=0;i<no_layers;i++){

            radius[i]=_radius[i];
         }

         epsilonr=new complex[no_layers];
         for(int i=0;i<no_layers;i++){

            epsilonr[i]=_epsilonr[i];
         }

     }

     wire(const wire& w)
             :no_layers(w.no_layers){

         memcpy(centre,w.centre,2*sizeof(double));

         radius=new double[no_layers];
         for(int i=0;i<no_layers;i++){

            radius[i]=w.radius[i];
         }

         epsilonr=new complex[no_layers];
         for(int i=0;i<no_layers;i++){

            epsilonr[i]=w.epsilonr[i];
         }
     }

     ~wire(){
         delete[] radius;
         delete[] epsilonr;
     }

     const wire& operator=(const wire& w) {

         if(this==&w) return(*this);

         no_layers=w.no_layers;

         memcpy(centre,w.centre,2*sizeof(double));

         radius=new double[no_layers];
         for(int i=0;i<no_layers;i++){

            radius[i]=w.radius[i];
         }

         epsilonr=new complex[no_layers];
         for(int i=0;i<no_layers;i++){

            epsilonr[i]=w.epsilonr[i];
         }

         return(*this);
     }

 };

//----------------------------------------------------------------------------------------

 struct wire_coupling_termer {

     const wire& from;

     const wire& to;

     const int max_harmonic;

     static vector<complex> T;

     static vector<complex> hank_from;

     static vector<complex> bes_to;

     static vector<double> temp;

        wire_coupling_termer(const wire& _to,const wire& _from,const int _max_harmonic,const complex kt)

             :from(_from),to(_to),max_harmonic(_max_harmonic){

         int nz,scale_len=1;

         static const complex jj(0.,1.);

         char scale='u';

         const double fnu(0.);

         int no_bes(2*max_harmonic+1),hank_kind(2),hank_its,ifail(0);

         T.reserve(no_bes*no_bes);

         const double dx(to.centre[0]-from.centre[0]);

         const double dy(to.centre[1]-from.centre[1]);

 //    y
 //    ^ phi - the angle to y axis, tan phi = x/y
 //    | /
 //    |/
 //    -----------> x
 //
         const double sig(sqrt(dx*dx+dy*dy));  //distance

         const double gampq(atan2(dx,dy));  //angle

         complex mjgampq=-jj*gampq;

        double zc[2];

         const int no_bes_needed(std::max(2,2*max_harmonic+1));

         temp.reserve(2*no_bes_needed);

                zc[0]=real(kt*from.radius[from.no_layers-1]);

                zc[1]=imag(kt*from.radius[from.no_layers-1]);

                 S17DLF(&hank_kind,&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&hank_its,&ifail);

         hank_from.reserve(no_bes_needed+1);

                 for(int k=0;k<no_bes_needed;k++) {hank_from[k]=complex(temp[2*k],temp[2*k+1]);}

         zc[0]=real(kt*to.radius[to.no_layers-1]);

         zc[1]=imag(kt*to.radius[to.no_layers-1]);
            //cout<<"\nto.radius: "<<to.radius[to.no_layers-1]<<endl;
            //cout<<"\nzc: "<<zc[0]<<" "<<zc[1]<<endl;

            S17DEF(&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&nz,&ifail);

            bes_to.reserve(no_bes_needed+1);

                 for(int k=0;k<no_bes_needed;k++) bes_to[k]=complex(temp[2*k],temp[2*k+1]);

         zc[0]=real(kt)*sig;

         zc[1]=imag(kt)*sig;


         S17DLF(&hank_kind,&fnu,zc,&no_bes_needed,&scale,scale_len,&temp[0],&hank_its,&ifail);

         for(int i=-max_harmonic;i<=max_harmonic;++i) { // to

             for(int j=-max_harmonic;j<=max_harmonic;++j) { // from

                 int m=j-i;

                 int abs_m=abs(m);

                 double sgn(1.);

                if(m<0&&(abs_m%2!=0)) sgn*=-1.; // Z<-n> = -Z<n> for odd n

                 if(i<0&&(abs(i)%2!=0)) sgn*=-1.;

                 if(j<0&&(abs(j)%2!=0)) sgn*=-1.;

                 T[(i+max_harmonic)+(j+max_harmonic)*(2*max_harmonic+1)]
                     =complex(temp[2*abs_m],temp[2*abs_m+1])*sgn*exp(double(m)*mjgampq)*bes_to[abs(i)]/hank_from[abs(j)];

                 //cout<<"Tnm: "<<T[(i+max_harmonic)+(j+max_harmonic)*(2*max_harmonic+1)]<<endl;

             }
         }
     }



     complex operator()(const int to_harmonic,const int from_harmonic){

         return(T[(to_harmonic+max_harmonic)+(from_harmonic+max_harmonic)*(2*max_harmonic+1)]);

     }

 };

 vector<complex> wire_coupling_termer::T;
 vector<complex> wire_coupling_termer::hank_from;
 vector<complex> wire_coupling_termer::bes_to;
 vector<double>  wire_coupling_termer::temp;

int get_determinant(vector<wire>& wires,
                    const double ko,
                    const complex beta,
                    const int max_harmonic,
                    vector<complex>& amps,
                    int which_amps,
                    vector<double>& sv);

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
                complex& hz,
                complex& ephi);

void plot_field(vector<wire>& wires,
                const double ko,
                const complex beta,
                const int max_harmonic,
                vector <complex>& amps);

void plot_field2(vector<wire>& wires,
                  const double ko,
                  const complex erd,
                  const complex erw,
                  const complex beta,
                  const int max_harmonic,
                  vector<complex>& amps);
#endif // CABLE_H_INCLUDED
