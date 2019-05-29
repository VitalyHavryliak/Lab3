#include <omp.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <iostream>
using namespace std;

int main() {
    int i, j, n;
    double a, b, p=0.0002, x, s, h, z;
    double t1, t2, t[5];

    ofstream fout("text_lab3.txt");

    a=0; b=2;
    //a=0.5; b=1.5;

//system("pause");
//____________________________________________________________
    t1 = omp_get_wtime();

    h=p;
    n=(b-a)/h;
    s=(a*a/pow(6.1+a,2)*pow(3.7+a,2)+a*a/pow(6.1+b,2)*pow(3.7+b,2)/2);
    //s=(a*exp(0.7*a)*cos(0.1*a)+b*exp(0.7*b)*cos(0.1*b)/2);

    for (i=1; i<=n; i++){
        x=a+h*i;
        z=x*x/pow(6.1+x,2)*pow(3.7+x,2);
        //z=x*exp(0.7*x)*cos(0.1*x);
        s=s+z;
    }
    s=s*h;
    cout<<s<<endl;

    t2 = omp_get_wtime();
    t[0]=t2-t1;
//____________________________________________________
    for (int k=1; k<=3; k++){
        int kk=pow(2, k-1);

        h=p;
        n=(b-a)/h;
        s=(a*a/pow(6.1+a,2)*pow(3.7+a,2)+a*a/pow(6.1+b,2)*pow(3.7+b,2)/2);
        //s=(a*exp(0.7*a)*cos(0.1*a)+b*exp(0.7*b)*cos(0.1*b)/2);
        t1 = omp_get_wtime();
#pragma omp parallel for num_threads(kk) private(x, z)
        for (i=1; i<=n; i++){
            x=a+h*i;
            z=x*x/pow(6.1+x,2)*pow(3.7+x,2);
            //z=x*exp(0.7*x)*cos(0.1*x);

#pragma omp atomic
            s=s+z;
        }
        s=s*h;
        t2 = omp_get_wtime();

        cout<<s<<endl;


        t[k]=t2-t1;
    }
//______________________________________________________________
    omp_lock_t lock;
    omp_init_lock(&lock);

#pragma omp parallel
    {
        omp_set_lock(&lock);
        fout<<"The beginning of the closed section (thread "<<omp_get_thread_num()<<")"<<endl;
        fout<<"The end of the closed section (thread "<<omp_get_thread_num()<<")"<<endl;
        omp_unset_lock(&lock);
    }
    omp_destroy_lock(&lock);

#pragma omp parallel
    {
#pragma omp master
        {
            fout<<"\nLab. robota #3"<<endl;
            //four<<<<endl;
            fout<<"KN-318"<<endl;
            fout<<"Havryliak V.T."<<endl;
            fout<<"Variant 5"<<endl;
        }
    }

//______________________________________________________________
    cout<<"\nt0 = "<<t[0]<<endl;
    cout<<"\nt1 = "<<t[1]<<"\np1 = "<<t[0]/t[1]<<"\ne1 = "<<t[0]/t[1]/1<<endl;
    cout<<"\nt2 = "<<t[2]<<"\np2 = "<<t[0]/t[2]<<"\ne2 = "<<t[0]/t[2]/2<<endl;
    cout<<"\nt4 = "<<t[3]<<"\np4 = "<<t[0]/t[3]<<"\ne4 = "<<t[0]/t[3]/4<<endl;
//	cout<<"\nt8 = "<<t[4]<<"\np8 = "<<t[0]/t[4]<<"\ne8 = "<<t[0]/t[4]/8<<endl;
}

// x*x/(0.1*x*x+0.2*x+0.3)
// x/(sin(0.95*x)*cos(0.95*x)
//		omp_lock_t lock;
//		omp_init_lock(&lock);
//		omp_set_lock(&lock);
//		omp_unset_lock(&lock);
//		omp_destroy_lock(&lock);