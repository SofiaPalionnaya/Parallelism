#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <omp.h>

using namespace std;

double func(double x, double k, double l, double m)
{
    return pow(x, 3) + k * pow(x, 2) + l * x + m;
}

double find(double infinum, double supremum, double epsilon, double k, double l, double m)
{
    while (fabs(supremum - infinum) > epsilon)
    {
        infinum = supremum - (supremum - infinum) * func(supremum,k,l,m) / (func(supremum,k,l,m) - func(infinum,k,l,m));
        supremum = infinum - (infinum - supremum) * func(infinum,k,l,m) / (func(infinum,k,l,m) - func(supremum,k,l,m));
    }
    if (supremum > 1)
    {
        return 1;
    }
    else if (supremum < 0)
    {
        return 0;
    }
    else
    {
        return supremum;
    }
}

double parab_random_variable(long double a, long double b)
{
    long double r = rand();
    r = r/RAND_MAX;
    long double denom;
    denom = pow(a,3) - pow(b,3) - 3*b*pow(a,2) + 3*a*pow(b,2);
    if (denom == 0){
        return 0;
    }
    long double S = -2*pow(a,3)/denom - pow(a,2)*(-3*a - 3*b)/denom - 6*b*pow(a,2)/denom - r;
    long double k = (-3*a - 3*b)/2;
    long double l = 3*a*b;
    long double m = S * denom / 2;
    return(find(a, b, 0.000000001, k, l, m));
}

int main()
{
    srand( time( 0 ) );
    long double my_time;
    int kol;
    kol = 100000;
    int kol2 = kol;
    long double *p = new long double[kol];
    my_time = 1000000;
    my_time = my_time*100000;

    long double a_l, b_l, a_m, b_m, lambda, mu;
    long double eps = 0.025;
    a_l = 0;
    b_l = 1;
    a_m = 0;
    b_m = 1;

    if ((b_l<a_l)or(b_m<a_m))
    {
        cout << "неверные входные параметры"<< endl;
        exit;
    }
    char h[kol];
    #pragma omp parallel for
    for(int i=0; i<kol; i++)
    {
        lambda = parab_random_variable(a_l,b_l);
        mu = parab_random_variable(a_m,b_m);
        if ((lambda == 7) or (mu == 7))
        {
            kol2 = kol2 - 1;
            p[i] = 0;
            continue;
        }
        long double b_t = 1;
        long double b_e = 1;
        long double a_e, a_t;

        a_t = 1-lambda-eps;
        a_e = mu-eps;
        b_t = a_t+eps;
        b_e = a_e+eps;

        if (a_t<0)
        {  a_t = 0; }
        if (a_e<0)
        {  a_e = 0; }
        if (b_t>1)
        {  b_t = 1; }
        if (b_e>1)
        {  b_e = 1; }

        p[i] = rand();
        p[i] = double(p[i]/RAND_MAX);
        long double etta;
        long double tetta;
        long double t = 0;
        do
        {
            t += rand();
            etta = rand();
            etta = (etta*(b_e - a_e))/RAND_MAX + a_e;
            tetta = rand();
            tetta = (tetta*(b_t - a_t))/RAND_MAX + a_t;
            p[i] = etta*p[i] + tetta*(1-p[i]);
        }
        while(t<my_time);
    }

    long double sr_nad = 0;
    long double sum = 0;
    for (int i=0; i<kol; i++)
    {
        if (p[i]<0)
        {p[i]=0;}
        sum = sum + p[i];
    }
    sr_nad = sum/kol2;
    cout << "aver_reliability= " << sr_nad << endl;
    cout << "clock= "<< clock()/1000000.0<< endl;
    delete []p;
    return 0;
}
