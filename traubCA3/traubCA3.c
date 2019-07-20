#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NN 1
#define NC 19
#define DT 0.05 // 10 mu-sec
#define T 1000 // ms
#define NT 40000 // T/DT

#define Cm (3.0) // mu-F/cm^2
#define Ri (0.1) // K Ohm-cm
#define Rm (10.0) // K Ohm-cm^2 (unused)
#define Beta (0.075)

#define V_leak (-60.0)
#define V_Na ((115.0+(V_leak)))
#define V_Ca ((140.0+(V_leak)))
#define V_K ((-15.0+(V_leak)))

double section(const double radius) { return M_PI*radius*radius; }

const double g_Na[NC] = {0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 15.0, 30.0, 15.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const double g_K_DR[NC] = {0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 5.0, 15.0, 5.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const double g_K_A[NC] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const double g_K_C[NC] = {0.0, 5.0, 5.0, 10.0, 10.0, 10.0, 5.0, 20.0, 10.0, 20.0, 5.0, 15.0, 15.0, 15.0, 15.0, 15.0, 5.0, 5.0, 0.0};
const double g_K_AHP[NC] = {0.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.0};
const double g_Ca[NC] = {0.0, 5.0, 5.0, 12.0, 12.0, 12.0, 5.0, 8.0, 4.0, 8.0, 5.0, 17.0, 17.0, 17.0, 10.0, 10.0, 5.0, 5.0, 0.0};
const double g_leak[NC] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
const double phi[NC] = {7769, 7769, 7769, 7769.0, 7769.0, 7769.0, 7769.0, 34530.0, 17402.0, 26404.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0};
const double radius[NC] = {2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 4.23e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4};
const double length[NC] = {1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.25e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2};
const double area[NC] = {2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 3.320e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5};

double v[NC], xi[NC];
double m[NC], s[NC], n[NC], c[NC], a[NC], h[NC], r[NC], b[NC], q[NC];
double I_inj[NC];
double g_compartment[NC][2];

double alpha_m(const double v){return 0.32*(13.1-(v-V_leak))/(exp((13.1-(v-V_leak))/4.0)-1);}
double alpha_s(const double v){return 1.6/(1+exp(-0.072*((v-V_leak)-65)));}
double alpha_n(const double v){return 0.016*(35.1-(v-V_leak))/(exp((35.1-(v-V_leak))/5.0)-1);}
double alpha_c(const double v){return (v<=(50+V_leak))? exp(((v-V_leak)-10)/11.0 - ((v-V_leak)-6.5)/27.0)/18.975 : 2*exp(-((v-V_leak)-6.5)/27.0);}
double alpha_a(const double v){return 0.02*(13.1-(v-V_leak))/(exp((13.1-(v-V_leak))/10.0)-1);}
double alpha_h(const double v){return 0.128*exp((17-(v-V_leak))/18.0);}
double alpha_r(const double v){return (v<=(0+V_leak))? 0.005 : exp(-(v-V_leak)/20.0)/200.0;}
double alpha_b(const double v){return 0.0016*exp((-13-(v-V_leak))/18.0);}
double alpha_q(const double x){return fmin((0.2e-4)*x, 0.01);}
double beta_m(const double v){return 0.28*((v-V_leak)-40.1)/(exp(((v-V_leak)-40.1)/5.0)-1);}
double beta_s(const double v){return 0.02*((v-V_leak)-51.1)/(exp(((v-V_leak)-51.1)/5.0)-1);}
double beta_n(const double v){return 0.25*exp((20-(v-V_leak))/40.0);}
double beta_c(const double v){return (v<=(50+V_leak))? 2*exp(-((v-V_leak)-6.5)/27.0) - alpha_c(v) : 0;}
double beta_a(const double v){return 0.0175*((v-V_leak)-40.1)/(exp(((v-V_leak)-40.1)/10.0)-1);}
double beta_h(const double v){return 4.0/(1+exp((40-(v-V_leak))/5.0));}
double beta_r(const double v){return (v<=(0+V_leak))? 0 : 0.005-alpha_r(v);}
double beta_b(const double v){return 0.05/(1+exp((10.1-(v-V_leak))/5.0));}
double beta_q(const double x){return 0.001;}
double tau_m(const double v){return 1.0/(alpha_m(v)+beta_m(v));}
double tau_s(const double v){return 1.0/(alpha_s(v)+beta_s(v));}
double tau_n(const double v){return 1.0/(alpha_n(v)+beta_n(v));}
double tau_c(const double v){return 1.0/(alpha_c(v)+beta_c(v));}
double tau_a(const double v){return 1.0/(alpha_a(v)+beta_a(v));}
double tau_h(const double v){return 1.0/(alpha_h(v)+beta_h(v));}
double tau_r(const double v){return 1.0/(alpha_r(v)+beta_r(v));}
double tau_b(const double v){return 1.0/(alpha_b(v)+beta_b(v));}
double tau_q(const double x){return 1.0/(alpha_q(x)+beta_q(x));}
double inf_m(const double v){return alpha_m(v)/(alpha_m(v)+beta_m(v));}
double inf_s(const double v){return alpha_s(v)/(alpha_s(v)+beta_s(v));}
double inf_n(const double v){return alpha_n(v)/(alpha_n(v)+beta_n(v));}
double inf_c(const double v){return alpha_c(v)/(alpha_c(v)+beta_c(v));}
double inf_a(const double v){return alpha_a(v)/(alpha_a(v)+beta_a(v));}
double inf_h(const double v){return alpha_h(v)/(alpha_h(v)+beta_h(v));}
double inf_r(const double v){return alpha_r(v)/(alpha_r(v)+beta_r(v));}
double inf_b(const double v){return alpha_b(v)/(alpha_b(v)+beta_b(v));}
double inf_q(const double x){return alpha_q(x)/(alpha_q(x)+beta_q(x));}

void initialize(void)
{
  for(int i = 0; i < NC; i++){
    v[i] = V_leak;
    m[i] = inf_m(V_leak);
    s[i] = inf_s(V_leak);
    n[i] = inf_n(V_leak);
    c[i] = inf_c(V_leak);
    a[i] = inf_a(V_leak);
    h[i] = inf_h(V_leak);
    r[i] = inf_r(V_leak);
    b[i] = inf_b(V_leak);
    double i_Ca = (g_Ca[i]*area[i])*s[i]*s[i]*r[i]*(v[i]-V_Ca);
    xi[i] = -i_Ca*phi[i]/Beta;
    q[i] = inf_q(xi[i]);
    I_inj[i] = 0;
  }
  //conductance between compartments
  //No.0
  g_compartment[0][0] = 0.0;
  g_compartment[0][1] = 2.0 / ((Ri*length[0])/section(radius[0]) + (Ri*length[1])/section(radius[1]));
  //No. 1-17
  for(int i = 0; i < 17; i++){
    g_compartment[i+1][0] = 2.0 / ((Ri*length[i])/section(radius[i]) + (Ri*length[i+1])/section(radius[i+1]));
    g_compartment[i+1][1] = 2.0 / ((Ri*length[i+1])/section(radius[i+1]) + (Ri*length[i+2])/section(radius[i+2]));
  }
  //No.18
  g_compartment[18][0] = 2.0 / ((Ri*length[17])/section(radius[17]) + (Ri*length[18])/section(radius[18]));
  g_compartment[18][1] = 0;
}
void euler(void)
{
  double I_ionic[NC], dv[NC], dxi[NC];
  double dm[NC], ds[NC], dn[NC], dc[NC], da[NC], dh[NC], dr[NC], db[NC], dq[NC];

  for(int i = 0; i < NC; i++){
    dm[i] = (DT/tau_m(v[i]))*(-m[i] + inf_m(v[i]));
    ds[i] = (DT/tau_s(v[i]))*(-s[i] + inf_s(v[i]));
    dn[i] = (DT/tau_n(v[i]))*(-n[i] + inf_n(v[i]));
    dc[i] = (DT/tau_c(v[i]))*(-c[i] + inf_c(v[i]));
    da[i] = (DT/tau_a(v[i]))*(-a[i] + inf_a(v[i]));
    dh[i] = (DT/tau_h(v[i]))*(-h[i] + inf_h(v[i]));
    dr[i] = (DT/tau_r(v[i]))*(-r[i] + inf_r(v[i]));
    db[i] = (DT/tau_b(v[i]))*(-b[i] + inf_b(v[i]));
    dq[i] = (DT/tau_q(xi[i]))*(-q[i] + inf_q(xi[i]));
    I_ionic[i] = (g_leak[i])*(v[i]-V_leak)
      +(g_Na[i])*m[i]*m[i]*h[i]*(v[i]-V_Na)
      +(g_Ca[i])*s[i]*s[i]*r[i]*(v[i]-V_Ca)
      +(g_K_DR[i])*n[i]*(v[i]-V_K)
      +(g_K_A[i])*a[i]*b[i]*(v[i]-V_K)
      +(g_K_AHP[i])*q[i]*(v[i]-V_K)
      +(g_K_C[i])*c[i]*fmin(1,xi[i]/250.0)*(v[i]-V_K)
      -I_inj[i];
    if (i == 0){
      dv[i] = (DT/Cm)*(-I_ionic[i]                                             + g_compartment[i][1]*(v[i+1]-v[i])/area[i]);
    }
    if (1 <= i && i < 18){
      dv[i] = (DT/Cm)*(-I_ionic[i] + g_compartment[i][0]*(v[i-1]-v[i])/area[i] + g_compartment[i][1]*(v[i+1]-v[i])/area[i]);
    }
    if (i == 18){
      dv[i] = (DT/Cm)*(-I_ionic[i] + g_compartment[i][0]*(v[i-1]-v[i])/area[i]);
    }
    double i_Ca = (g_Ca[i]*area[i])*s[i]*s[i]*r[i]*(v[i]-V_Ca);
    dxi[i] = DT*(-phi[i]*i_Ca - Beta*xi[i]);
  }
  for(int i = 0; i < NC; i++){
    v[i] += dv[i];
    xi[i] += dxi[i];
    m[i] += dm[i];
    s[i] += ds[i];
    n[i] += dn[i];
    c[i] += dc[i];
    a[i] += da[i];
    h[i] += dh[i];
    r[i] += dr[i];
    b[i] += db[i];
    q[i] += dq[i];
  }
}

void loop(void)
{
  for(int nt = 0; nt < NT; nt++){
    double t = DT*nt;
    I_inj[8] = 0.1e-3/area[8];
    printf("%f %f\n", t, v[8]);
    euler();
  }
}
int main(void)
{
  initialize();
  loop();

  return 0;
}

