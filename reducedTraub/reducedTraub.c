/*
  Title: An implementation of the reduced Traub model.

  Article:
  Paul F Pinsky, John Rinzel (1994) Intrinsic and network rhythmogenesis in a reduce Traub model
  for CA3 neurons. J Comput Neurosci 1, 39-60.

  Note:
  - Written in C99.
  - Potential numerical instability could be due to the forward Euler method.
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>

#define NN 100 // # of neurons
#define NS 2 // # of segments per neuron
#define DT 0.005 // 0.05 // Time step (ms)
#define T 1000.0 // Duration of simulation (ms)
#define NT ((int)((T)/(DT))) // # of steps 

#define Cm (3.0) // Capacitance (microF/cm^2) p59
#define gbar_L (0.1) // Maximal conductance (mS/cm^2)
#define gbar_Na (30) // Maximal conductance (mS/cm^2)
#define gbar_K_DR (15) // Maximal conductance (mS/cm^2)
#define gbar_Ca (10) // Maximal conductance (mS/cm^2)
#define gbar_K_AHP (0.8) // Maximal conductance (mS/cm^2)
#define gbar_K_C (15) // Maximal conductance (mS/cm^2)
#define gbar_NMDA (0.014) //(0.0) // Maximal conductance (mS/cm^2)
#define gbar_AMPA (0.014) //(0.0) // Maximal conductance (mS/cm^2)
#define V_Na (120) // Reversal potential (mV)
#define V_Ca (140) // Reversal potential (mV)
#define V_K (-15) // Reversal potential (mV)
#define V_L (0) // Reversal potential (mV)
#define V_Syn (60) // Reversal potential (mV)
#define p (0.5) // ratio of the compartment length 
#define g_c (2.1) // coupling conductance (mS/cm^2) p43 after a detailed discussion
#define Vmref (-60) // Reference potential (mV)
#define S_max (125)

enum { SOMA, DEND }; // S: soma; D: dendrite

static double xi(const double x) { return (x/250.0) < 1 ? (x/250.0) : 1; }
static double H(const double x) { return (x >= 0) ? 1 : 0; }

static double V[NN][NS], dV[NN][NS]; // membrane potential (mV)
static double Ca[NN], dCa[NN]; // Intracellular free calcium level 
static double n[NN], h[NN], s[NN], c[NN], q[NN];
static double dn[NN], dh[NN], ds[NN], dc[NN], dq[NN];
static double S[NN], dS[NN], W[NN], dW[NN]; // Synaptic inputs for NMDA and AMPA
static const double I_s = -0.5; // Applied current to the soma (microA/cm^2)
static double I_d; // Applied current to the dendrite (microA/cm^2)

double alpha_m(const double v)
{
  return (0.32*(13.1-v))/(exp((13.1-v)/4)-1);
}
double beta_m(const double v)
{
  return (0.28*(v-40.1))/(exp((v-40.1)/5)-1);
}
double alpha_n(const double v)
{
  return (0.016*(35.1-v))/(exp((35.1-v)/5)-1);
}
double beta_n(const double v)
{
  return 0.25*exp(0.5-0.025*v);
}
double alpha_h(const double v)
{
  return 0.128*exp((17-v)/18);
}
double beta_h(const double v)
{
  return (4)/(1+exp((40-v)/5));
}
double alpha_s(const double v)
{
  return (1.6)/(1+exp(-0.072*(v-65)));
}
double beta_s(const double v)
{
  return (0.02*(v-51.1))/(exp((v-51.1)/5)-1);
}
double alpha_c(const double v)
{
  if (v <= 50){
    return (exp((v-10)/11)-exp((v-6.5)/27))/(18.975);
  }else{
    return 2*exp((6.5-v)/27);
  }
}
double beta_c(const double v)
{
  if (v <= 50){
    return 2*exp((6.5-v)/27)-alpha_c(v);
  }else{
    return 0;
  }
}
double alpha_q(const double ca)
{
  return (0.00002*ca < 0.01) ? 0.00002*ca : 0.01;
}
double beta_q(const double ca)
{
  return 0.001;
}
void initialize(void)
{
  for(int i = 0; i < NN; i++){
    V[i][SOMA] = -4.6; // mV (p45 left col)
    V[i][DEND] = -4.5; // mV
    Ca[i] = 0.2;
    h[i] = 0.999;
    n[i] = 0.001;
    s[i] = 0.009;
    c[i] = 0.007;
    q[i] = 0.010;
    S[i] = 0;
    W[i] = 0;
  }
}
void finalize(void)
{
  return;
}
double m_inf(const double v)
{
  return alpha_m(v)/(alpha_m(v)+beta_m(v));
}
double tau_m(const double v)
{
  return 1.0/(alpha_m(v)+beta_m(v));
}
double n_inf(const double v)
{
  return alpha_n(v)/(alpha_n(v)+beta_n(v));
}
double tau_n(const double v)
{
  return 1.0/(alpha_n(v)+beta_n(v));
}
double h_inf(const double v)
{
  return alpha_h(v)/(alpha_h(v)+beta_h(v));
}
double tau_h(const double v)
{
  return 1.0/(alpha_h(v)+beta_h(v));
}
double s_inf(const double v)
{
  return alpha_s(v)/(alpha_s(v)+beta_s(v));
}
double tau_s(const double v)
{
  return 1.0/(alpha_s(v)+beta_s(v));
}
double c_inf(const double v)
{
  return alpha_c(v)/(alpha_c(v)+beta_c(v));
}
double tau_c(const double v)
{
  return 1.0/(alpha_c(v)+beta_c(v));
}
double q_inf(const double ca)
{
  return alpha_q(ca)/(alpha_q(ca)+beta_q(ca));
}
double tau_q(const double ca)
{
  return 1.0/(alpha_q(ca)+beta_q(ca));
}

double I_Leak(const double v)
{
  return gbar_L*(v - V_L);
}
double I_Na(const double v, const double h)
{
  double m = m_inf(v);
  return gbar_Na*m*m*h*(v - V_Na);
}
double I_K_DR(const double v, const double n)
{
  return gbar_K_DR*n*(v - V_K);
}
double I_Ca(const double v, const double s)
{
  return gbar_Ca*s*s*(v - V_Ca);
}
double I_K_C(const double v, const double c, const double ca)
{
  return gbar_K_C*c*xi(ca)*(v - V_K);
}
double I_K_AHP(const double v, const double q)
{
  return gbar_K_AHP*q*(v - V_K);
}
double I_NMDA(const double v, const double S)
{
  return gbar_NMDA*S*(v - V_Syn)/(1+0.28*exp(-0.062*(v-60)));
}
double I_AMPA(const double v, const double W)
{
  return gbar_AMPA*W*(v - V_Syn);
}
double I_Syn(const double v, const double S, const double W)
{
  return I_NMDA(v, S) + I_AMPA(v, W);
}
double dhdt(const double h, const double v_s)
{
  return (h_inf(v_s) - h)/tau_h(v_s);
}
double dndt(const double n, const double v_s)
{
  return (n_inf(v_s) - n)/tau_n(v_s);
}
double dsdt(const double s, const double v_d)
{
  return (s_inf(v_d) - s)/tau_s(v_d);
}
double dcdt(const double c, const double v_d)
{
  return (c_inf(v_d) - c)/tau_c(v_d);
}
double dqdt(const double q, const double ca)
{
  return (q_inf(ca) - q)/tau_q(ca);
}
double dCadt(const double ca, const double i_ca)
{
  return -0.13*i_ca - 0.075*ca;
}
double dSdt(const double S, double v[NN][NS]) // Couldn't set v as const
{
  double sum = 0;
  for(int j = 0; j < NN; j++){
    sum += H(v[j][SOMA] - (10-Vmref));
  }
  return sum - S/150;
}
double dWdt(const double W, double v[NN][NS]) // Couldn't set v as const
{
  double sum = 0;
  for(int j = 0; j < NN; j++){
    sum += H(v[j][SOMA] - (20-Vmref));
  }
  return sum - W/2;
}

void compute_dV(const double V[], const double h, const double n, const double s, const double q,
		const double c, const double ca, const double S, const double W, const double I_d,
		double dV[])
{
  dV[SOMA] = (DT/Cm)*(-I_Leak(V[SOMA])
		   -I_Na(V[SOMA],h)
		   -I_K_DR(V[SOMA],n)
		   +(g_c/p)*(V[DEND]-V[SOMA])
		   +I_s/p);
  dV[DEND] = (DT/Cm)*(-I_Leak(V[DEND])
		   -I_Ca(V[DEND],s)
		   -I_K_AHP(V[DEND],q)
		   -I_K_C(V[DEND],c,ca)
		   -I_Syn(V[DEND], S, W)/(1-p)
		   +(g_c/(1-p))*(V[SOMA]-V[DEND])
		   +I_d/(1-p));
}

void output(FILE *file, const double t, const double vs, const double vd)
{
  fprintf(file, "%f %f %f\n", t, vs+Vmref, vd+Vmref);
}

static struct timeval start, stop;

void timer_start(void)
{
  gettimeofday(&start, NULL);
}

double timer_elapsed(void)
{
  struct timeval elapsed;
  gettimeofday(&stop, NULL);
  timersub(&stop, &start, &elapsed);
  return (elapsed.tv_sec*1000.0 + elapsed.tv_usec/1000.0)/1000.0;
}

int main(void)
{
  FILE *file;

  file = fopen("out.dat", "w");
  initialize();

  timer_start();
  for(int nt = 0; nt < NT; nt++){
    double t = DT*nt;
    output(file, t, V[0][SOMA], V[0][DEND]);

    if (nt == (int)(50/DT)){
      I_d = 11000.0;
    }else{
      I_d = 0.0;
    }

    for(int i = 0; i < NN; i++){
      dh[i] = DT*dhdt(h[i], V[i][SOMA]);
      dn[i] = DT*dndt(n[i], V[i][SOMA]);
      ds[i] = DT*dsdt(s[i], V[i][DEND]);
      dc[i] = DT*dcdt(c[i], V[i][DEND]);
      dq[i] = DT*dqdt(q[i], Ca[i]);
      dCa[i] = DT*dCadt(Ca[i], I_Ca(V[i][DEND], s[i]));
      dS[i] = DT*dSdt(S[i], V);
      dW[i] = DT*dWdt(W[i], V);
      compute_dV(V[i], h[i], n[i], s[i], q[i], c[i], Ca[i], S[i], W[i], I_d, dV[i]);
    }

    for(int i = 0; i < NN; i++){
      h[i] += dh[i];
      n[i] += dn[i];
      s[i] += ds[i];
      c[i] += dc[i];
      q[i] += dq[i];
      Ca[i] += dCa[i];
      S[i] += dS[i];
      if (S[i] > S_max) { S[i] = S_max; }
      W[i] += dW[i];
      for(int j = 0; j < NS; j++){
	V[i][j] += dV[i][j];
      }
    }
  }

  double elapsedTime = timer_elapsed();
  fprintf(stderr, "Elapsed time = %f sec\n", elapsedTime);

  finalize();
  fclose(file);

  return 0;
}

