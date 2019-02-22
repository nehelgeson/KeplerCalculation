#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define clean_input() while (getchar() != '\n') { }

double O, Ok_0, c, Th, rfid, *z, *O_phi_z, *rr, h0, afin, dobs, B_obs, cHobs, *Cinv, *DM_obs, *H_obs, H0_av, s; 

void phiCDM(double* results, double* w, double t, double* zz);
void TwoDimDot(double* TwoD, double* OneD, int len, double* result);
double OneDimDot(double* first, double* second, int len);


int main()
{
	double test[3][3] = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
	double test2[3] = { 1, 2, 3 };
	double result[3];

	TwoDimDot(test, test2, 3, result);

	printf("%f %f %f", result[0], result[1], result[2]);

	return 0;
}

/*
Good to go.
*/
void phiCDM(double* results, double* w, double t, double* zz)
{
	double p, v, a, al, k, m, K;
	p = w[0];
	v = w[1];
	a = w[2];
	al = zz[0];
	k = zz[1];
	m = zz[2];
	K = zz[3];
	results[0] = v;
	//printf("Value: %f\n", ((((4.0 / (9 * a*a*a))) + (1.0 / 12)*(v*v + (k*m) / (pow(p, al))) - K / (a*a))));
	results[1] = -3 * v*sqrt((((4.0 / (9 * a*a*a))) + (1.0 / 12)*(v*v + (k*m) / (pow(p, al))) - K / (a*a))) + ((k*al*m) / 2) / (pow(p, al + 1));
	results[2] = sqrt(((4.0 / 9) / a) + (((a*a) / 12))*(v*v + (k*m)*(pow(p, -al))) - K);
	//printf("%f %f %f\n", results[0], results[1], results[2]);
}

/*
IN PROGRESS
*/
double rs(double H0, double O, double Th)
{
	double h = H0 / 100;
	double hSq = h * h;
	double OhSq = O * hSq;
	double Ob = 0.02221 / (hSq);
	double ObhSq = Ob * hSq;
	double Th_pow4 = pow(Th, -4);
	double b2 = 0.238*(pow(OhSq, 0.223));
	double b1 = 0.313*(pow(OhSq, -0.419))*(1 + (0.607*pow(OhSq, (0.674))));
	double zd = 1291 * ((pow(OhSq, 0.251)) / (1 + (0.659*pow((OhSq), (0.828)))))*(1 + (pow((b1*(ObhSq)), (b2))));
	double Rd = 31.5*(ObhSq)*(Th_pow4)*(1000 / zd);
	double zeq = 25000 * (OhSq)*(Th_pow4);
	double Req = 31.5*(ObhSq)*(Th_pow4)*(1000.0 / zeq);
	double keq = 0.0746*(OhSq)*(pow(Th, (-2)));
	double A = sqrt(1 + Rd);
	double B = sqrt(Rd + Req);
	double C = 1 + sqrt(Req);


	return (2.0 / (3 * keq))*(sqrt(6 / Req))*(log((A + B) / C));
}

double E(double O, double red, double Ok_0, double Ophiz)
{
	double temp = red + 1;
	return sqrt(O*(temp*temp*temp) + Ok_0*(temp*temp) + Ophiz);
}

double D_M(double H0, double Ok_0, double x, double c, double h0, double afin)
{
	double y;
	if (Ok_0 < 0)
	{
		y = (1 / (sqrt(-(Ok_0))))*(sin((sqrt(-(Ok_0)))*h0*afin*x));
	}
	if (Ok_0 > 0)
	{
		y = (1 / (sqrt((Ok_0))))*(sinh((sqrt((Ok_0)))*h0*afin*x));
	}
	if (Ok_0 == 0)
	{
		y = h0*afin*x;
	}
	return y * (c / H0);
}

double chi_sq(double H0, double O, double Ok_0, double c, double Th, double rfid, double* z, double* O_phi_z, double* rr, double h0, double afin, double dobs, double B_obs, double cHobs, double* Cinv, double* DM_obs, double* H_obs)
{
	double chi2B, chi2H, chi2d, chi2DV1, chi2DV2;
	chi2B = chi2H = chi2d = chi2DV1 = chi2DV2 = 0;
	double dH = c / H0;
	double r = rs(H0, O, Th) / rfid;
	double* DM_th = (double *)malloc(8 * sizeof(double));
	double* H_th = (double *)malloc(8 * sizeof(double));
	int q;

	for (q = 0; q < 8; q++)
	{
		double z1 = z[q];
		printf("Value of z1 for q value %d is %lf\n", q, z1);
		double H1 = H0*E(O, z1, Ok_0, O_phi_z[q]);
		double DM = D_M(H0, Ok_0, rr[q], c, h0, afin);
		double y = DM / dH;

		if (q <= 2)
		{
			DM_th[q] = DM;
			H_th[q] = H1;
		}
		if (z1 == 0.15) //weird floating point errors? Look here
		{
			double DV = dH*pow(((y*y)*z1) / (E(O, z1, Ok_0, O_phi_z[q])), (1.0/3));
			double DVobs1 = 664 * r;
			double unc1 = 25 * r;
			chi2DV1 = ((DV - DVobs1)*(DV - DVobs1)) / (unc1*unc1);
		}
		if (z1 == 0.106)
		{
			double DV = (c / H0)*pow(((y*y)*(z1)) / (E(O, z1, Ok_0, O_phi_z[q])), (1.0 / 3));
			double dth = (rs(H0, O, Th)) / DV;
			chi2d = ((dth - dobs)*(dth - dobs)) / (0.015*0.015);
		}
		if (z1 == 1.52)
		{
			double DV = dH*pow(((y*y)*z1) / (E(O, z1, Ok_0, O_phi_z[q])), (1.0 / 3));
			double DVobs2 = 3855 * r; //Ata
			double unc2 = 170 * r;
			chi2DV2 = ((DV - DVobs2)*(DV - DVobs2)) / (unc2*unc2);
		}
		if (z1 == 2.33)
		{
			double DH = c / H1;
			double F = pow(DH, (0.7));
			double G = pow(DM, (0.3));
			double B_th = (F*G) / (rs(H0, O, Th));
			chi2B = ((B_th - B_obs)*(B_th - B_obs)) / (0.35*0.35);
		}
		if (z1 == 2.36)
		{
			double cHth = (c / (rs(H0, O, Th)))*(1 / H1);
			chi2H = ((cHth - cHobs)*(cHth - cHobs)) / (0.3*0.3);
		}
	}

	double* result = (double *)malloc(6 * sizeof(double));
	double delta[6] = { DM_th[0] / r - DM_obs[0],
						r*H_th[0] - H_obs[0],
						DM_th[1] / r - DM_obs[1],
						r*H_th[1] - H_obs[1],
						DM_th[2] / r - DM_obs[2],
						r*H_th[2] - H_obs[2] };
	TwoDimDot(Cinv, delta, 6, result);
	double chi_sq_11 = OneDimDot(delta, result, 6);

	free(DM_th);
	free(H_th);
	free(result);
	
	printf("C: chi2B: %lf, chi2H: %lf, chi2d: %lf\nC: chi_sq_11: %lf, chi2DV1: %lf, chi2DV2: %lf\n", chi2B, chi2H, chi2d, chi_sq_11, chi2DV1, chi2DV2);
	return chi2B + chi2H + chi2d + chi_sq_11 + chi2DV1 + chi2DV2;
}

double OneDimDot(double* first, double* second, int len)
{
	int i;
	double total = 0.0;
	for (i = 0; i < len; i++) 
	{
		total += first[i] * second[i];
	}
	return total;
}

void TwoDimDot(double* TwoD, double* OneD, int len, double* result)
{
	int i, j;
	for (i = 0; i < len; i++)
	{
		double total = 0.0;
		for (j = 0; j < len; j++)
		{
			total += OneD[j] * TwoD[i*len+j];
		}
		result[i] = total;
	}
}

double IntegrateFunc(double H0)
{
	printf("Inside Integrate Func\n");
	double chi_sq_result = chi_sq(H0, O, Ok_0, c, Th, rfid, z, O_phi_z, rr, h0, afin, dobs, B_obs, cHobs, Cinv, DM_obs, H_obs);
	printf("Chi Sq result: %lf\n", chi_sq_result);
	double final_result = exp((-0.5 * chi_sq_result) + (-1.0 / (2*s*s)) * (H0 - H0_av)*(H0 - H0_av));
	printf("Final result: %lf\n", final_result);
	return final_result;
}

void SetGlobals(double pO, double pOk_0, double pc, double pTh, double prfid, double* pz, double* pO_phi_z, double* prr, double ph0, 
				double pafin, double pdobs, double pB_obs, double pcHobs, double* pCinv, double* pDM_obs, double* pH_obs, double pH0_av, 
				double ps, int pz_len, int pO_phi_z_len, int prr_len, int pCinv_len, int pDM_obs_len, int pH_obs_len)
{
	int i; 
	for (i=0; i < 8; i++)
	{
		printf("Value of z at %d is %lf\n", i, pz[i]);
	}
	O = pO;
	Ok_0 = pOk_0;
	c = pc;
	Th = pTh;
	rfid = prfid;
	h0 = ph0;
	afin = pafin;
	dobs = pdobs;
	B_obs = pB_obs;
	cHobs = pcHobs;
	H0_av = pH0_av;
	s = ps;
	
	z = (double*)malloc(pz_len*sizeof(double));
	O_phi_z = (double*)malloc(pO_phi_z_len*sizeof(double));
	rr = (double*)malloc(prr_len*sizeof(double));
	Cinv = (double*)malloc(pCinv_len*sizeof(double));
	DM_obs = (double*)malloc(pDM_obs_len*sizeof(double));
	H_obs = (double*)malloc(pH_obs_len*sizeof(double));
	
	memcpy(z, pz, pz_len*sizeof(double));
	memcpy(O_phi_z, pO_phi_z, pO_phi_z_len*sizeof(double));
	memcpy(rr, prr, prr_len*sizeof(double));
	memcpy(Cinv, pCinv, pCinv_len*sizeof(double));
	memcpy(DM_obs, pDM_obs, pDM_obs_len*sizeof(double));
	memcpy(H_obs, pH_obs, pH_obs_len*sizeof(double));
}


