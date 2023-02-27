#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define L 80.0			// Longitude [m]
#define D 0.5			// Diameter [m]
#define c 1.0		// wave speed [m/s]
#define simTime 10.0	// Simulation time [s]
#define nCells 100		    // Number of nodes
//#define dx (L/(imax-1))     // Spatial increment
#define CFL 1.0		// CFL number
#define g 9.81

#define INITIAL_CONDITION 1
#define MESH_REGULARITY 1

#define v0 5.0			// Initial velocity
#define PI (4*atan(1))

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define NormRANu (2.3283063671E-10F)
#define Frec_Max 100
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
int frec[Frec_Max];
extern float Random(void);
extern void ini_ran(int SEMILLA);

int PrintToFile (int *PrintCounter, double x[], double v[], double t);
double GlobalTimeStep (double dx[], double t, double *dtLocal);
double FluxFunction (int i, double v[], double vL[], int *FluxCounter);
double SquareWave (double x);
double Gaussian (double x);
void ErrorNorms (double v[], double dx[], double x[]);

int main () {

/*

FALTA INCORPORAR GESTIÓN DE FRONTERAS DE REGIONES DE DIFERENTE M-VALUE.
AHORA REPRODUCE LA INFORMACIÓN PERO NO ESTIMA BIEN LA VELOCIDAD DE PRO-
PAGACIÓN DE LA ONDA CUADRADA. ¡¡¡CHECK!!!

*/

    srand(time(NULL));
//    ini_ran(rand());

    int semilla = 32207;
    ini_ran(semilla);

    int i,iter=0,k=0,printFreq=1;
	double v[nCells];		// Velocity arrays (conserved variables)
	double vL[nCells];      // Interpolation points
	double Q[nCells], x[nCells], dx[nCells];			// Discharge, distance
	double t=0.0, dt, dtGlobal, dtLocal[nCells], tLocal=0., tFrozen[nCells];                // Time
	double flux[nCells];
	double A;                       // Area
	int PrintCounter=0, FluxCounter=0, IterCounter=0;             // Print counter
	int mLocal[nCells], mMax=0, mPow2Local[nCells];
	int p=0, pTotal=0;
//	bool allowsFrozenFlux[nCells];
	double log2 = log(2.);

    printf("\n   Linear convection");
	printf("\n   --------------------------------------------------");


// Reading simulation parameters

	printf("\n>> Simulation setup loaded");
	printf("\n   Simulation time: %.5lf s",simTime);
	printf("\n   CFL: %.2lf",CFL);
	printf("\n   Number of cells: %d",nCells);
	printf("\n   dx: %.2lf m",dx);
	printf("\n   Pipe length: %.2lf m\n",L);


	// Variable initialization
	for (i=0;i<nCells;i++){
        v[i] = 0.0;
        vL[i] = 0.0;
        Q[i] = 0.0;
        x[i] = 0.0;
        flux[i] = 0.0;
        tFrozen[i] = 0.0;
	}

//	Irregular mesh

#if MESH_REGULARITY == 1

	dx[0]=1.0;
	x[0]=dx[0]/2;
	for (i=1;i<35;i++){
        dx[i] = 1.0;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (i=35;i<40;i++){
        dx[i] = 0.5;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (i=40;i<60;i++){
        dx[i] = 0.25;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (i=60;i<65;i++){
        dx[i] = 0.5;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (i=65;i<nCells;i++){
        dx[i] = 1.0;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}

#elif MESH_REGULARITY == 2

    double to100m = 100./nCells;
    dx[0] = 2.*to100m*Random();
    x[0] = dx[0]/2.;
    for (i=0; i<nCells; i++){
        dx[i] = 2.*to100m*Random();
        x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i]);
    }

#else

// Undefined

#endif // MESH_REGULARITY

// INITIAL CONDITION SETUP

#if INITIAL_CONDITION == 1

        for (i=0; i<nCells; i++) v[i] = SquareWave(x[i]);

#else

        for (i=0; i<nCells; i++) v[i] = Gaussian(x[i]);

#endif // INITIAL_CONDITION

//  Testing mesh initialization
//    for (i=0;i<nCells;i++){
//        printf("%d\t%.3f\t%.3f\n",i,dx[i],x[i]);
//    }

//  Testing velocity initialization
//    for (i=0;i<nCells;i++){
//        printf("%d\t%.3f\n",i,v[i]);
//    }
//  getchar();

    PrintToFile(&PrintCounter, x, v, t);

    while (t < simTime){

//      Initiating Global Time Step

//        printf("Iniciando Global Time Step...\n");

        dt = GlobalTimeStep(dx,t,dtLocal);

        if (t+dt > simTime){
            dt = simTime - t;
        }

//        printf("dt=%lf\n",dt);
//        getchar();

        mLocal[0] = (int) log(dtLocal[0]/dt)/log2;
        mMax = mLocal[0];
        for (i=1; i<nCells; i++){
            mLocal[i] = (int) log(dtLocal[i]/dt)/log2;
            mPow2Local[i] = (int) pow(2,mLocal[i]);
            mMax = MAX(mMax,mLocal[i]);
        }
        pTotal = (int) pow(2,mMax);

//        printf("dtGlobal=%lf\n",dtGlobal);
//        printf("kMax=%d\n",kMax);
//        getchar();

//        if (t+dtGlobal > simTime){
//            dtGlobal =  simTime - t;
//        }

        for (i=0; i<nCells; i++){
//            allowsFrozenFlux[i] = (dtLocal[i] >= 2.*dt);
//            printf(allowsFrozenFlux[i] ? "allowsFrozenFlux[i]=true\n" : "allowsFrozenFlux[i]=false\n");

            vL[i] = v[i];
//            printf("dtLocal[i]=%lf flux[i]=%lf\n",dtLocal[i],flux[i]);
        }

//        getchar();

//        tLocal = t + dtGlobal;

//        while (t < tLocal ){
////            printf("dt=%lf\n",dt);
//            t += dt;
////            printf("tLocal=%lf\n", tLocal);
////            getchar();
//
//            for (i=0; i<nCells; i++){
//                v[i] = v[i] - dt/dx[i]*flux[i];
//            }
//
//            dt = GlobalTimeStep(dx,tLocal,dtLocal);
//            if(t + dt > tLocal){
//                dt = tLocal - t;
//            }
//
//            for (i=0; i<nCells; i++){
//                if (allowsFrozenFlux[i]){
//
////                    printf("i=%d\n",i);
//
//                    if (t + dt > tFrozen[i]){
//                        flux[i] = FluxFunction(i, v, &FluxCounter);
//                        allowsFrozenFlux[i] = (dtLocal[i] >= 2.*dt);
//                        tFrozen[i] = t + dtLocal[i];
////                        printf("i=%d",i);
//
////                        printf("i=%d, t=%lf\n",i,t);
//
//                    }
//
//                } else {
//                    //allowsFrozenFlux[i] = (dtLocal[i] > 2.*dt);
//                    flux[i] = FluxFunction(i, v, &FluxCounter);
//                }
//            }

            for (p=0; p<pTotal; p++){

                for (i=0; i<nCells; i++){

                    if (p % mPow2Local[i] == 0){

                        vL[i] = v[i];
                        v[i] = v[i] - dt * mPow2Local[i] / dx[i] * FluxFunction(i, v, vL, &FluxCounter);

                    }

                }

                t = t + dt;

            }

            IterCounter++;

//            getchar();

        }


//    t = t + dtGlobal;



    PrintToFile(&PrintCounter, x, v, t);
    printf("Flux function evaluations:%d, iterations=%d",FluxCounter,IterCounter);
    ErrorNorms(v,dx,x);
        return 0;
}

int PrintToFile(int *PrintCounter, double x[], double v[], double t){

    int i;
    char filename[40];

    FILE *output;

//	Write the numbering in the file name with the corresponding extension
    sprintf(filename,"OutputFiles/conveccionLinealLocalTimeStep%d.dat",*PrintCounter);
    output=fopen(filename,"w");

    fprintf(output,"#t=%lf\n",t);
    for(i=0;i<nCells;i++){
       fprintf(output,"%f %lf\n",x[i],v[i]);
    }

    fclose(output);
    *PrintCounter=*PrintCounter+1;

    return 0;

}

double GlobalTimeStep(double dx[], double t, double *dtLocal){

    int i = 0;
    double dt = 0.0;

    dtLocal[0] = dx[0]/c;
    dt = dtLocal[0];

    for (i = 1; i < nCells; i++){
        dtLocal[i] = dx[i]/c;
        dt = MIN(dt,dtLocal[i]);
    }

    return dt;
}

double FluxFunction(int i, double v[], double vL[], int *FluxCounter){

    *FluxCounter = *FluxCounter + 1;

    if (i==0){
        return 0.;
    }

    return c * (v[i]- 0.5 * (v[i-1] + vL[i-1]));

}

double SquareWave (double x){

    if (x < 25. || x > 55.) return 0.;

    return v0;
}

double Gaussian (double x){

    double sigma = 5.0;

    return v0*exp(-(x-40.)*(x-40.)/2./sigma/sigma);

}

void ErrorNorms (double v[], double dx[], double x[]){

    double L1=0., L2=0., Linf=0.;
    int i;

#if INITIAL_CONDITION == 1

    for (i=0; i<nCells; i++){
        L1 = L1 + fabs(v[i]-SquareWave(x[i]-c*simTime))*dx[i];
        L2 = L2 + pow(fabs(v[i]-SquareWave(x[i]-c*simTime)),2)*dx[i];
        Linf = MAX(Linf, fabs(v[i]-SquareWave(x[i]-c*simTime)));
    }

#else

    for (i=0; i<nCells; i++){
            L1 = L1 + fabs(v[i]-Gaussian(x[i]-c*simTime))*dx[i];
            L2 = L2 + pow(fabs(v[i]-Gaussian(x[i]-c*simTime)),2)*dx[i];
            Linf = MAX(Linf, fabs(v[i]-Gaussian(x[i]-c*simTime)));
    }

#endif

    printf("\n\nError norms:\n");
    printf("\tL1=%lf",L1);
    printf("\tL2=%lf",sqrt(L2));
    printf("\tLinf=%lf",Linf);

}


float Random(void)
{
float r;
ig1=ind_ran-24;
ig2=ind_ran-55;
ig3=ind_ran-61;
irr[ind_ran]=irr[ig1]+irr[ig2];
ir1=(irr[ind_ran]^irr[ig3]);
ind_ran++;
r=ir1*NormRANu;
//printf("r=%f\n",r);
return r;
}

void ini_ran(int SEMILLA)
{
printf("%d",SEMILLA);
int INI,FACTOR,SUM,i;
srand(SEMILLA);
INI=SEMILLA;
FACTOR=67397;
SUM=7364893;
for(i=0;i<256;i++)
{
INI=(INI*FACTOR+SUM);
irr[i]=INI;
}
ind_ran=ig1=ig2=ig3=0;
}









