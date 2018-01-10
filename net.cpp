#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
//include's de c++
#include <string>
#include <sstream>
using namespace std;

template<class T>
inline string to_string(const T& t){
    stringstream ss;
    ss << t;
    return ss.str();
}

#define N 100
#define TMAX 500
#define dt 0.01
#define PI 3.14159
#define w0 2
#define DELTAMAX 0.5

double** allocMatrix(int n,int m){

	double** matrix=(double**)malloc(n * sizeof(double*));
	int i;
	for(i=0;i<m;i++){
	    matrix[i] = (double*)malloc(m * sizeof(double));		
	}
	return matrix;
} 
void freeMatrix(double** matrix,int m){
	int i;
	for(i=0;i<m;i++){
	    free(matrix[i]);		
	}
	free(matrix);
} 

void RK2(double *v,double *w,double **matrix,double epsilon1,double epsilon2,double phi){

    double *vaux;
    double *K1,*K2;
    int i;

    vaux = (double*)malloc(N*sizeof(double));
    K1 = (double*)malloc(N*sizeof(double));
    K2 = (double*)malloc(N*sizeof(double));

    for(i=0;i<N;i++)
       	vaux[i]=v[i];
    K1=Coupling(vaux,w,matrix,epsilon1,epsilon2,phi);
   
    for(i=0;i<N;i++)
       	vaux[i]=v[i]+dt*K1[i];
    K2=Coupling(vaux,w,matrix,epsilon1,epsilon2,phi);
    
    for(i=0;i<N;i++){
       	v[i]=v[i]+0.5*(dt*K1[i]+dt*K2[i]);
	if (v[i]>2*PI){
		v[i]=v[i]-2*PI;
		}
	}

    free(vaux);
    free(K1);
    free(K2);
}

void ConnectivityAlltoAll(double** matrix){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			matrix[i][j]=1;
			if(i==j)
				matrix[i][j]=0;		
		}
	}
}


double* Coupling(double* v,double* w,double** matrix,double epsilon1,double epsilon2,double phi){

	double* x;
    	x = (double*)malloc(N*sizeof(double));
	int i,j;

	for(i=0;i<N;i++){
		x[i]=w[i];
		for(j=0;j<N;j++){
			x[i]+=matrix[i][j]*InternalCoupling(epsilon1,v[j],v[i]);		
		}
		x[i]+=ExternalCoupling(epsilon2,phi,v[i]);  
	}	
	
	return x;
}

void InitialConditions(double *v){

	int i;
	for(i=0;i<N;i++){
		//v[i]=2*PI*(1.0*rand()/RAND_MAX);
		v[i]=(0.1*rand()/RAND_MAX);
	}
}

void PrintMatrix(double** matrix,int n,int m){

	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<m;j++)
			printf("%lf\t",matrix[i][j]);
		printf("\n");
	}
}

int main(){

	string directorio1 = to_string("/home/melisa/Desktop/HippocampusProject/Results/");
	FILE *file1,*file2,*file3,*file0,*file4;
	string output0=directorio1+"synchr_vs_delta_gauss.dat";
	if((file0 = fopen(output0.c_str(), "w")) == NULL){
		printf("No puedo abrir el archivo %s.\n", output0.c_str());
		exit(1);
	}
	/*string output4=directorio1+"phase_vs_delta_2osc.dat";
	if((file4 = fopen(output4.c_str(), "w")) == NULL){
		printf("No puedo abrir el archivo %s.\n", output4.c_str());
		exit(1);
	}*/

	srand(time(NULL));
	double* neuron=(double*)malloc(N*sizeof(double));
	double** conn=allocMatrix(N,N);
	double* freq=(double*)malloc(N*sizeof(double));
	double* max=(double*)malloc(N*sizeof(double));
	int* counter=(int*)malloc(N*sizeof(int));
	double* freq_final=(double*)malloc(N*sizeof(double));
	double eps1=1,eps2=0;
	double phase=0;
	double t;
	int i,j;

	for (i=0;i<N;i++)
		freq[i]=0;

	Connectivity_1(conn);
	//ConnectivityAlltoAll(conn);
	//PrintMatrix(conn,N,N);
	//OscillatorsGaussianProfile(freq);
	InitialConditions(neuron);


	double delta=0;
	while(delta<DELTAMAX){
		//printf("%lf\n",delta);
		//OscillatorsUniformProfile(freq,w0-delta,w0+delta);
		//OscillatorsFreqProfile(freq,w0-delta,w0+delta);//here we put the linear variation along longitudinal axes
		OscillatorsGaussianProfile(freq,w0-delta,w0+delta,1);//here we add gaussian noise
		/*string output1=directorio1+"oscillators"+"_"+to_string(delta)+".dat";
		string output2=directorio1+"synchronization"+"_"+to_string(delta)+".dat";
		string output3=directorio1+"frequency"+"_"+to_string(delta)+".dat";
		if((file1 = fopen(output1.c_str(), "w")) == NULL){
			printf("No puedo abrir el archivo %s.\n", output1.c_str());
			exit(1);
		}
		if((file2 = fopen(output2.c_str(), "w")) == NULL){
			printf("No puedo abrir el archivo %s.\n", output2.c_str());
			exit(1);
		}
		if((file3 = fopen(output3.c_str(), "w")) == NULL){
			printf("No puedo abrir el archivo %s.\n", output3.c_str());
			exit(1);
		}
		for(i=0;i<N;i++){
			max[i]=0;
			counter[i]=0;
		}*/
		t=0;
		double sinchr=0;	
		while(t<TMAX){
			RK2(neuron,freq,conn,eps1,eps2,phase);
			//fprintf(file1, "%lf\t",t);
			double x=0,y=0;
			for(i=0;i<N;i++){
				//printf(file1, "%lf\t",neuron[i]);
				x+=cos(neuron[i]);
				y+=sin(neuron[i]);
				/*if(neuron[i]>max[i])
					max[i]=neuron[i];
				if(neuron[i]<max[i]){
					counter[i]+=1;
					max[i]=neuron[i];				
					} */
			}
			//fprintf(file1, "\n");
			sinchr=1.0*sqrt(x*x+y*y)/N;
			//fprintf(file2,"%lf\t%lf\n",t,sinchr);
			t+=dt;
		}

		/*for(i=0;i<N;i++){
			freq_final[i]=2.0*PI*counter[i]/TMAX;			
			fprintf(file3, "%lf\n",freq_final[i]);
			}*/

		//fclose(file1);
		//fclose(file2);
		//fclose(file3);
		fprintf(file0,"%lf\t%lf\n",delta,sinchr);
		printf("%lf\t%lf\n",delta,sinchr);
		//printf("%lf\t%lf\n",2.0*delta,sinchr);
		/*fprintf(file4, "%lf\t",delta);
		for(i=0;i<N;i++){
			fprintf(file4, "%lf\t",neuron[i]);
			}
		fprintf(file4,"\n");*/
		delta+=0.05;
		
	}

	free(max);
	free(counter);
	free(neuron);
	free(freq_final);
	freeMatrix(conn,N);
	free(freq);
	fclose(file0);
	return 0;
}
