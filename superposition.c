#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ein.h"
#define ARRSIZE 5000

void tabulargenerator(int max_secs){
	int fr, dur, totdur=0;
	FILE *tabulargenerated;
	tabulargenerated = fopen("tabulargenerated.csv", "w+");
	while(1){
		fr = (rand() % 250) * 5;
		dur = (rand() % 10) * 86400; //to have a minimum flow period of 1 day
		fprintf(tabulargenerated,"%d,%d,%d\n",totdur,totdur+dur,fr);
		totdur += dur;
		if(totdur>max_secs){break;}
	}
	fclose(tabulargenerated);
}

void superposition(int end[],int t_size, int flow[],double pi,double B,double mu,double k,double h,double phi,double ct,double r,int max_dur){
	FILE *pressure_storage;
	pressure_storage = fopen("pressure.csv", "w+");
	fprintf(pressure_storage,"Time(s),Pressure(psi)\n");
	int t_idx;
	int t;
	for(t=0;t<max_dur;++t){
		int idx;
		for(idx=0;idx<t_size;++idx){
			if(end[idx]>=t){
				t_idx = idx;
				break;
			}
			t_idx = t_size;
		}
		double pressure_drop = 0;
		int i;
		for(i=0;i<t_idx;++i){
			pressure_drop += (70.6*(flow[i+1]-flow[i])*B*mu/(k*h)) * -Exponential_Integral_Ei(-(948*phi*mu*ct*pow(r,2))/(k*(t-end[i])/3600));
		}
	fprintf(pressure_storage, "%d,%f\n",t,pi-pressure_drop);
	}
	fclose(pressure_storage);
}

int main(){
	srand(time(NULL));
	FILE *tabular;
	//double pi, B, mu, k, h, phi, ct, r;
	int i, var1, var2, var3, max_dur, decision;
	int start_arr[ARRSIZE], end_arr[ARRSIZE], flow_arr[ARRSIZE];
	start_arr[0] = 0;
	end_arr[0] = 0;
	flow_arr[0] = 0;

	printf("Pressure response generator program v1.0\n");
	printf("This program will create 1Hz pressure response data for IARF.\n");
	printf("Please note that IARF behavior is not realistic due to the finite radius. ");
	printf("The aim of this program is to generate a data that ");
	printf("has the characteristics of a pressure response data for various applications.\n");
	printf("Please select tabular data:\n");
	printf("1-Use Example case\n");
	printf("2-Generate Random data\n");
	scanf("%d",&decision);
	printf("\n");
	printf("Please enter dataset max. duration in seconds:\t");
	scanf("%d",&max_dur);

	if(decision==1){tabular = fopen("tabularexample.csv", "r+");}
	else if(decision==2){
		tabulargenerator(max_dur);
		tabular = fopen("tabulargenerated.csv", "r");
	}
	else{
		printf("Please select by typing an integer. Program will be terminated.\n");
		return 0;
	}

	for(i=1; i<ARRSIZE-1; ++i){
		fscanf(tabular,"%d,%d,%d",&var1,&var2,&var3);
		start_arr[i] = var1;
		end_arr[i] = var2;
		flow_arr[i] = var3;
		if(i==ARRSIZE-2){flow_arr[i+1]=var3;} //appending last flow rate twice (see explanation)
	}
	fclose(tabular);

	int t_size = sizeof end_arr / sizeof *end_arr; // get the length (size) of time list
	superposition(end_arr,t_size,flow_arr,2906,1.27,0.1,20,30,0.28,0.000003,100,max_dur);

	return 0;
}
