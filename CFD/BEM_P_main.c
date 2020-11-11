#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<time.h>
#include<complex.h>
#include<unistd.h>
#include"solver.h" //Header file of skyline_method
#include"BEM_P.h"  //Header file of BEM
#define GNUPLOT_PATH "C:/gnuplot/bin/gnuplot.exe" //PATH for gnuplot.exe

int main(void){
    FILE *fp;
    /*--Constant--*/
    double U0=1.0;
    char OutputFile_name[64]="Output.txt";
    char title[64]="Cm vs alpha";
    /*--Variables--*/
    double _Complex *z_foil;
    char foil_name[256],file_name[256],file_format[]=".txt";
    double alpha,amin,amax,da;
    int nodes;
    int n;

	/*--Input foil name--*/
	printf("Input foil name\n>>");
	fgets(foil_name,sizeof(foil_name),stdin);
	foil_name[strlen(foil_name)-1]='\0';
	sprintf(file_name,"%s%s",foil_name,file_format); //Concatinate foil_name & file_format
    printf("Input amin,amax,da\n>>");
    scanf("%lf,%lf,%lf",&amin,&amax,&da);

    Check_FileSize(file_name,&nodes);
    z_foil=(double _Complex *)malloc(sizeof(double _Complex)*(nodes+1)); //Allocate matrix
    Input_Airfoil(file_name,nodes,z_foil);

    if((fp=fopen(OutputFile_name,"w"))==NULL){exit(EXIT_FAILURE);}else{fclose(fp);} //Renew output file

    n=0;
    do{
        alpha=amin+da*n;
        Analize_foil(OutputFile_name,nodes,z_foil,U0,alpha);
        n++;
    }while(alpha<amax);

    plot_polar(title,OutputFile_name);

    return 0;
}