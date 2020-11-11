#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include"solver.h"
#define GNUPLOT_PATH "C:/gnuplot/bin/gnuplot.exe" //PATH for gnuplot.exe

void plot_results(int nmax,int nmin[nmax+1]);

int main(void){
    clock_t s1,s2,e1,e2;
    int nmax=32,amax=10;
    double b[nmax+1],b2[nmax+1],x[nmax+1];
    double *L,*L2,**A,*Ap,*b1;
    int nd[nmax+1],nmin[nmax+1]; //n_diagonal,n_min
    int i,j,k,j0,ssize; //skyline_size
    int ij,ji,n;

    A=(double **)malloc(sizeof(double *)*(nmax+1));Ap=(double *)malloc(sizeof(double)*(nmax+1)*(nmax+1));for(n=0;n<=nmax;n++){A[n]=Ap+(nmax+1)*n;}
    b1=(double *)malloc(sizeof(double)*(nmax+1));

    for(i=0;i<=nmax;i++){for(j=0;j<=nmax;j++){A[i][j]=0;}} //Initialize A

    //Create random matrix
    ssize=0;
    srand((int)time(NULL));
    for(i=0;i<=nmax;i++){
        j0=rand() % (i+1); 
        for(j=j0;j<=i;j++){
            A[i][j]=rand() % (amax+1)+1;
            ssize++;
        }
        b[i]=rand() % (amax+1);
    }
    ssize--;
    //Make A symmetry
    for(i=0;i<=nmax;i++){
        for(j=0;j<i;j++){
            A[j][i]=A[i][j];
        }
    }
    L=malloc(sizeof(double)*(ssize+1)); //Allocate L
    L2=malloc(sizeof(double)*(ssize+1)); //Allocate L2

    //Convert A to L
    k=0;
    for(i=0;i<=nmax;i++){
        for(j=i;j>=0;j--){
            nmin[i]=0;
            if(A[i][j]==0){
                //Make nmin
                nmin[i]=j+1;
                break;
            }else{
                L[k]=A[i][j];
            }
            if(i==j){nd[i]=k;}
            k++;
        }
    }
    for(i=0;i<=ssize;i++){L2[i]=L[i];} //Copy L to L2

    //Show A, nd and nmin
    //printf("\n");for(i=0;i<=nmax;i++){for(j=0;j<=nmax;j++){printf("%2.0f ",A[i][j]);}printf(" nd[%2d]=%3d nmin[%2d]=%2d\n",i,nd[i],i,nmin[i]);}printf("\n");

    for(i=0;i<=nmax;i++){b1[i]=b[i];b2[i]=b[i];} //Copy b to b1, b2
    s1=clock();
    Gaussian_Elimination(nmax,A,b1);
    e1=clock();s2=clock();
    skyline_method(nmax,ssize,nd,nmin,L,b2);
    e2=clock();

    for(j=0;j<=nmax;j++){printf("b1[%3d]=%7.3f b2[%3d]=%7.3f\n",j,b1[j],j,b2[j]);}
    
    //Show reslut
    printf("\nGaussian Elimination\nMatrix size >> %d\nTime >> %f\n",(nmax+1)*(nmax+1),(e1-s1)/(double)CLOCKS_PER_SEC);
    printf("\nSkyline method\nMatrix size >> %d\nTime >> %f\n",ssize,(e2-s2)/(double)CLOCKS_PER_SEC);

    plot_results(nmax,nmin);

    free(L);
    return 0;
}

void plot_results(int nmax,int nmin[nmax+1]){
    FILE *gp; //gp:gnuplot pointer
    double xy_contour[2][2],z0,zmin=0,zmax=0;
    char plot1[64],plot2[64],plot3[64],title[32]="Skyline method";
    int i,j,k,n,kmax=30;
    int *ip,*jp;

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){
        fprintf(stderr,"Not Found %s.",GNUPLOT_PATH);
        exit(EXIT_FAILURE);
    }

    //Send comands to gnuplot
    fprintf(gp,"set title '%s'\n",title);
    fprintf(gp,"set xrange [%f:%f]\n",0,(double)(nmax+1));
    fprintf(gp,"set yrange [%f:%f]\n",0,(double)(nmax+1));
    fprintf(gp,"unset xtics\n");fprintf(gp,"unset ytics\n");
    fprintf(gp,"set style fill solid\n");
    fprintf(gp,"unset key\n"); //Hide legend
    fprintf(gp,"set size square\n");

    //Plot graphic
    strcpy(plot1,"'-' with boxes linetype rgbcolor 'red'");
    strcpy(plot2,"'-' with boxes linetype rgbcolor 'white'");
    fprintf(gp, "plot %s,%s\n",plot1,plot2);

    //Plot Skyline
    for(k=0;k<=nmax;k++){
        fprintf(gp,"%f\t%f\n",(double)k+0.5,(double)((nmax+1)-nmin[k]));
        fprintf(gp,"\n");
    }
    fprintf(gp,"e\n"); //End of array
    for(k=0;k<=nmax;k++){
        fprintf(gp,"%f\t%f\n",(double)k+0.5,(double)((nmax+1)-(k+1)));
        fprintf(gp,"\n");
    }
    fprintf(gp,"e\n"); //End of array

    fflush(gp);
    printf("\n");
    system("pause");
    fprintf(gp, "exit\n");
    _pclose(gp);

    return;
}