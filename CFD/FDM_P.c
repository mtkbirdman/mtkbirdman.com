#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#define GNUPLOT_PATH "C:/gnuplot/bin/gnuplot.exe" //PATH for gnuplot.exe

void contour(int imax,int jmax,double,double z[imax+1][jmax+1],double xy[imax+1][jmax+1][2],double cntr[][2],int *,int *);
void plot_results(int imax,int jmax,int,int,int,double,double,double,double,double z[imax+1][jmax+1],double xy[imax+1][jmax+1][2]);

int main(void){
    int imax=200,jmax=100,i1=80,i2=120,jbox=40,n_max=10000;
    double xmin=0,xmax=2,ymin=0,ymax=1;
    double xy[imax+1][jmax+1][2];
    double psi_max=1.0,eps=0.00001,e_max,psi[imax+1][jmax+1],ds,s,dpsi;
    int i,j,k,n;

    //Create mesh
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            xy[i][j][0]=(xmax-xmin)*((double)i/(double)imax);
            xy[i][j][1]=(ymax-ymin)*((double)j/(double)jmax);
        }
    }

    //Initialize stream function
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            psi[i][j]=0;
        }
    }

    //Input boundry conditions
    for(i=0;i<=imax;i++){ psi[i][0]=0; } //bottom wall
    for(j=0;j<=jbox;j++){ psi[i1][j]=0; } //box left-side wall
    for(i=i1;i<=i2;i++){ psi[i][jbox]=0; } //box top wall
    for(j=0;j<=jbox;j++){ psi[i2][j]=0; } //box right-side wall
    for(i=0;i<=imax;i++){ psi[i][jmax]=psi_max; } //free stream
    ds=psi_max/jmax;
    for(j=0;j<=jmax;j++){
        s=ds*j;
        psi[0][j]=s; //left-side inflow
        psi[imax][j]=s; //right-side outflow
    }

    //Gauss-Seidel method
    n=0;
    do{
        e_max=eps;
        for(i=1;i<imax;i++){
            if(i<i1||i>i2){
                for(j=1;j<jmax;j++){
                    dpsi=(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1])/4-psi[i][j];
                    psi[i][j]=psi[i][j]+dpsi;
                    if(fabs(dpsi)>e_max){e_max=fabs(dpsi);} //update maximum error
                }    
            }else{
                for(j=jbox+1;j<jmax;j++){
                    dpsi=(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1])/4-psi[i][j];
                    psi[i][j]=psi[i][j]+dpsi;
                    if(fabs(dpsi)>e_max){e_max=fabs(dpsi);} //update maximum error
                }
            }
        }
        n++;
        fprintf(stderr,"\r[%3d][%5.4f]",n,e_max);
    }while(e_max>eps && n<n_max);
    printf("\nGauss-Seidel method iteration>>%d\n",n);

    plot_results(imax,jmax,i1,i2,jbox,xmin,xmax,ymin,ymax,psi,xy);

    return 0;
}

void plot_results(int imax,int jmax,int i1,int i2,int jbox,double xmin,double xmax,double ymin,double ymax,double psi[imax+1][jmax+1],double xy[imax+1][jmax+1][2]){
    FILE *gp; //gp:gnuplot pointer
    double xy_contour[2][2],psi0;
    char plot1[256],plot2[256];
    int i,j,k,n;
    int *ip,*jp;

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){ //Start gnuplot with pipe
        fprintf(stderr,"Not Found %s.",GNUPLOT_PATH); //Standard error output
        exit(EXIT_FAILURE); //Program abnormal termination
    }

    //Send comands to gnuplot
    fprintf(gp,"set xrange [%f:%f]\n",xmin,xmax); //Set ranges
    fprintf(gp,"set yrange [%f:%f]\n",ymin,ymax);
    fprintf(gp,"set xtics 0.1\nset ytics 0.1\n");
    fprintf(gp,"set grid \n");
    fprintf(gp,"set size ratio %f\n",(ymax-ymin)/(xmax-xmin)); //Set aspect ratio

    //Plot graphic
    strcpy(plot1,"'-' with lines linewidth 2");
    strcpy(plot2,"'-' with lines linewidth 2");
    fprintf(gp, "plot %s, %s\n",plot1,plot2);

    //Plot Wall
    fprintf(gp,"%f\t%f\n",xy[0][0][0],xy[0][0][1]);
    fprintf(gp,"%f\t%f\n",xy[i1][0][0],xy[i1][0][1]);
    fprintf(gp,"%f\t%f\n",xy[i1][jbox][0],xy[i1][jbox][1]);
    fprintf(gp,"%f\t%f\n",xy[i2][jbox][0],xy[i2][jbox][1]);
    fprintf(gp,"%f\t%f\n",xy[i2][0][0],xy[i2][0][1]);
    fprintf(gp,"%f\t%f\n",xy[imax][0][0],xy[imax][0][1]);
    fprintf(gp,"e\n"); //End of array

    //Plot streamline
    ip=&i;jp=&j;
    for(k=1;k<=10;k++){
        i=0;j=0;
        psi0=0.1*(double)k;
        while(i<imax || j<jmax){
            contour(imax,jmax,psi0,psi,xy,xy_contour,ip,jp);
            fprintf(gp,"%f\t%f\n",xy_contour[0][0],xy_contour[0][1]);
            fprintf(gp,"%f\t%f\n",xy_contour[1][0],xy_contour[1][1]);
        }
        fprintf(gp,"\n");
    }
    fprintf(gp,"e\n"); //End of array


    fflush(gp); //Spit out the data stored in the buffer (required)
    system("pause");
    fprintf(gp, "exit\n"); //Terminate gnuplot
    _pclose(gp); //Close Gnuplot

    return;
}

void contour(int imax,int jmax,double z0,double z[imax+1][jmax+1],double xy[imax+1][jmax+1][2],double cntr[][2],int *ip,int *jp){ //Draw contour of z0
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    double x_tmp[4],y_tmp[4];
    bool b1,b2,b3,b4;
    int i,j,k,n,num1=0,num2=0;
    
    //Scan the Edges
    for(i=0;i<imax;i++){
        if(i>=*ip){
            for(j=0;j<jmax;j++){
                if(!(i==*ip && j<=*jp)){
                    x1=xy[i][j][0];y1=xy[i][j][1];z1=z[i][j];
                    x2=xy[i+1][j][0];y2=xy[i+1][j][1];z2=z[i+1][j];
                    x3=xy[i+1][j+1][0];y3=xy[i+1][j+1][1];z3=z[i+1][j+1];
                    x4=xy[i][j+1][0];y4=xy[i][j+1][1];z4=z[i][j+1];
                    b1=(z1<=z0 && z0<=z2) || (z2<=z0 && z0<=z1);
                    b2=(z2<=z0 && z0<=z3) || (z3<=z0 && z0<=z2);
                    b3=(z3<=z0 && z0<=z4) || (z4<=z0 && z0<=z3);
                    b4=(z4<=z0 && z0<=z1) || (z1<=z0 && z0<=z4);
                    if(b1){
                        x_tmp[0]=x1+((z0-z1)/(z2-z1))*(x2-x1);
                        y_tmp[0]=y1;
                    }
                    if(b2){
                        x_tmp[1]=x2;
                        y_tmp[1]=y2+((z0-z2)/(z3-z2))*(y3-y2);
                    }
                    if(b3){
                        x_tmp[2]=x3+((z0-z3)/(z4-z3))*(x4-x3);
                        y_tmp[2]=y3;
                    }
                    if(b4){
                        x_tmp[3]=x4;
                        y_tmp[3]=y4+((z0-z4)/(z1-z4))*(y1-y4);
                    }
                    if(b1 && b2){num1=0;num2=1;}
                    if(b2 && b3){num1=1;num2=2;}
                    if(b3 && b4){num1=2;num2=3;}
                    if(b4 && b1){num1=3;num2=0;}
                    if(b1 && b3){num1=0;num2=2;}
                    if(b4 && b2){num1=3;num2=1;}
                    if(num1!=num2){
                        cntr[0][0]=x_tmp[num1];
                        cntr[0][1]=y_tmp[num1];
                        cntr[1][0]=x_tmp[num2];
                        cntr[1][1]=y_tmp[num2];
                        goto OUT;
                    }
                }
            }
        }
    }
    OUT:
    *ip=i;*jp=j;

    return;
}