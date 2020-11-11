#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<time.h>
#include"solver.h" //Header file of skyline_method
#define GNUPLOT_PATH "C:/gnuplot/bin/gnuplot.exe" //PATH for gnuplot.exe

void Input_data1(int *imax,int *jmax,int *ib,int *jb,int *dof_elm,double *xmin,double *xmax,double *ymin,double *ymax,int *nodes,int *elm){
    FILE *input;

    input=fopen("FEM_data.dat","r");
    fscanf(input,"%d,%d,%d,%d",imax,jmax,ib,jb);
    fscanf(input,"%lf,%lf,%lf,%lf",xmin,xmax,ymin,ymax);
    fscanf(input,"%d,%d,%d,%d",nodes,elm,dof_elm);
    fclose(input);

    return;
}

void Input_data2(int nodes,int elm,double **xy,int **cnct,double *f_S1,bool *b_S1,double *f_S2){
    FILE *input;
    char tmp[64];
    int i,j,k,n,e;

    input=fopen("FEM_data.dat","r");
    for(i=0;i<3;i++){fgets(tmp,sizeof(tmp),input);} //Skip rows 1~3
    for(n=0;n<=nodes;n++){fscanf(input,"%lf,%lf",&xy[n][0],&xy[n][1]);}
    for(e=0;e<=elm;e++){fscanf(input,"%d,%d,%d",&cnct[e][0],&cnct[e][1],&cnct[e][2]);}
    for(n=0;n<=nodes;n++){fscanf(input,"%lf,%d,%lf",&f_S1[n],&b_S1[n],&f_S2[n]);}
    fclose(input);

    return;
}

void Make_nd_nmin(int nodes,int elm,int dof_elm,int *Ksize,int **cnct,int *nd,int *nmin){
    int i,j,k,n,e;

    //Make nmin
    for(e=0;e<=elm;e++){
        for(i=0;i<dof_elm;i++){
            for(j=0;j<dof_elm;j++){
                if(nmin[cnct[e][i]]>cnct[e][j]){
                    nmin[cnct[e][i]]=cnct[e][j];
                }
            }
        }
    }
    //Make Ksize & nd
    k=0;
    for(i=0;i<=nodes;i++){
        for(j=i;j>=nmin[i];j--){
            if(i==j){nd[i]=k;}
            k++;
        }
    }
    k--;
    *Ksize=k;

    return;
}

void Make_K(int nodes,int elm,int dof_elm,int Ksize,int **cnct,int *nd,int *nmin,double **xy,double *K){
    double Ke[elm][dof_elm][dof_elm];
    double xx[dof_elm],yy[dof_elm],Area;
    double x1,x2,x3,y1,y2,y3;
    int i,j,k,n,e,Ksize0,Ki,Kj;

    //Make Ke
    for(e=0;e<=elm;e++){
        //x1~x3,y1~y3
        x1=xy[cnct[e][0]][0];x2=xy[cnct[e][1]][0];x3=xy[cnct[e][2]][0];
        y1=xy[cnct[e][0]][1];y2=xy[cnct[e][1]][1];y3=xy[cnct[e][2]][1];
        //xx,yy
        xx[0]=-(x2-x3);xx[1]=-(x3-x1);xx[2]=-(x1-x2);
        yy[0]=y2-y3;yy[1]=y3-y1;yy[2]=y1-y2;
        Area=0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y3));
        //Ke       
        for(i=0;i<dof_elm;i++){
            for(j=0;j<dof_elm;j++){
                Ke[e][i][j]=(xx[i]*xx[j]+yy[i]*yy[j])/(4*Area);
            }
        }
    }    

    //Sum Ke -> Make K
    for(e=0;e<=elm;e++){
        for(i=0;i<dof_elm;i++){
            for(j=0;j<dof_elm;j++){
                Ki=cnct[e][i];Kj=cnct[e][j];
                if(Ki>=Kj){
                    n=nd[Ki]+Ki-Kj;
                    K[n]=K[n]+Ke[e][i][j];
                }
            }
        }
    }

    return;
}

void Set_BoundaryCondition(int nodes,int *nd,int *nmin,bool *b_S1,double *K,double *f_S1,double *f_S2,double *(*F)){
    int i,j,k,ji,ij,n;
    
    for(i=0;i<=nodes;i++){ //Loop for row
        if(b_S1[i]){ //f_S1[i] is Essential B.C.
            for(j=i+1;j<=nodes;j++){ //Loop for Lji (j>i)
                ji=nd[j]+j-i;
                if(nmin[j]<=i){
                    f_S2[j]=f_S2[j]-f_S1[i]*K[ji];
                }else{
                    f_S2[j]=f_S2[j];
                }
            }
            for(j=nmin[i];j<i;j++){ //Loop for Uji=Lij (j<i)
                ij=nd[i]+i-j;
                f_S2[j]=f_S2[j]-f_S1[i]*K[ij];
            }
            //Set the component of a row/column where a Essential boundary condition exists to 0
            for(j=i+1;j<=nodes;j++){ //Loop for Lji (j>i)
                ji=nd[j]+j-i;
                if(nmin[j]<=i){
                    K[ji]=0;
                }
            }
            for(j=nmin[i];j<i;j++){ //Loop for Uji=Lij (j<i)
                ij=nd[i]+i-j;
                K[ij]=0;
            }
            K[nd[i]]=1;
            f_S2[i]=f_S1[i];
        }
    }

    *F=f_S2;

    return;
}

void Make_uv(int elm,int dof_elm,int **cnct,double **xy,double *F,double **uv){
    double **uv0,*uv_base;
    double xx[dof_elm],yy[dof_elm],Area;
    double x1,x2,x3,y1,y2,y3,p1,p2,p3;
    int n1,n2,n3;
    int e;
    
    //Make uv
    for(e=0;e<=elm;e++){
        //n1~n3,x1~x3,y1~y3,p1~p3
        n1=cnct[e][0];n2=cnct[e][1];n3=cnct[e][2];
        x1=xy[n1][0];x2=xy[n2][0];x3=xy[n3][0];
        y1=xy[n1][1];y2=xy[n2][1];y3=xy[n3][1];
        p1=F[n1];p2=F[n2];p3=F[n3];
        //xx,yy
        xx[0]=-(x2-x3);xx[1]=-(x3-x1);xx[2]=-(x1-x2);
        yy[0]=y2-y3;yy[1]=y3-y1;yy[2]=y1-y2;
        Area=0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y3));
        //uv
        uv[e][0]=(1/(2*Area))*(xx[0]*p1+xx[1]*p2+xx[2]*p3); //u
        uv[e][1]=-(1/(2*Area))*(yy[0]*p1+yy[1]*p2+yy[2]*p3); //v
    }

    return;
}

void Output_results(int imax,int jmax,int ib,int jb,int nodes,int elm,int dof_elm,int **cnct,double xmin,double xmax,double ymin,double ymax,double *F,double **xy,double **uv){
    FILE *output;
    int n,e;
    
    //Output mesh data
    output=fopen("FEM_result.dat","w");
    fprintf(output,"%d,%d,%d,%d\n",imax,jmax,ib,jb);
    fprintf(output,"%f,%f,%f,%f\n",xmin,xmax,ymin,ymax);
    fprintf(output,"%d,%d,%d\n",nodes,elm,dof_elm);
    for(n=0;n<=nodes;n++){fprintf(output,"%f,%f,%f\n",xy[n][0],xy[n][1],F[n]);}
    for(e=0;e<=elm;e++){fprintf(output,"%d,%d,%d,%f,%f\n",cnct[e][0],cnct[e][1],cnct[e][2],uv[e][0],uv[e][1]);}
    fclose(output);

    return;
}

void contour(int elm,int **cnct,double z0,double *z,double **xy,double cntr[][2],int *ep){ //Draw contour of z0
    double x1,x2,x3,y1,y2,y3,z1,z2,z3;
    double x_tmp[3],y_tmp[3];
    bool b1,b2,b3;
    int n1,n2,n3;
    int i,j,k,n,e,num1=0,num2=0;
    
    //Scan the Edges
    for(e=*ep;e<=elm;e++){
        n1=cnct[e][0];n2=cnct[e][1];n3=cnct[e][2];
        x1=xy[n1][0];y1=xy[n1][1];z1=z[n1];
        x2=xy[n2][0];y2=xy[n2][1];z2=z[n2];
        x3=xy[n3][0];y3=xy[n3][1];z3=z[n3];
        b1=(z1<=z0 && z0<=z2) || (z2<=z0 && z0<=z1);
        b2=(z2<=z0 && z0<=z3) || (z3<=z0 && z0<=z2);
        b3=(z3<=z0 && z0<=z1) || (z1<=z0 && z0<=z3);
        if(b1){
            x_tmp[0]=x1+((z0-z1)/(z2-z1))*(x2-x1);
            y_tmp[0]=y1+((z0-z1)/(z2-z1))*(y2-y1);
        }
        if(b2){
            x_tmp[1]=x2+((z0-z2)/(z3-z2))*(x3-x2);
            y_tmp[1]=y2+((z0-z2)/(z3-z2))*(y3-y2);
        }
        if(b3){
            x_tmp[2]=x3+((z0-z3)/(z1-z3))*(x1-x3);
            y_tmp[2]=y3+((z0-z3)/(z1-z3))*(y1-y3);
        }
        if(b1 && b2){num1=0;num2=1;}
        if(b2 && b3){num1=1;num2=2;}
        if(b3 && b1){num1=2;num2=0;}
        if(num1!=num2){
            cntr[0][0]=x_tmp[num1];
            cntr[0][1]=y_tmp[num1];
            cntr[1][0]=x_tmp[num2];
            cntr[1][1]=y_tmp[num2];
            break;
        }
    }
    
    *ep=e+1;
    return;
}

void plot_results(int imax,int jmax,int ib,int jb,int nodes,int elm,int **cnct,double xmin,double xmax,double ymin,double ymax,double *z,double **xy){
    FILE *gp; //gp:gnuplot pointer
    double xy_contour[2][2],z0,zmin=0,zmax=0;
    char plot1[256],plot2[256],title[64]="Potential flow oever backward facing step";
    int i,j,k,n,e,nb,eb,j0,kmax=20;
    int *ep;

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){exit(EXIT_FAILURE);}

    //Send comands to gnuplot
    fprintf(gp,"set title '%s'\n",title);
    fprintf(gp,"set xrange [%f:%f]\n",xmin,xmax); //Set ranges
    fprintf(gp,"set yrange [%f:%f]\n",ymin,ymax);
    fprintf(gp,"set xtics 0.1\nset ytics 0.1\n");
    fprintf(gp,"set grid \n");
    fprintf(gp,"unset key\n"); //Hide legend
    fprintf(gp,"set size ratio %f\n",(ymax-ymin)/(xmax-xmin)); //Set aspect ratio

    //Plot graphic
    strcpy(plot1,"'-' with lines linewidth 2");
    strcpy(plot2,"'-' with lines linewidth 2");
    fprintf(gp, "plot %s, %s\n",plot1,plot2);

    //Plot Wall
    for(i=0;i<=ib;i++){
        nb=ib*jb;if(i<ib){nb=(i+1)*jb;}
        n=i*jmax+i+jb-nb;
        fprintf(gp,"%f\t%f\n",xy[n][0],xy[n][1]);
    }
    fprintf(gp,"\n");
    for(j=0;j<=jb;j++){
        n=ib*jmax+ib+j-ib*jb;
        fprintf(gp,"%f\t%f\n",xy[n][0],xy[n][1]);
    }
    fprintf(gp,"\n");
    for(i=ib;i<=imax;i++){
        n=i*jmax+i+0-ib*jb;
        fprintf(gp,"%f\t%f\n",xy[n][0],xy[n][1]);
    }
    fprintf(gp,"e\n"); //End of array

    //Plot contour
    for(n=0;n<=nodes;n++){
        if(zmin>z[n]){zmin=z[n];}
        if(zmax<z[n]){zmax=z[n];}
    }
    ep=&e;
    for(k=1;k<kmax;k++){
        e=0;    
        z0=zmin+(zmax-zmin)*((double)k/(double)kmax);
        while(e<=elm){
            contour(elm,cnct,z0,z,xy,xy_contour,ep);
            fprintf(gp,"%f\t%f\n",xy_contour[0][0],xy_contour[0][1]);
            fprintf(gp,"%f\t%f\n",xy_contour[1][0],xy_contour[1][1]);
            fprintf(gp,"\n");
        }
    }
    fprintf(gp,"e\n"); //End of array

    fflush(gp); //Spit out the data stored in the buffer (required)
    system("pause");
    fprintf(gp, "exit\n"); //Terminate gnuplot
    _pclose(gp); //Close Gnuplot

    return;
}

void plot_sparse(int nmax,int *nmin){
    FILE *gp; //gp:gnuplot pointer
    double xy_contour[2][2],z0,zmin=0,zmax=0;
    char plot1[64],plot2[64],plot3[64],title[64]="Sparse Matrix";
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
    system("pause");
    fprintf(gp,"exit\n");
    _pclose(gp);

    return;
}

void plot_vector(int imax,int jmax,int ib,int jb,int nodes,int elm,int dof_elm,int **cnct,double xmin,double xmax,double ymin,double ymax,double **xy,double **uv){
    FILE *gp; //gp:gnuplot pointer
    double coef=0.01;
    double xy_contour[2][2],xg,yg;
    char plot1[256],plot2[256],title[64]="Potential flow oever backward facing step";
    int i,j,k,n,e,nb,eb,j0,kmax=20;
    int *ep;

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){exit(EXIT_FAILURE);}

    //Send comands to gnuplot
    fprintf(gp,"set title '%s'\n",title);
    fprintf(gp,"set xrange [%f:%f]\n",xmin,xmax); //Set ranges
    fprintf(gp,"set yrange [%f:%f]\n",ymin,ymax);
    fprintf(gp,"set xtics 0.1\nset ytics 0.1\n");
    fprintf(gp,"set grid \n");
    fprintf(gp,"unset key\n"); //Hide legend
    fprintf(gp,"set size ratio %f\n",(ymax-ymin)/(xmax-xmin)); //Set aspect ratio

    //Plot vector
    for(e=0;e<=elm;e++){
        xg=0;yg=0;
        for(n=0;n<dof_elm;n++){
            xg=xg+xy[cnct[e][n]][0];
            yg=yg+xy[cnct[e][n]][1];
        }
        xg=xg/(double)dof_elm;yg=yg/(double)dof_elm;
        printf("set arrow from %f,%f to %f,%f\n",xg,yg,xg+coef*uv[e][0],yg+coef*uv[e][1]);
        fprintf(gp,"set arrow from %f,%f to %f,%f\n",xg,yg,xg+coef*uv[e][0],yg+coef*uv[e][1]);
    }

    //Plot graphic
    strcpy(plot1,"'-' with lines linewidth 2");
    fprintf(gp, "plot %s\n",plot1);

    //Plot Wall
    for(i=0;i<=ib;i++){
        nb=ib*jb;if(i<ib){nb=(i+1)*jb;}
        n=i*jmax+i+jb-nb;
        fprintf(gp,"%f\t%f\n",xy[n][0],xy[n][1]);
    }
    fprintf(gp,"\n");
    for(j=0;j<=jb;j++){
        n=ib*jmax+ib+j-ib*jb;
        fprintf(gp,"%f\t%f\n",xy[n][0],xy[n][1]);
    }
    fprintf(gp,"\n");
    for(i=ib;i<=imax;i++){
        n=i*jmax+i+0-ib*jb;
        fprintf(gp,"%f\t%f\n",xy[n][0],xy[n][1]);
    }
    fprintf(gp,"e\n"); //End of array

    fflush(gp); //Spit out the data stored in the buffer (required)
    system("pause");
    fprintf(gp, "exit\n"); //Terminate gnuplot
    _pclose(gp); //Close Gnuplot

    return;
}

int main(void){
    FILE *input;
    clock_t start,end;
    int imax,jmax,ib,jb;
    double xmin,xmax,ymin,ymax;
    int nodes,elm,dof_elm,Ksize; 
    double **xy,*xy_base,*f_S1,*f_S2;
    double *K,*F,**uv,*uv_base;
    bool *b_S1;
    int **cnct,*cnct_base,*nd,*nmin;
    int i,j,k,n,e;

    start=clock();
    Input_data1(&imax,&jmax,&ib,&jb,&dof_elm,&xmin,&xmax,&ymin,&ymax,&nodes,&elm);

    /*--Allocate Matrix--*/
    xy=(double **)malloc(sizeof(double *)*(nodes+1));xy_base=(double *)malloc(sizeof(double)*(nodes+1)*2);for(n=0;n<=nodes;n++){xy[n]=xy_base+2*n;}
    cnct=(int **)malloc(sizeof(int *)*(elm+1));cnct_base=(int *)malloc(sizeof(int)*(elm+1)*(dof_elm));for(e=0;e<=elm;e++){cnct[e]=cnct_base+dof_elm*e;}
    uv=(double **)malloc(sizeof(double *)*(elm+1));uv_base=(double *)malloc(sizeof(double)*(elm+1)*2);for(e=0;e<=elm;e++){uv[e]=uv_base+2*e;}
    f_S1=(double *)malloc(sizeof(double)*(nodes+1));
    f_S2=(double *)malloc(sizeof(double)*(nodes+1));
    b_S1=(bool *)malloc(sizeof(bool)*(nodes+1));
    nmin=(int *)malloc(sizeof(int)*(nodes+1));for(n=0;n<=nodes;n++){nmin[n]=n;}
    nd=(int *)malloc(sizeof(int)*(nodes+1));
    
    Input_data2(nodes,elm,xy,cnct,f_S1,b_S1,f_S2);
    Make_nd_nmin(nodes,elm,dof_elm,&Ksize,cnct,nd,nmin);

    /*--Allocate Matrix--*/
    K=(double *)malloc(sizeof(double)*(Ksize+1));for(n=0;n<=Ksize;n++){K[n]=0;}

    Make_K(nodes,elm,dof_elm,Ksize,cnct,nd,nmin,xy,K);
    Set_BoundaryCondition(nodes,nd,nmin,b_S1,K,f_S1,f_S2,&F);
    skyline_method(nodes,Ksize,nd,nmin,K,F);
    Make_uv(elm,dof_elm,cnct,xy,F,uv);
    end=clock();
    
    /*--Output resluts--*/
    Output_results(imax,jmax,ib,jb,nodes,elm,dof_elm,cnct,xmin,xmax,ymin,ymax,F,xy,uv);
    plot_results(imax,jmax,ib,jb,nodes,elm,cnct,xmin,xmax,ymin,ymax,F,xy);
    plot_vector(imax,jmax,ib,jb,nodes,elm,dof_elm,cnct,xmin,xmax,ymin,ymax,xy,uv);
    plot_sparse(nodes,nmin);
    printf("\n/--Specification--/\nnodes=%3d element=%3d Ksize=%3d\n",nodes,elm,Ksize);
    printf("Time >> %f\n",(end-start)/(double)CLOCKS_PER_SEC);
    
    return 0;    
}