#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<time.h>
#include<complex.h>
#include <unistd.h>
#include"solver.h" //Header file of skyline_method
#define GNUPLOT_PATH "C:/gnuplot/bin/gnuplot.exe" //PATH for gnuplot.exe

void Check_FileSize(char *file_name,int *nodes){
	FILE *fp; //fp:file pointer
	char foil_data[256];
	int i,j;

	//Open file
    *nodes=-1;
	if((fp=fopen(file_name,"r"))==NULL){ //Failed
		fprintf(stderr,"Not Found %s.",file_name);
		exit(EXIT_FAILURE);
	}else{ //Success
		while(fgets(foil_data,sizeof(foil_data),fp)!=NULL){
            *nodes=*nodes+1;
        } //Read the data to the last row
	}
	fclose(fp);
    *nodes=*nodes-1;
	
    return;
}

void Input_Airfoil(char *file_name,int nodes,double _Complex *z_foil){
	FILE *fp; //fp:file pointer
	char foil_data[256],*ch;
	int i,n;
	double x,y,val[1024][4];

	//Open file
    n=0;
	if((fp=fopen(file_name,"r"))==NULL){ //Failed
		fprintf(stderr,"Not Found %s.",file_name); //Standard error output
		exit(EXIT_FAILURE); //Program abnormal termination
	}else{ //Success
		while(fgets(foil_data,sizeof(foil_data),fp)!=NULL){ //Read the data to the last row
			ch=strtok(foil_data," \n"); //Separate the first value by a delimiter
			for(i=0;i<256;i++){
				if(ch==NULL){ //End of string
					break;
				}else{
					val[n][i]=atof(ch); //Convert a read value from char to double
				}
				ch=strtok(NULL," \n"); //Separate the next value
			}
			n++; //Update file size
		}
	}
	fclose(fp); //Close file

    for(n=0;n<=nodes;n++){
        x=val[n+1][0];y=val[n+1][1];
        z_foil[n]=x+I*y;
    }
	
    return;
}

void Make_S(int nodes,double *l,double _Complex *z_foil,double _Complex **S1,double _Complex **S2){
    double _Complex z_ref[nodes];
    int i,j,k,n;

    //Make l & z_ref
    for(i=0;i<nodes;i++){
        l[i]=cabsf(z_foil[i+1]-z_foil[i]); //eq(4.28)
        z_ref[i]=(z_foil[i]+z_foil[i+1])/2; //eq(4.29)
    }

    //Make S1 & S2 eq(4.34)
    for(k=0;k<nodes;k++){
        for(j=0;j<nodes;j++){
            S1[k][j]=((I*l[j])/(2*M_PI*(z_foil[j+1]-z_foil[j])))*(-1+(z_foil[j+1]-z_ref[k])/(z_foil[j+1]-z_foil[j])*clogf((z_foil[j+1]-z_ref[k])/(z_foil[j]-z_ref[k])));
            S2[k][j]=((I*l[j])/(2*M_PI*(z_foil[j+1]-z_foil[j])))*(+1-(z_foil[j]-z_ref[k])/(z_foil[j+1]-z_foil[j])*clogf((z_foil[j+1]-z_ref[k])/(z_foil[j]-z_ref[k])));
        }
    }

    return;
}

void Make_A(int nodes,double *l,double **A,double _Complex *z_foil,double _Complex **S1,double _Complex **S2){
    double _Complex z_n[nodes+1]; //n:normal
    int i,j,k,n;

    //Make z_n
    for(i=0;i<nodes;i++){
        z_n[i]=((z_foil[i+1]-z_foil[i])/l[i])*(1/I); //eq(4.29)
    }

    //Make A eq(4.40)
    for(k=0;k<nodes;k++){
        for(j=0;j<=nodes;j++){
            if(j==0){
                A[k][j]=crealf(S1[k][j]*z_n[k]);
            }else if(j==nodes){
                A[k][j]=crealf(S2[k][j-1]*z_n[k]);
            }else{
                A[k][j]=crealf((S2[k][j-1]+S1[k][j])*z_n[k]);
            }
        }
    }
    A[nodes][0]=1;A[nodes][nodes]=1; //Kutta condition

    return;
}

void Make_B(int nodes,double *l,double *B,double _Complex U_inf,double _Complex *z_foil){
    double _Complex z_n[nodes+1]; //n:normal
    int i,j,k,n;

    //Make z_n
    for(i=0;i<nodes;i++){
        z_n[i]=((z_foil[i+1]-z_foil[i])/l[i])*(1/I); //eq(4.29)
    }

    //Make B eq(4.40d)
    for(k=0;k<nodes;k++){
        B[k]=-crealf(U_inf*z_n[k]);
    }
    B[nodes]=0; //Kutta Condition

    return;
}

void Make_uv(int nodes,double *l,double *gamma,double _Complex U_inf,double _Complex *uv,double _Complex z_ref,double _Complex *z_foil){
    double _Complex S1_uv[nodes+1],S2_uv[nodes+1];
    int k;

    //eq.(4.43)
    for(k=0;k<nodes;k++){
        S1_uv[k]=((I*l[k])/(2*M_PI*(z_foil[k+1]-z_foil[k])))*(-1+(z_foil[k+1]-z_ref)/(z_foil[k+1]-z_foil[k])*clogf((z_foil[k+1]-z_ref)/(z_foil[k]-z_ref)));
        S2_uv[k]=((I*l[k])/(2*M_PI*(z_foil[k+1]-z_foil[k])))*(+1-(z_foil[k]-z_ref)/(z_foil[k+1]-z_foil[k])*clogf((z_foil[k+1]-z_ref)/(z_foil[k]-z_ref)));
    }
    //eq.(4.42)
    *uv=U_inf+S1_uv[0]*gamma[0]+S2_uv[nodes-1]*gamma[nodes];
    for(k=1;k<nodes;k++){
        *uv=*uv+(S2_uv[k-1]+S1_uv[k])*gamma[k];
    }

    return;
}

void Make_AerodynamicCoefficient(int nodes,double alpha,double *Cl,double *Cm,double *l,double *gamma,double *Cp,double _Complex U_inf,double _Complex *z_foil){
    /*--Constant--*/
    double scale=1.0+pow(10,-10),rad=M_PI/180;
    double _Complex z_qtr=0.25+I*0; //qtr:quarter
    /*--Variables--*/    
    double _Complex z_ref[nodes];
    double _Complex uv;
    double U0=cabs(U_inf),U;
    double theta,pt,pn; //n:normal t:tangent
    double x,y;
    int i,j,k,n;

    //Make z_ref
    for(n=0;n<nodes;n++){
        z_ref[n]=(z_foil[n]+z_foil[n+1])/2; //eq(4.29)
        z_ref[n]=z_ref[n]*scale;
    }

    //Calculation Cp
    for(n=0;n<nodes;n++){
        Make_uv(nodes,l,gamma,U_inf,&uv,z_ref[n],z_foil);
        U=cabs(uv);
        Cp[n]=(1-(U/U0)*(U/U0));
    }

    //Calculation Cl, Cm
    *Cl=0;*Cm=0;
    for(n=0;n<nodes;n++){
        theta=carg(z_foil[n+1]-z_foil[n]);
        pn=Cp[n]*l[n]*cos(theta);pt=-Cp[n]*l[n]*sin(theta);
        x=crealf(z_ref[n]-z_qtr);y=cimagf(z_ref[n]-z_qtr);
        *Cl=*Cl+pn*cos(rad*alpha)-pt*sin(rad*alpha);
        *Cm=*Cm-pn*x+pt*y;
    }

    return;
}

void Output_resluts(char *OutputFile_name,double alpha,double Cl,double Cm){
    FILE *fp;

    if((fp=fopen(OutputFile_name,"a"))==NULL){ //Failed
		fprintf(stderr,"Not Found %s.",OutputFile_name);exit(EXIT_FAILURE);
	}else{ //Success
        fprintf(fp,"%f\t%f\t%f\n",alpha,Cl,Cm);
	}
	fclose(fp);
    printf("alpha=%6.3f Cl=%6.3f Cm=%6.3f\n",alpha,Cl,Cm);

}

void plot_vector(char *title,int imax,int jmax,int nodes,double xmin,double xmax,double ymin,double ymax,double *l,double *gamma,double _Complex U_inf,double _Complex *z_foil){
	FILE *gp; //gp:gnuplot pointer
    /*--Constant--*/
    double coef=0.01;
    /*--Variables--*/
    char plot1[256],plot2[256];
    double _Complex z_ref,S1_uv[nodes+1],S2_uv[nodes+1];
    double _Complex uv;
    double dx=(xmax-xmin)/(imax),dy=(ymax-ymin)/(jmax);
    double x,y,u,v;
    int i,j,k,n;

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){exit(EXIT_FAILURE);}

	//Send comands to gnuplot
    fprintf(gp,"set title '%s'\n",title);
	fprintf(gp,"set xrange [%f:%f]\n",xmin,xmax); //Set ranges
	fprintf(gp,"set yrange [%f:%f]\n",ymin,ymax);
    fprintf(gp,"unset xtics\nunset ytics\n");
    fprintf(gp,"unset key\n"); //Hide legend
	fprintf(gp,"set size ratio %f\n",(ymax-ymin)/(xmax-xmin)); //Set aspect ratio

    //Plot Vector
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            x=xmin+dx*i;y=ymin+dy*j;
            z_ref=x+y*I;
            Make_uv(nodes,l,gamma,U_inf,&uv,z_ref,z_foil);
            u=crealf(uv),v=-cimagf(uv);
            fprintf(gp,"set arrow from %f,%f to %f,%f\n",x,y,x+u*coef,y+v*coef);
        }
    }

	//Plot graphic
    strcpy(plot1,"'-' with lines linewidth 2");
	fprintf(gp, "plot %s\n",plot1);

    //Plot Airfoil
    for(n=0;n<=nodes;n++){
        fprintf(gp,"%f\t%f\n",crealf(z_foil[n]),cimagf(z_foil[n]));
    }
    fprintf(gp,"e\n"); //End of array

	fflush(gp);system("pause");fprintf(gp, "exit\n");_pclose(gp); //Close Gnuplot

    return;
}

void plot_contour(char *title,int kmax,int nodes,double xmin,double xmax,double ymin,double ymax,double *l,double *gamma,double _Complex U_inf,double _Complex *z_foil){
	FILE *gp; //gp:gnuplot pointer
    /*--Constant--*/
    double dt0=(1.0/cabs(U_inf))*0.01,dtheta_max=100000;
    double rad=M_PI/180;
    /*--Variables--*/
    double _Complex z_ref,S1_uv[nodes+1],S2_uv[nodes+1];
    double _Complex uv;
    double alpha=atan(-cimagf(U_inf)/creal(U_inf)); //AOA
    double dl1,dl2,x0,y0,x1,x2,y1,y2;
    double dt; //time step
    double theta,theta0,dtheta; //angle of stream line
    double x,y,u,v;
    char plot1[256],plot2[256];
    int i,j,k,n;

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){exit(EXIT_FAILURE);}

	//Send comands to gnuplot
    fprintf(gp,"set title '%s'\n",title);
	fprintf(gp,"set xrange [%f:%f]\n",xmin,xmax); //Set ranges
	fprintf(gp,"set yrange [%f:%f]\n",ymin,ymax);
    fprintf(gp,"unset xtics\nunset ytics\n");
    fprintf(gp,"unset key\n"); //Hide legend
	fprintf(gp,"set size ratio %f\n",(ymax-ymin)/(xmax-xmin)); //Set aspect ratio

	//Plot graphic
    strcpy(plot1,"'-' with lines linewidth 2");
    strcpy(plot2,"'-' with lines linewidth 2");
	fprintf(gp, "plot %s,%s\n",plot1,plot2);

    //Plot Airfoil
    for(n=0;n<=nodes;n++){
        fprintf(gp,"%f\t%f\n",crealf(z_foil[n]),cimagf(z_foil[n]));
    }
    fprintf(gp,"e\n"); //End of array
    
    //Plot Contour
    if(alpha>0){
        x0=xmin;y0=ymin;
        dl1=fabs((ymax-ymin)*cos(alpha));dl2=fabs((xmax-xmin)*sin(alpha));
    }else{
        x0=xmin;y0=ymax;
        dl1=fabs((xmax-xmin)*sin(alpha));dl2=fabs((ymax-ymin)*cos(alpha));
    }
    x1=x0-dl1*sin(alpha);y1=y0+dl1*cos(alpha);
    x2=x0+dl2*sin(alpha);y2=y0-dl2*cos(alpha);
    for(k=0;k<=kmax;k++){
        //Set x,y
        x=x1*((double)(k)/(double)(kmax))+x2*((double)(kmax-k)/(double)(kmax));
        y=y1*((double)(k)/(double)(kmax))+y2*((double)(kmax-k)/(double)(kmax));
        fprintf(gp,"%f\t%f\n",x,y);
        //calculate stream line
        theta0=alpha;
        do{
            z_ref=x+y*I;
            Make_uv(nodes,l,gamma,U_inf,&uv,z_ref,z_foil);
            u=crealf(uv),v=-cimagf(uv);
            //Adjust dt
            theta=atan(v/u);
            dtheta=fabs((theta-theta0)/dt);if(dtheta>dtheta_max){dtheta=dtheta_max;}
            dt=dt0*(1/(1+dtheta*dtheta));
            //calculate x,y
            x=x+u*dt;y=y+v*dt;
            fprintf(gp,"%f\t%f\n",x,y);
            theta0=theta;
        }while(x<xmax);
        fprintf(gp,"\n"); //End of array
    }
    fprintf(gp,"e\n"); //End of array

	fflush(gp);system("pause");fprintf(gp, "exit\n");_pclose(gp); //Close Gnuplot

    return;
}

void plot_Cp_x(char *title,int nodes,double *Cp,double _Complex *z_foil){
	FILE *gp; //gp:gnuplot pointer
    /*--Constant--*/
    double xmin=0,xmax=1.0;
    /*--Variables--*/
    double _Complex z_ref[nodes];
    double ymin,ymax;
    char plot1[256],plot2[256];
    int i,j,k,n;

    //Make z_ref, ymin, ymax
    for(n=0;n<nodes;n++){
        z_ref[n]=(z_foil[n]+z_foil[n+1])/2; //eq(4.29)
        if(ymin>Cp[n]){ymin=Cp[n];}
        if(ymax<Cp[n]){ymax=Cp[n];}
    }

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){exit(EXIT_FAILURE);}

	//Send comands to gnuplot
    fprintf(gp,"set title '%s'\n",title);
	fprintf(gp,"set xrange [%f:%f]\n",xmin,xmax); //Set ranges
	fprintf(gp,"set yrange [%f:%f] reverse\n",ymax,ymin);
    fprintf(gp,"unset key\n"); //Hide legend
    fprintf(gp,"set zeroaxis\n");

	//Plot graphic
    strcpy(plot1,"'-' with lines linewidth 2");
	fprintf(gp, "plot %s\n",plot1);

    //Plot Cp vs x
    for(n=0;n<nodes;n++){
        fprintf(gp,"%f\t%f\n",crealf(z_ref[n]),Cp[n]);
    }
    fprintf(gp,"e\n"); //End of array

	fflush(gp);system("pause");fprintf(gp, "exit\n");_pclose(gp); //Close Gnuplot

    return;
}

void plot_polar(char *title,char *file_name){
	FILE *gp; //gp:gnuplot pointer
    /*--Constant--*/
    double xmin=-5,xmax=15,ymin=-0.14,ymax=0.01;
    char plot1[256]="using 1:3 with lines linewidth 2";
    /*--Variables--*/

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){exit(EXIT_FAILURE);}

	//Send comands to gnuplot
    fprintf(gp,"set title '%s'\n",title);
    fprintf(gp,"unset key\n"); //Hide legend
    fprintf(gp,"set zeroaxis\n");
	//fprintf(gp,"set xrange [%f:%f]\n",xmin,xmax); //Set ranges
	fprintf(gp,"set yrange [%f:%f]\n",ymin,ymax);
	fprintf(gp,"set size square\n"); //Set aspect ratio

	//Plot graphic
	fprintf(gp,"plot '%s' %s\n",file_name,plot1);

	fflush(gp);system("pause");fprintf(gp, "exit\n");_pclose(gp); //Close Gnuplot

    return;
}

void Analize_foil(char *OutputFile_name,int nodes,double _Complex *z_foil,double U0,double alpha){
    /*--Constant--*/
    int kmax=40; //kmax:number of streamline
    double xmin=-0.5,xmax=1.5,ymin=-0.5,ymax=0.5; //graph range
    char title[64];
    /*--Variables--*/
    double _Complex **S1,**S2;
    double _Complex *S1_base,*S2_base;
    double _Complex U_inf=U0*cos(alpha*(M_PI/180))-I*U0*sin(alpha*(M_PI/180)); //Uniform flow
    double **A,*B;
    double *gamma,*l,*Cp; //gamma:circulation, l:length, Cp:pressure coefficient
    double *A_base;
    double Cl,Cm; //Cl:Lift coefficient, Cm:piching moment coefficient
    int i,j,k,n;

    /*--Allocate Matrix--*/
    S1=(double _Complex **)malloc(sizeof(double _Complex *)*(nodes+1));S1_base=(double _Complex *)malloc(sizeof(double _Complex)*(nodes+1)*(nodes+1));for(n=0;n<=nodes;n++){S1[n]=S1_base+(nodes+1)*n;}
    S2=(double _Complex **)malloc(sizeof(double _Complex *)*(nodes+1));S2_base=(double _Complex *)malloc(sizeof(double _Complex)*(nodes+1)*(nodes+1));for(n=0;n<=nodes;n++){S2[n]=S2_base+(nodes+1)*n;}
    l=(double *)malloc(sizeof(double)*nodes); //value at z_ref
    Cp=(double *)malloc(sizeof(double)*nodes); //value at z_ref
    A=(double **)malloc(sizeof(double *)*(nodes+1));A_base=(double *)malloc(sizeof(double)*(nodes+1)*(nodes+1));for(n=0;n<=nodes;n++){A[n]=A_base+(nodes+1)*n;}
    B=(double *)malloc(sizeof(double)*(nodes+1));
    gamma=(double *)malloc(sizeof(double)*(nodes+1));

    Make_S(nodes,l,z_foil,S1,S2);
    Make_A(nodes,l,A,z_foil,S1,S2);
    Make_B(nodes,l,B,U_inf,z_foil);
    Gaussian_Elimination(nodes,A,B);for(n=0;n<=nodes;n++){gamma[n]=B[n];}
    Make_AerodynamicCoefficient(nodes,alpha,&Cl,&Cm,l,gamma,Cp,U_inf,z_foil);

    Output_resluts(OutputFile_name,alpha,Cl,Cm);

    strcpy(title,"Cp vs x");
    plot_Cp_x(title,nodes,Cp,z_foil);
    strcpy(title,"Stream line");
    plot_contour(title,kmax,nodes,xmin,xmax,ymin,ymax,l,gamma,U_inf,z_foil);

    return;
}