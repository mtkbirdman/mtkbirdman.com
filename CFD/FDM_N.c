#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#define GNUPLOT_PATH "C:/gnuplot/bin/gnuplot.exe" //PATH for gnuplot.exe

void Initialize(int imax,int jmax,int i1,int i2,int jbox,double psi_max,double psi[imax+1][jmax+1],double zeta[imax+1][jmax+1]){
    FILE *input;
    double val[4][1024];
    char str[32],*ch;
    int i,j,j0,k,n;

    //Initialize
    printf("Initialize by data file? (y/n)\n>>");
    fgets(str,sizeof(str),stdin); //Input from keyboard
    if(strcmp(str,"y\n")==0){
        strcpy(str,"\0");
        printf("Input file name\n>>");
        fgets(str,sizeof(str),stdin); //Input from keyboard
	    str[strlen(str)-1]='\0'; //Remove "\n"
        if((input=fopen(str,"r"))==NULL){ //Failed
            fprintf(stderr,"Not Found %s.",str); //Standard error output
            exit(EXIT_FAILURE); //Program abnormal termination
        }else{ //Success
            strcpy(str,"\0");
            while(fgets(str,sizeof(str),input)!=NULL){ //Read the data to the last row
                ch=strtok(str," \n"); //Separate the first value by a delimiter
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
        fclose(input); //Close file
        for(i=0;i<=imax;i++){
            for(j=0;j<=jmax;j++){
                psi[i][j]=val[i*jmax+i+j][0];
                zeta[i][j]=val[i*jmax+i+j][1];
            }
        }
    }else{
        for(i=0;i<=imax;i++){
            j0=0;
            if(i1<i&&i<i2){j0=jbox;}
            for(j=j0;j<=jmax;j++){
                psi[i][j]=psi_max*((double)j)/((double)jmax);
                zeta[i][j]=0;
            }
        }
    }

    return;
}

void Create_Mesh(int imax,int jmax,double dx,double dy,double xy[imax+1][jmax+1][2]){
    int i,j,k,n;

    //Create mesh
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            xy[i][j][0]=dx*i;
            xy[i][j][1]=dy*j;
        }
    }

    return;
}

void Output(int imax,int jmax,double psi[imax+1][jmax+1],double zeta[imax+1][jmax+1]){
    FILE *output;
    double val[4][1024];
    char str[32],*ch;
    int i,j,k,n;

    //Output data file
    output=fopen("output_FDM_N.txt","w");
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            fprintf(output,"%f\t%f\n",psi[i][j],zeta[i][j]);
        }
    }
    fclose(output);

    return;
}

void Fixed_BoundaryCondition(int imax,int jmax,int i1,int i2,int jbox,double psi_max,double psi[imax+1][jmax+1],double zeta[imax+1][jmax+1]){
    int i,j;

    //Input fixed boundry conditions
    for(i=0;i<=imax;i++){ psi[i][0]=0; } //bottom wall
    for(j=0;j<=jbox;j++){ psi[i1][j]=0; } //box left-side wall
    for(i=i1;i<=i2;i++){ psi[i][jbox]=0; } //box top wall
    for(j=0;j<=jbox;j++){ psi[i2][j]=0; } //box right-side wall
    for(i=0;i<=imax;i++){ psi[i][jmax]=psi_max; } //free stream
    for(i=0;i<=imax;i++){ zeta[i][jmax]=0; } //free stream

    return;
}

void Variable_BoundaryCondition(int imax,int jmax,int i1,int i2,int jbox,double dx,double dy,double psi[imax+1][jmax+1],double zeta[imax+1][jmax+1]){
    double dx2=dx*dx,dy2=dy*dy,dx2dy2=dx2*dy2;
    int i,j;

    //Input Variable boundry conditions
    for(i=1;i<i1;i++){zeta[i][0]=(2.0/dy2)*(psi[i][0]-psi[i][1]);}
    for(j=1;j<jbox;j++){zeta[i1][j]=(2.0/dy2)*(psi[i1][j]-psi[i1-1][j]);}
    for(i=i1;i<=i2;i++){zeta[i][jbox]=(2.0/dy2)*(psi[i][jbox]-psi[i][jbox+1]);}
    for(j=1;j<jbox;j++){zeta[i2][j]=(2.0/dy2)*(psi[i2][j]-psi[i2-1][j]);}
    for(i=i2+1;i<imax;i++){zeta[i][0]=(2.0/dy2)*(psi[i][0]-psi[i][1]);}
    for(j=0;j<=jmax;j++){
        psi[0][j]=psi[1][j];zeta[0][j]=zeta[1][j];
        psi[imax][j]=psi[imax-1][j];zeta[imax][j]=zeta[imax-1][j];
    }

    return;
}

void Solve_VorticityTransportEquation(int imax,int jmax,int i1,int i2,int jbox,double nu,double dt,double dx,double dy,double *e_max,double psi[imax+1][jmax+1],double zeta[imax+1][jmax+1]){
    double dx2=dx*dx,dy2=dy*dy,dx2dy2=dx2*dy2;
    double dt4dxdy=dt/(4*dx*dy),ndtdx2=(nu*dt)/dx2,ndtdy2=(nu*dt)/dy2;
    double zeta0[imax+1][jmax+1],dzeta;
    int i,j,j0;

    //Copy zeta to zeta0
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            zeta0[i][j]=zeta[i][j];
        }
    }
    //Solve vorticity transport equation (eq.2.26)
    for(i=1;i<imax;i++){
        j0=1;
        if(i1<=i&&i<=i2){j0=jbox+1;}
        for(j=j0;j<jmax;j++){
            dzeta=-dt4dxdy*((psi[i][j+1]-psi[i][j-1])*(zeta0[i+1][j]-zeta0[i-1][j])-(psi[i+1][j]-psi[i-1][j])*(zeta0[i][j+1]-zeta0[i][j-1]))+ndtdx2*(zeta0[i-1][j]-2*zeta0[i][j]+zeta0[i+1][j])+ndtdy2*(zeta0[i][j-1]-2*zeta0[i][j]+zeta0[i][j+1]);
            zeta[i][j]=zeta0[i][j]+dzeta;
            if(fabs(dzeta)>*e_max){*e_max=fabs(dzeta);}
        }
    }

    return;
}

void GauseSeidelMethod(int imax,int jmax,int i1,int i2,int jbox,int np_max,double eps_p,double psi0,double dx,double dy,double psi[imax+1][jmax+1],double zeta[imax+1][jmax+1]){
    double dx2=dx*dx,dy2=dy*dy,dx2dy2=dx2*dy2;
    double r2dx2pdy2=1/(2*(dx2+dy2));
    double dpsi,e_max;
    int i,j,j0,n=0;

    //Gauss-Seidel method
    do{
        e_max=0;
        for(j=0;j<=jmax;j++){
            psi[0][j]=psi[1][j];zeta[0][j]=zeta[1][j];
            psi[imax][j]=psi[imax-1][j];zeta[imax][j]=zeta[imax-1][j];
        }
        for(i=1;i<imax;i++){
            j0=1;
            if(i1<=i&&i<=i2){j0=jbox+1;}
            for(j=j0;j<jmax;j++){
                dpsi=r2dx2pdy2*((psi[i][j-1]+psi[i][j+1])*dx2+(psi[i-1][j]+psi[i+1][j])*dy2+zeta[i][j]*dx2dy2)-psi[i][j];
                psi[i][j]=psi[i][j]+dpsi;
                if(fabs(dpsi)>e_max){e_max=fabs(dpsi);} //update maximum error
            }
        }
        n++;
        e_max=e_max/psi0;
    }while(e_max>eps_p && n<np_max);
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

void plot_results(int imax,int jmax,int i1,int i2,int jbox,double xmin,double xmax,double ymin,double ymax,double z[imax+1][jmax+1],double xy[imax+1][jmax+1][2]){
    FILE *gp; //gp:gnuplot pointer
    double xy_contour[2][2],z0,zmin=0,zmax=0;
    char plot1[256],plot2[256],title[32]="Streamlines at Re=2";
    int i,j,k,n,kmax=30;
    int *ip,*jp;

    //Start Gnuplot comand
    if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){ //Start gnuplot with pipe
        fprintf(stderr,"Not Found %s.",GNUPLOT_PATH); //Standard error output
        exit(EXIT_FAILURE); //Program abnormal termination
    }

    //Send comands to gnuplot
    fprintf(gp,"set title '%s'\n",title);
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

    //Plot contour
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            if(zmin>z[i][j]){zmin=z[i][j];}
            if(zmax<z[i][j]){zmax=z[i][j];}
        }
    }
    ip=&i;jp=&j;
    for(k=1;k<=kmax;k++){
        i=0;j=0;
        z0=zmin+(zmax-zmin)*((double)k/(double)kmax);
        while(i<imax || j<jmax){
            contour(imax,jmax,z0,z,xy,xy_contour,ip,jp);
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

int main(void){
    int imax=20,jmax=10,i1=8,i2=12,jbox=4,n_max=5000,np_max=50000;
    double xmin=0,xmax=2,ymin=0,ymax=1;
    double psi_max=0.5,eps_p=0.001,eps_z=0.001,e_max;
    double cou=0.1; //cou:Courant number
    double nu=0.01; //nu:Kinematic viscosity coefficient
    double dx=(xmax-xmin)/(double)imax,dy=(ymax-ymin)/(double)jmax;
    double xy[imax+1][jmax+1][2];
    double psi[imax+1][jmax+1]; //stream-function
    double zeta[imax+1][jmax+1]; //vorticity
    double V0=psi_max/(dy*jmax); //Uniform flow velocity
    double Re=V0*(dy*jbox)/nu; //Reynolds number
    double dt=cou*(dx/V0); //Time step
    double psi0=psi_max,zeta0=V0/(dy*jmax);
    double dx2=dx*dx,dy2=dy*dy,dx2dy2=dx2*dy2;
    double ez_max,*pt=&ez_max;
    int i,j,k,n=0;
    // Re=2 -> psi_max=0.05
    // Re=20 -> psi_max=0.5

    printf("Re=%f\n",Re);
    printf("dt=%f\n",dt);
    Create_Mesh(imax,jmax,dx,dy,xy);
    Initialize(imax,jmax,i1,i2,jbox,psi_max,psi,zeta);
    Fixed_BoundaryCondition(imax,jmax,i1,i2,jbox,psi_max,psi,zeta);
    do{
        printf("\riteration %d",n);
        ez_max=0;
        Variable_BoundaryCondition(imax,jmax,i1,i2,jbox,dx,dy,psi,zeta);
        Solve_VorticityTransportEquation(imax,jmax,i1,i2,jbox,nu,dt,dx,dy,pt,psi,zeta);
        GauseSeidelMethod(imax,jmax,i1,i2,jbox,np_max,eps_p,psi0,dx,dy,psi,zeta);
        ez_max=ez_max/zeta0;
        n++;
    }while(ez_max>eps_p && n<n_max);
    printf("\n");
    Output(imax,jmax,psi,zeta);
    plot_results(imax,jmax,i1,i2,jbox,xmin,xmax,ymin,ymax,psi,xy);

    return 0;
}
