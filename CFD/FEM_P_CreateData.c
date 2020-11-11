#include<stdio.h>
#include<stdbool.h>

int main(void){
    FILE *output;
    /*--Constance--*/
    int imax=6,jmax=4,ib=2,jb=2; //i:couner for x, j:counter for y
    int dof_elm=3; //dof:Degree of Freedom, elm:elments
    double xmin=0,xmax=0.6,ymin=0,ymax=0.4;
    double psi_max=1.0,psi_min=0;
    /*--Variables--*/
    int nodes=(imax+1)*(jmax+1)-1-ib*jb,elm=(imax*jmax-ib*jb)*2-1; //b:box
    int n_S1=2*imax+1+(jb-1),n_S2=(jmax+1)*2-jb; //S1:Essential boundary, S2:Natural boundary
    double xy[nodes+1][2],f_S1[nodes+1],f_S2[nodes+1];
    double f,l; //f:dpsi/dn at S2, l:length of S2
    bool b_S1[nodes+1]; //If S1 is fixed,boolean->True
    int cnct[elm+1][dof_elm]; //cnct:Connectivity
    double dx=(xmax-xmin)/(double)imax,dy=(ymax-ymin)/(double)jmax;
    int i,j,j0,k,n,e,nb,eb;

    //Create Mesh
    n=0;
    for(i=0;i<=imax;i++){
        j0=0;
        if(i<ib){j0=jb;}
        for(j=j0;j<=jmax;j++){
            xy[n][0]=dx*i;
            xy[n][1]=dy*j;
            n++;
        }
    }

    //Create connectivity
    for(i=0;i<imax;i++){
        nb=ib*jb;eb=2*ib*jb;j0=jmax+1;
        if(i<ib){nb=(i+1)*jb;eb=2*(i+1)*jb;}
        if(i<ib-1){j0=(jmax+1-jb);}
        for(j=0;j<jmax;j++){
            e=2*i*jmax+2*j-eb;
            n=i*jmax+i+j-nb;
            if(!(i<ib&&j<jb)){
                cnct[e][0]=n;cnct[e][1]=n+j0;cnct[e][2]=n+1;
                cnct[e+1][0]=n+j0+1;cnct[e+1][1]=n+1;cnct[e+1][2]=n+j0;
            }
        }
    }
    
    //Create Essential boundary conditions(S1) & Natural boundary conditions(S2)
    for(n=0;n<=nodes;n++){f_S1[n]=0;f_S2[n]=0;b_S1[n]=false;}
    for(i=0;i<=imax;i++){
        j0=0;nb=ib*jb;eb=2*ib*jb;
        if(i<ib){nb=(i+1)*jb;eb=2*(i+1)*jb;j0=jb;}
        for(j=j0;j<=jmax;j++){
            //S1
            n=i*jmax+i+j-nb;
            if(j==j0){f_S1[n]=psi_min;b_S1[n]=true;} //psi=0
            if(i==ib&&(j0<j&&j<=jb)){f_S1[n]=psi_min;b_S1[n]=true;} //psi=0
            if(j==jmax){f_S1[n]=psi_max;b_S1[n]=true;} //psi=psi_max*/
            //S2
            e=2*i*jmax+2*j-eb;
            f=0;
            l=dy;
            if((i==0||i==imax-1)&&j<jmax){
                f_S2[cnct[e][2]]=f_S2[cnct[e][2]]+(f*l)/2.0;
                f_S2[cnct[e][0]]=f_S2[cnct[e][0]]+(f*l)/2.0;
            }
        }
    }

    //Output mesh data
    output=fopen("FEM_data.dat","w");
    fprintf(output,"%d,%d,%d,%d\n",imax,jmax,ib,jb);
    fprintf(output,"%f,%f,%f,%f\n",xmin,xmax,ymin,ymax);
    fprintf(output,"%d,%d,%d\n",nodes,elm,dof_elm);
    for(n=0;n<=nodes;n++){fprintf(output,"%f,%f\n",xy[n][0],xy[n][1]);}
    for(e=0;e<=elm;e++){fprintf(output,"%d,%d,%d\n",cnct[e][0],cnct[e][1],cnct[e][2]);}
    for(n=0;n<=nodes;n++){fprintf(output,"%f,%d,%f\n",f_S1[n],b_S1[n],f_S2[n]);}
    fclose(output);

    return 0;
}