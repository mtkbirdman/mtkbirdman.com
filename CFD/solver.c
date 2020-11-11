#include<stdio.h>
#include<math.h>
#include <stdlib.h>

void Gaussian_Elimination(int nmax,double A[nmax+1][nmax+1],double b[nmax+1]){
    double pivot,p;
    int i,j,ii,jj;

    for(i=0;i<=nmax;i++){
        pivot=A[i][i];
        for(j=i;j<=nmax;j++){
            A[i][j]=A[i][j]/pivot;
        }
        b[i]=b[i]/pivot;
        for(ii=i+1;ii<=nmax;ii++){
            pivot=A[ii][i];
            for(jj=i;jj<=nmax;jj++){
                A[ii][jj]=A[ii][jj]-pivot*A[i][jj];
            }
            b[ii]=b[ii]-pivot*b[i];
        }
    }

    for(i=nmax-1;i>=0;i--){
        for(j=i+1;j<=nmax;j++){
            b[i]=b[i]-A[i][j]*b[j];
        }
    }

    return;
}

int max(int a,int b){
    if(a>=b){return a;}
    if(a<b){return b;}
}

void skyline_method(int nmax,int ssize,int nd[nmax+1],int nmin[nmax+1],double L[ssize+1],double b[nmax+1]){
    double tmp;
    int i,j,k,ik,jk,ij,ki;

    //LDU decomposition
    for(i=0;i<=nmax;i++){
        if(i==0){
            L[0]=L[0]; //D11=K11
        }else if(i==1){
            if(nd[2]==3){
                L[2]=L[2]/L[0]; //L21=K21/D11
                L[1]=L[1]-L[2]*L[0]*L[2]; //D22=K22-L21*D11*U12
            }else{
                L[1]=L[1];
            }
        }else{ //i>2
            for(j=nmin[i];j<i;j++){
                tmp=0;
                for(k=max(nmin[i],nmin[j]);k<j;k++){ //i>j
                    ik=nd[i]+i-k;
                    jk=nd[j]+j-k;
                    tmp=tmp+L[ik]*L[jk]; //SUM{(Lik*Dkk)*Ljk}
                }
                ij=nd[i]+i-j;
                L[ij]=L[ij]-tmp; //(Lij*Djj)=Kij-SUM{(Lik*Dkk)*Ljk}
            }
            for(j=nmin[i];j<i;j++){
                ij=nd[i]+i-j;
                L[ij]=L[ij]/L[nd[j]]; //Lij=(Lij*Djj)/Djj
            }
            tmp=0;
            for(k=nmin[i];k<i;k++){
                ik=nd[i]+i-k;
                tmp=tmp+L[ik]*L[ik]*L[nd[k]]; //SUM{Lik*Lik*Dkk}
            }
            L[nd[i]]=L[nd[i]]-tmp; //Dii=Kii-SUM{Lik*Lik*Dkk}
        }
    }
    
    //Forward substitution
    for(i=0;i<=nmax;i++){
        if(i==0){
            b[i]=b[i]; //y0=b0
        }else{
            for(k=nmin[i];k<i;k++){
                ik=nd[i]+i-k;
                b[i]=b[i]-L[ik]*b[k]; //yi=bi-SUM{Lik*xk}
            }
        }
    }

    //Backward elimination
    for(i=0;i<=nmax;i++){
        b[i]=b[i]/L[nd[i]]; //xi=yi/Dii
    }
    for(i=nmax;i>=0;i--){
        if(i==nmax){
            b[i]=b[i]; //bn=xn
        }else{
            for(k=nmax;k>i;k--){
                ki=nd[k]+k-i;
                if(nmin[k]<=i){
                    b[i]=b[i]-L[ki]*b[k]; //xj=xj-Lij*bi
                }else{
                    b[i]=b[i];
                }
            }
        }
    }

    return;
}

void Transpose(int imax,int jmax,double **A,double **AT){
    int i,j;

    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            AT[j][i]=A[i][j];
        }
    }

    return;
}

void Inverse(int nmax,double **A,double **Ainv){
    double pivot,tmp;
    double B[nmax+1][2*(nmax+1)];
    int i,j,k,ip;


    for(i=0;i<=nmax;i++){
        for(j=0;j<=nmax;j++){
            //Put the matrix you want to inverse in the left half of the matrix B
            B[i][j]=A[i][j];
            //put a unit matrix in the right half of the matrix B
            B[i][(nmax+1)+j]=0;
            if(i==j){B[i][(nmax+1)+j]=1;}
        }
    }

    //Basic operation of the sweeping method
    for(k=0;k<=nmax;k++){
        pivot=0;ip=0;
        //Find the row with the highest absolute value after row k
        for(i=k;i<=nmax;i++){
            if(fabs(B[i][k])>pivot){
                pivot=fabs(B[i][k]);
                ip=i;
            }
        }
        //Swapping the components of row k and row ip
        if(ip!=k){
            for(j=k;j<=2*nmax+1;j++){
                tmp=B[k][j];
                B[k][j]=B[ip][j];
                B[ip][j]=tmp;
            }
        }
        //Divide the diagonal component after row k by the diagonal component of row k
        tmp=B[k][k];
        for(j=k;j<=2*nmax+1;j++){
            B[k][j]=B[k][j]/tmp;
        }
        for(i=0;i<=nmax;i++){
            if(i!=k){
                tmp=B[i][k];
                for(j=k;j<=2*nmax+1;j++){
                    B[i][j]=B[i][j]-tmp*B[k][j];
                }
            }
        }
    }

    for(i=0;i<=nmax;i++){
        for(j=0;j<=nmax;j++){
            Ainv[i][j]=B[i][(nmax+1)+j];
        }
    }

}

void Multiplication(int imax,int jmax,int kmax,double **A,double **B,double **AB){
    int i,j,k;

    for(i=0;i<=imax;i++){for(k=0;k<=kmax;k++){AB[i][k]=0;}}

    for(i=0;i<=imax;i++){
        for(k=0;k<=kmax;k++){
            for(j=0;j<=jmax;j++){
                AB[i][k]=AB[i][k]+A[i][j]*B[j][k];
            }
        }
    }

    return;
}

void LeastSquareMethod(int n_approx,int nmax,double **A,double **B, double **x){
    double **M,**MT,**MTM,**MTMinv,**MTMinvMT;
    double *M_base,*MT_base,*MTM_base,*MTMinv_base,*MTMinvMT_base;
    int i,j,k,n;

    M=(double **)malloc(sizeof(double *)*(nmax+1));M_base=(double *)malloc(sizeof(double)*(nmax+1)*(n_approx+1));for(n=0;n<=nmax;n++){M[n]=M_base+(n_approx+1)*n;}
    MT=(double **)malloc(sizeof(double *)*(n_approx+1));MT_base=(double *)malloc(sizeof(double)*(n_approx+1)*(nmax+1));for(n=0;n<=n_approx;n++){MT[n]=MT_base+(nmax+1)*n;}
    MTM=(double **)malloc(sizeof(double *)*(n_approx+1));MTM_base=(double *)malloc(sizeof(double)*(n_approx+1)*(n_approx+1));for(n=0;n<=n_approx;n++){MTM[n]=MTM_base+(n_approx+1)*n;}
    MTMinv=(double **)malloc(sizeof(double *)*(n_approx+1));MTMinv_base=(double *)malloc(sizeof(double)*(n_approx+1)*(n_approx+1));for(n=0;n<=n_approx;n++){MTMinv[n]=MTMinv_base+(n_approx+1)*n;}
    MTMinvMT=(double **)malloc(sizeof(double *)*(n_approx+1));MTMinvMT_base=(double *)malloc(sizeof(double)*(n_approx+1)*(nmax+1));for(n=0;n<=n_approx;n++){MTMinvMT[n]=MTMinvMT_base+(nmax+1)*n;}

    for(i=0;i<=nmax;i++){
        for(j=0;j<=n_approx;j++){
            M[i][j]=pow(A[i][0],j);
        }
    }

    Transpose(nmax,n_approx,M,MT);
    Multiplication(n_approx,nmax,n_approx,MT,M,MTM);
    Inverse(n_approx,MTM,MTMinv);
    Multiplication(n_approx,n_approx,nmax,MTMinv,MT,MTMinvMT);
    Multiplication(n_approx,nmax,0,MTMinvMT,B,x);

    //for(i=0;i<=nmax;i++){for(j=0;j<=n_approx;j++){printf("%8.5f ",M[i][j]);}printf("\n");}
    //for(i=0;i<=n_approx;i++){for(j=0;j<=n_approx;j++){printf("%8.3f ",MTMinv[i][j]);}printf("\n");}


    return;
}

