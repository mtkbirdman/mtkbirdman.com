void Gaussian_Elimination(int nmax,double **A,double *b);
void skyline_method(int nmax,int ssize,int nd[nmax+1],int jmin[nmax+1],double L[ssize+1],double b[nmax+1]);
void Transpose(int imax,int jmax,double **A,double **AT);
void Inverse(int nmax,double **A,double **Ainv);
void Multiplication(int nmax,double **A,double **B,double **AB);
void LeastSquareMethod(int n_approx,int nmax,double **A,double **B, double **x);