void Check_FileSize(char *file_name,int *nodes);
void Input_Airfoil(char *file_name,int nodes,double _Complex *z_foil);
void Make_S(int nodes,double *l,double _Complex *z_foil,double _Complex **S1,double _Complex **S2);
void Make_A(int nodes,double *l,double **A,double _Complex *z_foil,double _Complex **S1,double _Complex **S2);
void Make_B(int nodes,double *l,double *B,double _Complex U_inf,double _Complex *z_foil);
void Make_uv(int nodes,double *l,double *gamma,double _Complex U_inf,double _Complex *uv,double _Complex z_ref,double _Complex *z_foil);
void Make_AerodynamicCoefficient(int nodes,double alpha,double *Cl,double *Cm,double *l,double *gamma,double *Cp,double _Complex U_inf,double _Complex *z_foil);
void plot_vector(char *title,int imax,int jmax,int nodes,double xmin,double xmax,double ymin,double ymax,double *l,double *gamma,double _Complex U_inf,double _Complex *z_foil);
void plot_contour(char *title,int imax,int jmax,int kmax,int nodes,double xmin,double xmax,double ymin,double ymax,double *l,double *gamma,double _Complex U_inf,double _Complex *z_foil);
void plot_Cp_x(char *title,int nodes,double *Cp,double _Complex *z_foil);
void plot_polar(char *title,char *file_name);
void Analize_foil(char *OutputFile_name,int nodes,double _Complex *z_foil,double U0,double alpha);