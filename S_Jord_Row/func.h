void InputMatrix (int k, int n, double *a);
double Enteringfunction(int n, int k, int i, int j);
int fileInputMatrix(int n, char *name, double *a);
void InputVector (int n, double *a, double *b);
void OutputMatrix (int l, int n, int m, double *a);
double SolutionError(int n, double* a, double* b, double* x);
double SolutionAccuracy(int n, double* x);
int SolveSystem (int n, int *mass, double *a, double *x, double *answer);