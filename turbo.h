void turbo(double r_serial[],int r_l,double u_hat[],int u_l,double sigma2);
void bcjr(double r[][N],double E_matrix[][2]);

void depuncture(double serial[],int serial_l,double r1[][N],double r2[][N],int u_l);
void propagate(int interleaver_array[],double q1[][2],double q2[][2],bool direct);
void initialize_q_array(double q[][2]);
void marginalization(double q_arr[][2],double E_matrix[][2],double u_hat[]);
double calculate_diff_squared_norm(double r[],double coded[]);

double B(int l,unsigned char s,double r[][N]);
double F(int l,unsigned char s,double r[][N]);
double E(int l,unsigned char u,double r[][N]);

double q(int l, unsigned char u);
double g(int l,unsigned char y,double r[][N]);

double log_sum(double a,double b);
double log_approx(double x);
double mymax(double x,double y);
