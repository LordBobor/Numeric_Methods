import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.Math.exp;
import static java.lang.Math.sin;

public class Parabolic {
    static double a, b, c;
    static int n, k;
    static double Tk=1, Xn=Math.PI/2;
    static double h, tau;
    static double alpha;
    static double beta;
    static double gamma;
    static double delta;

    double[] x;
    double[] res;

    public Parabolic(){
        alpha = 1; beta = 0;
        gamma = 0; delta = 1;
    }

    public Parabolic(double a, double b, double c, int n, int k) {
        this.a=a; this.b=b; this.c=c;
        this.n=n; this.k=k;
        alpha = 1; beta = 0;
        gamma = 0; delta = 1;

    }

    static double f(double x){return 0;}
    static double fi_0(double t){return exp((c-a)*t);}
    static double fi_l(double t){return exp((c-a)*t);}
    static double ksi(double x){return sin(x);}



    public double[] Solve(double param){
        h= Xn/n;
        tau = Tk/k;
        n++;
        double[] u = new double[n];
        initSimple(u);
        x = new double[n];

        for(int i =0; i<n; ++i){
            x[i]=i*h;
            System.out.format("%2d) x=%.2f Exact=%.2f u[0]=%.2f delta=%.2f\n", i, i*h,
                    exactSolution(i*h, 0), u[i], exactSolution(i*h, 0)-u[i]);
        }

        System.out.println("n=" + n + " h=" + h);
        System.out.println("k=" + k + " tau=" + tau );
        if(tau <= ((h*h)/(2*a)) ) System.out.println("tau < " + ((h*h)/(2*a)));

        if(param == 0) res = solve_explicit(u);
        else if(param>0 && param <= 1) res = solve_implicit(u, param);
        else throw new IllegalArgumentException();


        for(int i =0; i<n; ++i){ System.out.format("%.3f ", x[i]); }
        System.out.format("\n");
        print(res, k);

        return res;
    }
    /*  public static void main(String[] args){
        a=0.02; b=0; c=-0.04;
        n=10;
        k=10;
        h= Xn/n;
        //tau = h*h/2*a;
        tau = Tk/k;
        alpha = 1; beta = 0;
        gamma = 0; delta = 1;
        double[] u = new double[n+1];
        initSimple(u);

        for(int i =0; i<=n; ++i){
            System.out.format("%2d) x=%.2f Exact=%.2f u[0]=%.2f delta=%.2f\n", i, i*h,
                    exactSolution(i*h, 0), u[i], exactSolution(i*h, 0)-u[i]);
        }

        System.out.println("n=" + n + " h=" + h);
        System.out.println("k=" + k + " tau=" + tau );
        if(tau <= ((h*h)/(2*a)) ) System.out.println("tau < " + ((h*h)/(2*a)));

        double[] res = solve_implicit(u);
    }            */

    public static double exactSolution(double x, double t){return exp((c-a)*t)*sin(x);}
    public static void initSimple(double[] u){for(int i =0; i< u.length; ++i) u[i]=ksi(i*h);}



    public static double leftApproximation2pI_explicit(double u[], double t){
        return (-1)*(( ( alpha/h )/( beta - (alpha/h) ) )*u[1]) + ( fi_0(t)/( beta - (alpha/h) ));
    }
    public static double rightApproximation2pI_explicit(double u[], double t){
        return ((gamma/h)/(delta + gamma/h))*u[u.length - 2] + (fi_l(t)/(delta + gamma/h));
    }

    public static double leftApproximation3pII_explicit(double u[], double t){
        return ((( alpha/(2*h) )/( beta - ((3*alpha)/(2*h)) ) )*(u[2]-4*u[1]) ) + ( fi_0(t)/( beta - ((3*alpha)/(2*h)) ));
    }
    public static double rightApproximation3pII_explicit(double u[], double t){
        return ((gamma/(2*h))/(delta + ((3*gamma)/(2*h))))*(4*u[u.length - 2]-u[u.length-3]) + (fi_l(t)/(delta + (3*gamma)/(2*h)));
    }
    public static double leftApproximation2pII_explicit(double u[], double u_prev[], double t){
        return (1/((2*a*alpha/h) + (h*alpha/tau) - c*h*alpha - beta*(2*a-b*h)))*(u[1]*(2*a*alpha/h)
                + u_prev[0]*(h*alpha/tau) - fi_0(t)*(2*a-b*h) + f(0)*h*alpha);
    }
    public static double rightApproximation2pII_explicit(double u[], double u_prev[], double t){
        return (1/((2*a*gamma/h) + (h*gamma/tau) - c*h*gamma - delta*(2*a-b*h)))*(u[u.length-2]*(2*a*gamma/h)
                + u_prev[u_prev.length-1]*(h*gamma/tau) - fi_l(t)*(2*a-b*h) + f(Xn)*h*alpha);
    }

    public static double[] finiteDifferenceExplicit(double u_prev[], double t){
        double u[] = new double[u_prev.length];
        for(int i=1; i < u_prev.length-1; ++i){
               u[i]=u_prev[i] + (a*tau/(h*h))*(u_prev[i+1]-2*u_prev[i]+u_prev[i-1]) +
                       (b*tau/2*h)*(u_prev[i+1]-u_prev[i-1]) + c*tau*u_prev[i] + tau*f(i*h);
        }
        //u[0]= leftApproximation3pII(u, t);
        //u[u.length-1]=rightApproximation3pII(u, t);
        u[0]= leftApproximation2pII_explicit(u, u_prev, t);
        u[u.length-1]= rightApproximation2pII_explicit(u, u_prev, t);
        return u;
    }

    public static double[] solve_explicit(double[] u){
        print(u, 0);
        for (int i =1; i<= k; ++i){
            u = finiteDifferenceExplicit(u, i * tau);
            print(u, i);
        }
        return u;
    }








    public static void leftApproximation2pI_implicit(double u[][], double d[][], double t){
        u[0][0]=beta-(alpha/h);
        u[0][1]=alpha/h;
        d[0][0]=fi_0(t);
    }
    public static void rightApproximation2pI_implicit(double u[][], double d[][], double t){
        u[n-1][n-2]=(-1)*(gamma/h);
        u[n-1][n-1]=(delta+(gamma/h));
        d[n-1][0]=fi_l(t);
    }

    public static void leftApproximation3pII_implicit(double u[][], double d[][], double t){
        u[0][0]= (beta - ((3*alpha)/(2*h))) + ((alpha/(2*h))*(u[1][0]/u[1][2])) ;
        u[0][1]= (2*alpha/h) + ((u[1][1]/u[1][2])*(alpha/(2*h)));
        d[0][0]= fi_0(t)+(d[1][0]/u[1][2])*(alpha/(2*h));
    }
    public static void rightApproximation3pII_implicit(double u[][], double d[][], double t){
        u[n-1][n-2]= ((-2)*gamma/h) - (gamma/(2*h))*(u[n-2][n-2]/u[n-2][n-3]);
        u[n-1][n-1]=  delta +(3*gamma)/(2*h) - (gamma/(2*h))*(u[n-2][n-1]/u[n-2][n-3]);
        d[n-1][0]= fi_l(t)-(gamma/(2*h))*(d[n-2][0]/u[n-2][n-3]);
    }

    public static double[] solve_implicit(double[] u_prev, double tetta){
        double[][] u = create_Matrix(n, n);
        double[][] d = create_Matrix(n, 1);
        tetta=1;
        print(u_prev, 0);

        for(int j = 1; j<=k; ++j){

            for (int i =1; i<n-1; ++i){
                u[i][i-1]= tetta*tau*((a/(h*h)) - (b/(2*h)));
                System.out.println("tau=" + tau + " c=" + c + " a=" + a + " h=" + h);
                u[i][i]= -1 + tetta*tau*(c - (2*a)/(h*h));
                u[i][i+1] = tetta*tau*((a/(h*h))+(b/(2*h)));
                d[i][0]=(-1)*u_prev[i]+ (tetta-1)*tau*( (a*(u_prev[i+1]-2*u_prev[i]+u_prev[i-1])/(h*h))
                        + (b*(u_prev[i+1]-u_prev[i-1])/(2*h)) + c*u_prev[i] ) +tau*f(i*h);
            }
            leftApproximation3pII_implicit(u, d, j * tau);
            rightApproximation3pII_implicit(u, d, j * tau);


            for(int i = 0; i<n; ++i){
                for(int z = 0; z<n; ++z){
                    System.out.format("%.3f ", u[i][z]);
                }
                System.out.format("    %.3f ", d[i][0]);
                System.out.format("\n");
            }
            System.out.format("\n");

            u_prev = Sweep_method(u, d);

            print(u_prev, j);
        }
        return u_prev;
    }




    public static void print(double[] u, int j){
        System.out.format("    u[%d] | ", j);
        for(int i =0; i<u.length; ++i){
            System.out.format("%.3f ", u[i]);
        }
        System.out.format("\nExact[%d] | ", j);
        for (int i=0; i<u.length; ++i) System.out.format("%.3f ", exactSolution(i*h, j*tau));
        System.out.format("\ndelta[%d] | ", j);
        for (int i=0; i<u.length; ++i) System.out.format("%.3f ", exactSolution(i*h, j*tau)-u[i]);
        System.out.println("\n");
    }

    public static double[][] create_Matrix(int n, int m){
        double[][] a = new double[n][m];
        for(int i=0; i<a.length; ++i){
            a[i] = new double[m];
        }
        return a;
    }

    public static double[] Sweep_method(double[][] A, double[][] B){
        double[][] P=create_Matrix(n, 1);
        double[][] Q=create_Matrix(n, 1);

        P[0][0]= A[0][1] /(-1)*A[0][0];
        Q[0][0]=-B[0][0] /(-1)*A[0][0];

        for(int i=1; i<n; i++){
            if(i<n-1)
                P[i][0]= A[i][i+1] /(-A[i][i]- A[i][i-1]*P[i-1][0]);
            Q[i][0]= (A[i][i-1]*Q[i-1][0] - B[i][0])/(-A[i][i]- A[i][i-1]*P[i-1][0]);
        }

        double[][] X=create_Matrix(n, 1);
        X[n-1][0]= Q[n-1][0];
        for(int i=n-1; i>0; i--){
            X[i-1][0]=P[i-1][0]*X[i][0]	+ Q[i-1][0];
        }

        double[] res = new double[n];
        for (int i = 1; i<n; ++i) res[i]=X[i][0];

        return res;
    }



    public static int getN(){return n;}
    public static int getK(){return k;}
    public double[] getX() {return x;}
    public double[] getRes() {return res;}

    Pattern p = Pattern.compile("(\\d+(?:\\.\\d+))");


    public void setA(String a) {
        Matcher m = p.matcher(a);
        while(m.find()) {
            this.a = Double.parseDouble(m.group(1));
        }
        //Scanner sc = new Scanner(a);
       // this.a = sc.nextDouble();
    }
    public void setB(String b) {
        Matcher m = p.matcher(b);
        while(m.find()) {
            this.b = Double.parseDouble(m.group(1));
        }
    }
    public void setC(String c) {
        Matcher m = p.matcher(c);
        while(m.find()) {
            this.c = Double.parseDouble(m.group(1));
        }
    }
    public void setN(String n) {
        Scanner sc = new Scanner(n);
        this.n = sc.nextInt();
    }
    public void setK(String k) {
        Scanner sc = new Scanner(k);
        this.k = sc.nextInt();
    }

}
