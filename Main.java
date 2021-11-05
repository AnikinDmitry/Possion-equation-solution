import java.io.FileWriter;
import java.io.IOException;
import java.util.Formatter;

public class Main {
    public static void main(String[] args) throws IOException {
        FileWriter out = new FileWriter("out.txt");
        Formatter s = new Formatter();

        s.format("First scheme:\n");
        for (int n = 5; n <= 80; n *= 2) {
            s.format("%d\t%d\t%.15f\t%.15f\n", n, n, schemeError(1, n, n),
                    schemeConverge(1, n, n));
        }
        s.format("Second scheme:\n");
        for (int n = 5; n <= 80; n *= 2) {
            s.format("%d\t%d\t%.15f\t%.15f\n", n, n, schemeError(2, n, n),
                    schemeConverge(2, n, n));
        }
        out.write(String.valueOf(s));
        out.close();
    }

    // Exact solution
    public static double u(double x, double y) {
        return Math.cos(x + y) * Math.sin(x * y);
    }

    // y-differential of exact solution
    public static double u_y(double x, double y) {
        return x * Math.cos(x + y) * Math.cos(x * y) - Math.sin(x + y) * Math.sin(x * y);
    }

    // Right part
    public static double f(double x, double y) {
        return -Math.cos(x + y) * Math.sin(x * y) * (2 + x * x + y * y) -
                2 * Math.sin(x + y) * Math.cos(x * y) * (x + y);
    }

    // A - square matrix
    public static double norm(double[][] A) {
        double n = 0;

        for (double[] A_i : A) {
            for (int j = 0; j < A.length; j++) {
                if (Math.abs(A_i[j]) > n) {
                    n = Math.abs(A_i[j]);
                }
            }
        }

        return n;
    }

    // Ax = f, A - square three-diagonal matrix
    public static double[] sweepMethod(double[][] A, double[] f) {
        int n = A.length;
        double[] a = new double[n];
        double[] b = new double[n];

        a[0] = -A[0][1] / A[0][0];
        b[0] = f[0] / A[0][0];
        for (int i = 1; i < n - 1; i++) {
            a[i] = -A[i][i + 1] / (A[i][i] + a[i - 1] * A[i][i - 1]);
            b[i] = (f[i] - b[i - 1] * A[i][i - 1]) /
                    (A[i][i] + a[i - 1] * A[i][i - 1]);
        }

        double[] x = new double[n];
        x[n - 1] = (f[n - 1] - b[n - 2] * A[n - 1][n - 2]) /
                (A[n - 1][n - 1] + a[n - 2] * A[n - 1][n - 2]);
        for (int i = n - 2; i >= 0; i--) x[i] = a[i] * x[i + 1] + b[i];

        return x;
    }

    // Longitudinal-transverse sweep scheme
    public static double[][] firstScheme(int N_x, int N_y) {
        double e = 0.0001; // stop error
        double a = Math.PI / 2, b = 3 * Math.PI / 2; // c = a, d = b
        double h_x = (b - a) / N_x, h_y = (b - a) / N_y;
        double t = 0.025;
        double[][] A = new double[N_x + 1][N_x + 1]; // matrix for first sweep
        double[] g = new double[N_x + 1]; // right part for first sweep
        double[][] u = new double[N_x + 1][N_y + 1]; // u^n, u^(n+1) and first iteration
        double[][] u_1 = new double[N_x + 1][N_y + 1]; // u^(n+0.5)
        double[][] B = new double[N_x - 1][N_y - 1]; // matrix for norm and error
        double[][] C = new double[N_y + 1][N_y + 1]; // matrix for second sweep
        double[][] C_clone = new double[N_y + 1][N_y + 1];
        double[] h = new double[N_y + 1]; // right part for first sweep

        A[0][0] = 1;
        A[N_x][N_x] = 1;
        for (int m = 1; m < N_x; m++) {
            A[m][m - 1] = 1 / (h_x * h_x);
            A[m][m] = -2 * (1 / (h_x * h_x) + 1 / t);
            A[m][m + 1] = 1 / (h_x * h_x);
        }

        C[0][0] = -3 / (2 * h_y);
        C[0][1] = 2 / h_y;
        C[0][2] = -1 / (2 * h_y);
        C[N_y][N_y - 2] = 1 / (2 * h_y);
        C[N_y][N_y - 1] = -2 / h_y;
        C[N_y][N_y] = 3 / (2 * h_y);
        for (int j = 1; j < N_y; j++) {
            C[j][j - 1] = 1 / (h_y * h_y);
            C[j][j] = -2 * (1 / (h_y * h_y) + 1 / t);
            C[j][j + 1] = 1 / (h_y * h_y);
        }

        u_1[0][0] = u(a, a);
        u_1[0][N_y] = u(a, b);
        u_1[N_x][0] = u(b, a);
        u_1[N_x][N_y] = u(b, b);

        do {
            for (int j = 1; j < N_y; j++) {
                for (int m = 1; m < N_x; m++) {
                    g[m] = f(a + m * h_x, a + j * h_y) + 2 * (1 / (h_y * h_y) - 1 / t) *
                            u[m][j] - (u[m][j + 1] + u[m][j - 1]) / (h_y * h_y);
                }
                g[0] = u(a, a + j * h_y);
                g[N_x] = u(b, a + j * h_y);

                double[] d = sweepMethod(A, g);
                for (int m = 0; m <= N_x; m++) {
                    u_1[m][j] = d[m];
                }
            }

            for (int m = 1; m < N_x; m++) {
                u_1[m][0] = (2 * h_y * u_y(a + m * h_x, a) - 4 * u_1[m][1] + u_1[m][2]) / -3;
                u_1[m][N_y] = (2 * h_y * u_y(a + m * h_x, b) +
                        4 * u_1[m][N_y - 1] - u_1[m][N_y - 2]) / 3;
            }

            for (int m = 1; m < N_x; m++) {
                for (int i = 0; i <= N_y; i++) {
                    for (int j = 0; j <= N_y; j++) {
                        C_clone[i][j] = C[i][j];
                    }
                }

                for (int j = 1; j < N_y; j++) {
                    h[j] = f(a + m * h_x, a + j * h_y) -
                            (u_1[m - 1][j] + u_1[m + 1][j]) / (h_x * h_x) +
                            2 * u_1[m][j] * (1 / (h_x * h_x) - 1 / t);
                }
                h[0] = u_y(a + m * h_x, a);
                h[N_y] = u_y(a + m * h_x, b);

                double coef_0 = C_clone[0][2] / C_clone[1][2];
                double coef_N_y = C_clone[N_y][N_y - 2] / C_clone[N_y - 1][N_y - 2];
                for (int i = 0; i <= N_y; i++) {
                    C_clone[0][i] -= C_clone[1][i] * coef_0;
                    C_clone[N_y][i] -= C_clone[N_y - 1][i] * coef_N_y;
                }
                h[0] -= h[1] * coef_0;
                h[N_y] -= h[N_y - 1] * coef_N_y;

                double[] d = sweepMethod(C_clone, h);
                for (int j = 0; j <= N_y; j++) {
                    u[m][j] = d[j];
                }
            }

            for (int j = 0; j <= N_y; j++) {
                u[0][j] = u(a, a + j * h_y);
                u[N_x][j] = u(b, a + j * h_y);
            }

            for (int m = 1; m < N_x; m++) {
                for (int j = 1; j < N_y; j++) {
                    B[m - 1][j - 1] = (u[m - 1][j] - 2 * u[m][j] + u[m + 1][j]) / (h_x * h_x) +
                            (u[m][j - 1] - 2 * u[m][j] + u[m][j + 1]) / (h_y * h_y) -
                            f(a + m * h_x, a + j * h_y);
                }
            }
        } while (norm(B) > e);

        return u;
    }

    // Alternately-triangular method
    public static double[][] secondScheme(int N_x, int N_y) {
        double e = 0.0001; // stop error
        double a = Math.PI / 2, b = 3 * Math.PI / 2; // c = a, d = b
        double h_x = (b - a) / N_x, h_y = (b - a) / N_y;
        double t = 0.025;
        double[][] u = new double[N_x + 1][N_y + 1]; // u^n and first iteration
        double[][] u_1 = new double[N_x + 1][N_y + 1]; // u^(n+0.5)
        double[][] u_2 = new double[N_x + 1][N_y + 1]; // u^(n+1)
        double[][] B = new double[N_x - 1][N_y - 1]; // matrix for norm and error
        double coef = 1 + t / (h_x * h_x) + t / (h_y * h_y);

        for (int j = 0; j <= N_y; j++) {
            u_2[0][j] = u(a, a + j * h_y);
            u_2[N_x][j] = u(b, a + j * h_y);
        }

        do {
            for (int m = 1; m < N_x; m++) {
                double alf = coef;
                double bet = -t / (h_x * h_x);
                double gam = -t / (h_y * h_y);
                double alf_1 = -3 / (2 * h_y);
                double bet_1 = 2 / h_y;
                double gam_1 = -1 / (2 * h_y);
                double g_m1 = (u[m - 1][1] - 2 * u[m][1] + u[m + 1][1]) / (h_x * h_x) +
                        (u[m][0] - 2 * u[m][1] + u[m][2]) / (h_y * h_y) -
                        f(a + m * h_x, a + h_y) - bet * u_1[m - 1][1];
                double g_m2 = (u[m - 1][2] - 2 * u[m][2] + u[m + 1][2]) / (h_x * h_x) +
                        (u[m][1] - 2 * u[m][2] + u[m][3]) / (h_y * h_y) -
                        f(a + m * h_x, a + 2 * h_y) - bet * u_1[m - 1][2];

                u_1[m][0] = ((gam_1 * gam / (alf * alf) - bet_1 / alf) *
                        g_m1 - (gam_1 / alf) * g_m2) / (alf_1 - bet_1 * gam / alf + gam_1 * gam *
                        gam / (alf * alf));
                u_1[m][1] = (g_m1 - gam * u_1[m][0]) / alf;
                u_1[m][2] = (g_m2 - gam * u_1[m][1]) / alf;

                for (int j = 3; j < N_y; j++) {
                    u_1[m][j] = ((u[m - 1][j] - 2 * u[m][j] + u[m + 1][j]) / (h_x * h_x) +
                            (u[m][j - 1] - 2 * u[m][j] + u[m][j + 1]) / (h_y * h_y) -
                            f(a + m * h_x, a + j * h_y) + t * u_1[m - 1][j] / (h_x * h_x) +
                            t * u_1[m][j - 1] / (h_y * h_y)) / coef;
                }

                u_1[m][N_y] = (4 * u_1[m][N_y - 1] - u_1[m][N_y - 2]) / 3;
            }

            for (int m = N_x - 1; m > 0 ; m--) {
                double alf = coef / t;
                double bet = -1 / (h_x * h_x);
                double gam = -1 / (h_y * h_y);
                double alf_1 = 3 / (2 * h_y);
                double bet_1 = -2 / h_y;
                double gam_1 = 1 / (2 * h_y);
                double g_mNy_1 = u_1[m][N_y - 1] + bet * u[m + 1][N_y - 1] +
                        gam * u[m][N_y] + alf * u[m][N_y - 1] - bet * u_2[m + 1][N_y - 1];
                double g_mNy_2 = u_1[m][N_y - 2] + bet * u[m + 1][N_y - 2] +
                        gam * u[m][N_y - 1] + alf * u[m][N_y - 2] - bet * u_2[m + 1][N_y - 2];;

                u_2[m][N_y] = (u_y(a + m * h_x, b) + (gam_1 * gam / (alf * alf) - bet_1 / alf) *
                        g_mNy_1 - (gam_1 / alf) * g_mNy_2) / (alf_1 - bet_1 * gam / alf + gam_1 * gam *
                        gam / (alf * alf));
                u_2[m][N_y - 1] = (g_mNy_1 - gam * u_2[m][N_y]) / alf;
                u_2[m][N_y - 2] = (g_mNy_2 - gam * u_2[m][N_y - 1]) / alf;

                for (int j = N_y - 3; j > 0; j--) {
                    u_2[m][j] = (u_1[m][j] + bet * (u[m + 1][j] - u_2[m + 1][j]) +
                            gam * (u[m][j + 1] - u_2[m][j + 1]) + alf * u[m][j]) / alf;
                }

                u_2[m][0] = (2 * h_y * u_y(a + m * h_x, a) - 4 * u_2[m][1] + u_2[m][2]) / -3;
            }

            for (int m = 0; m <= N_x; m++) {
                for (int j = 0; j <= N_y; j++) {
                    u[m][j] = u_2[m][j];
                }
            }

            for (int m = 1; m < N_x; m++) {
                for (int j = 1; j < N_y; j++) {
                    B[m - 1][j - 1] = (u[m - 1][j] - 2 * u[m][j] + u[m + 1][j]) / (h_x * h_x) +
                            (u[m][j - 1] - 2 * u[m][j] + u[m][j + 1]) / (h_y * h_y) -
                            f(a + m * h_x, a + j * h_y);
                }
            }
        } while (norm(B) > e);

        return u;
    }

    public static double schemeError(int sn, int N_x, int N_y) {
        double e1 = 0, e2 = 0.0001;
        double a = Math.PI / 2, b = 3 * Math.PI / 2; // c = a, d = b
        double h_x = (b - a) / N_x, h_y = (b - a) / N_y;
        double[][] u;

        if (sn == 1) u = firstScheme(N_x, N_y);
        else u = secondScheme(N_x, N_y);

        for (int m = 0; m <= N_x; m++) {
            for (int j = 0; j <= N_y; j++) {
                if (Math.abs(u[m][j] - u(a + m * h_x, a + j * h_y)) > e1) {
                    e1 = Math.abs(u[m][j] - u(a + m * h_x, a + j * h_y));
                }

                if (Math.abs(u(a + m * h_x, a + j * h_y)) > e2) {
                    e2 = Math.abs(u(a + m * h_x, a + j * h_y));
                }
            }
        }

        return e1 / e2;
    }

    public static double schemeConverge(int sn, int N_x, int N_y) {
        return Math.log(schemeError(sn, N_x, N_y) / schemeError(sn, 2 * N_x, 2 * N_y)) /
                Math.log(2);
    }
}