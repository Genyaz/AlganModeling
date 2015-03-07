package com.company;

import java.util.Arrays;
import java.util.Random;

/**
 * System is represented by an array of functions that returns 0 on the answer.
 */
public class EquationSystem {

    private final static Random random = RandomHolder.random;
    private final static int SEARCH_SEGMENT_DIVISION = 10;
    private final static int SEARCH_ITERATIONS = 10;
    private final static double GRADIENT_DESCENT_PRECISION = 1e-6;

    /**
     * Finds argument of rough function minimum
     * @param left left search border
     * @param right right search border
     * @param f function
     * @param m search segment division
     * @param n search iterations
     * @return x| f(x) is rough minimum
     */
    public static double superSearch(double left, double right, Function f, int m, int n) {
        double[] arg = new double[1];
        double answer = 1;
        arg[0] = 1;
        double min = f.calculate(arg);
        for (int i = 0; i < m; i++) {
            double l = left + i * (right - left) / m;
            arg[0] = l;
            double cur = f.calculate(arg);
            if (cur < min) {
                min = cur;
                answer = l;
            }
            double r = left + (i + 1) * (right - left) / m;
            arg[0] = r;
            cur = f.calculate(arg);
            if (cur < min) {
                min = cur;
                answer = r;
            }
            double p, q, fp, fq;
            for (int j = 0; j < n; j++) {
                p = (2 * l + r) / 3;
                arg[0] = p;
                fp = f.calculate(arg);
                if (fp < min) {
                    min = fp;
                    answer = p;
                }
                q = (l + 2 * r) / 3;
                arg[0] = q;
                fq = f.calculate(arg);
                if (fq < min) {
                    min = fq;
                    answer = q;
                }
                if (fp < fq) {
                    r = q;
                } else {
                    l = p;
                }
            }
        }
        return answer;
    }

    /**
     * Finds argument of rough function minimum
     * @param f function
     * @param x0 initial argument
     * @param initialStep initial dx
     * @param precision minimum dx
     * @return x| f(x) is rough minimum
     */
    public static double gradientDescent(Function f, double x0, double initialStep, double precision) {
        double[] arg = new double[1];
        arg[0] = x0;
        double derivative = f.totalDerivative(arg)[0];
        double step = initialStep, x = x0, x1, min = f.calculate(arg);
        while (step > precision) {
            if (derivative < 0) {
                x1 = x + step;
            } else {
                x1 = x - step;
            }
            arg[0] = x1;
            double cur = f.calculate(arg);
            if (cur < min) {
                min = cur;
                x = x1;
                arg[0] = x;
                derivative = f.totalDerivative(arg)[0];
            } else {
                step /= 2;
            }
        }
        return x;
    }

    private Function[] functions;
    private int n;

    public EquationSystem(Function[] functions) {
        this.functions = functions;
        n = functions.length;
    }

    /**
     * @param x argument
     * @return discrepancy = sum fi(x)*fi(x), i = 0..n-1
     */
    public double discrepancy(double[] x) {
        double discrepancy = 0, q;
        for (int i = 0; i < n; i++) {
            q = functions[i].calculate(x);
            discrepancy += q * q;
        }
        return discrepancy;
    }

    /**
     * @param x0
     * @param d
     * @param t
     * @return {@link com.company.EquationSystem#discrepancy}(x0 + t * d)
     */
    public double discrepancy(double[] x0, double[] d, double t) {
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = x0[i] + t * d[i];
        }
        return discrepancy(x);
    }

    /**
     * @param x initial function argument
     * @return array[0..n], where array[0..n-1] is dx, and array[n] is discrepancy of fi(x)
     */
    public double[] linearDerivativeSolution(double[] x) {
        double[] b = new double[n], dx;
        double[][] matrix = new double[n][];
        double discrepancy = 0;
        for (int i = 0; i < n; i++) {
            b[i] = -functions[i].calculate(x);
            discrepancy += b[i] * b[i];
            matrix[i] = functions[i].totalDerivative(x);
        }
        Matrix m = new Matrix(matrix);
        dx = Arrays.copyOf(m.gaussMethod(b), n + 1);
        dx[n] = discrepancy;
        return dx;
    }

    /**
     * @param eps allowable discrepancy
     * @param maxIterations maximum iterations count
     * @return argument x, discrepancy(x) < eps
     */
    public double[] newtonMethod(double eps, long maxIterations) {
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = random.nextDouble();
        }
        for (int q = 0; q < maxIterations; q++) {
            double [] dx = linearDerivativeSolution(x);
            if (dx[n] < eps) break;
            for (int i = 0; i < n; i++) {
                x[i] += dx[i];
            }
        }
        return x;
    }

    /**
     * Finds local discrepancy minimum on line x(t) = x + t * direction
     * @param x initial argument
     * @param d line direction
     * @return t| x(t) = x + t * direction is local discrepancy minimum
     */
    public double localMinimum(final double[] x, final double[] d) {
        double cur = discrepancy(x);
        cur = Math.min(cur, discrepancy(x, d, 1));
        double r = 1, dr;
        do {
            r *= 2;
            dr = discrepancy(x, d, r);
            cur = Math.min(cur, dr);
        } while (dr <= cur);
        Function f = new Function() {
            @Override
            public double calculate(double[] arg) {
                return discrepancy(x, d, arg[0]);
            }
        };
        //return superSearch(0, r, f, SEARCH_SEGMENT_DIVISION, SEARCH_ITERATIONS);
        return gradientDescent(f, 1, 0.5, GRADIENT_DESCENT_PRECISION);
    }

    /**
     * @param eps allowable discrepancy
     * @param maxIterations maximum iterations count
     * @return argument x, discrepancy(x) < eps
     */
    public double[] universalMethod(double eps, long maxIterations) {
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = random.nextDouble();
        }
        for (int q = 0; q < maxIterations; q++) {
            double [] dx = linearDerivativeSolution(x);
            if (dx[n] < eps) break;
            double k = localMinimum(x, dx);
            for (int i = 0; i < n; i++) {
                x[i] += k * dx[i];
            }
        }
        return x;
    }
}
