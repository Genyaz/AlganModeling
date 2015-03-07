package com.company;

import java.util.Arrays;

public abstract class Function {
    private static double eps = 1e-6;

    //public abstract int arity();
    public abstract double calculate(double[] x);

    public double[] totalDerivative(double[] x) {
        double[] xn = Arrays.copyOf(x, x.length);
        double y = calculate(x);
        double[] res = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            xn[i] += eps;
            res[i] = (calculate(xn) - y) / eps;
            xn[i] = x[i];
        }
        return res;
    }
}
