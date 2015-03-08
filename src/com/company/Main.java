package com.company;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

public class Main {

    public static final double ALLOWED_DISCREPANCY = 1000;

    public static void solveAlClx(Map<String, Double> pressure, double T, double delta, PrintWriter out) {
        final String[] chemicalAgent = new String[5];
        chemicalAgent[0] = "HCl";
        chemicalAgent[1] = "AlCl";
        chemicalAgent[2] = "AlCl2";
        chemicalAgent[3] = "AlCl3";
        chemicalAgent[4] = "H2";
        final double k1 = DataHolder.getK(T, 1);
        final double k2 = DataHolder.getK(T, 2);
        final double k3 = DataHolder.getK(T, 3);
        Function[] functions = new Function[5];
        functions[0] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] - k1 * x[1] * x[1] * x[4];
            }
        };
        functions[1] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] - k2 * x[2] * x[4];
            }
        };
        functions[2] = new Function() {
            @Override
            public double calculate(double[] x) {
                return Math.pow(x[0], 6) - k3 * x[3] * x[3] * x[4] * x[4] * x[4];
            }
        };
        final double[] p = new double[5];
        final double[] d = new double[5];
        for (int i = 0; i < 5; i++) {
            p[i] = pressure.get(chemicalAgent[i]);
            d[i] = DataHolder.getD(chemicalAgent[i], T);
        }
        functions[3] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[0] * (p[0] - x[0]) + 2 * d[4] * (p[4] - x[4]);
            }
        };
        functions[4] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[1] * (p[1] - x[1]) + 2 * d[2] * (p[2] - x[2]) + 3 * d[3] * (p[3] - x[3]) + d[0] * (p[0] - x[0]);
            }
        };
        EquationSystem equationSystem = new EquationSystem(functions);
        boolean correct = false;
        double[] x = null;
        while (!correct) {
            x = equationSystem.universalMethod(1e-12, 1000000);
            correct = true;
            for (int i = 0; i < 5; i++) {
                correct &= (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.PRESSURE + ALLOWED_DISCREPANCY);
            }
        }
        System.out.println("T = " + T + "K");
        out.println("T = " + T + "K");
        for (int i = 0; i < 5; i++) {
            System.out.println("Pe(" + chemicalAgent[i] + ") = " + x[i]);
            out.println("Pe(" + chemicalAgent[i] + ") = " + x[i]);
        }
        double[] g = new double[5];
        for (int i = 0; i < 5; i++) {
            g[i] = d[i] * (p[i] - x[i]) / (8314 * T * delta);
            System.out.println("G(" + chemicalAgent[i] + ") = " + g[i]);
            out.println("G(" + chemicalAgent[i] + ") = " + g[i]);
        }
        double v = (g[1] + g[2] + g[3]) * (DataHolder.getDouble("mu", "Al") / DataHolder.getDouble("density", "Al")) * 1000000000;
        System.out.println("Ve(Al) = " + v);
        out.println("Ve(Al) = " + v);
    }

    public static void solveGaClx(Map<String, Double> pressure, double T, double delta, PrintWriter out) {
        final String[] chemicalAgent = new String[5];
        chemicalAgent[0] = "HCl";
        chemicalAgent[1] = "GaCl";
        chemicalAgent[2] = "GaCl2";
        chemicalAgent[3] = "GaCl3";
        chemicalAgent[4] = "H2";
        final double k4 = DataHolder.getK(T, 4);
        final double k5 = DataHolder.getK(T, 5);
        final double k6 = DataHolder.getK(T, 6);
        Function[] functions = new Function[5];
        functions[0] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] - k4 * x[1] * x[1] * x[4];
            }
        };
        functions[1] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] - k5 * x[2] * x[4];
            }
        };
        functions[2] = new Function() {
            @Override
            public double calculate(double[] x) {
                return Math.pow(x[0], 6) - k6 * x[3] * x[3] * x[4] * x[4] * x[4];
            }
        };
        final double[] p = new double[5];
        final double[] d = new double[5];
        for (int i = 0; i < 5; i++) {
            p[i] = pressure.get(chemicalAgent[i]);
            d[i] = DataHolder.getD(chemicalAgent[i], T);
        }
        functions[3] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[0] * (p[0] - x[0]) + 2 * d[4] * (p[4] - x[4]);
            }
        };
        functions[4] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[1] * (p[1] - x[1]) + 2 * d[2] * (p[2] - x[2]) + 3 * d[3] * (p[3] - x[3]) + d[0] * (p[0] - x[0]);
            }
        };
        EquationSystem equationSystem = new EquationSystem(functions);
        boolean correct = false;
        double[] x = null;
        while (!correct) {
            x = equationSystem.universalMethod(1e-12, 1000000);
            correct = true;
            for (int i = 0; i < 5; i++) {
                correct &= (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.PRESSURE + ALLOWED_DISCREPANCY);
            }
        };
        System.out.println("T = " + T + "K");
        out.println("T = " + T + "K");
        for (int i = 0; i < 5; i++) {
            System.out.println("Pe(" + chemicalAgent[i] + ") = " + x[i]);
            out.println("Pe(" + chemicalAgent[i] + ") = " + x[i]);
        }
        double[] g = new double[5];
        for (int i = 0; i < 5; i++) {
            g[i] = d[i] * (p[i] - x[i]) / (8314 * T * delta);
            System.out.println("G(" + chemicalAgent[i] + ") = " + g[i]);
            out.println("G(" + chemicalAgent[i] + ") = " + g[i]);
        }
        double v = (g[1] + g[2] + g[3]) * (DataHolder.getDouble("mu", "Ga") / DataHolder.getDouble("density", "Ga")) * 1000000000;
        System.out.println("Ve(Ga) = " + v);
        out.println("Ve(Ga) = " + v);
    }

    public static void solveAlGaN(Map<String, Double> pressure, double T, double delta, PrintWriter out) {
        final String[] chemicalAgent = new String[5];
        chemicalAgent[0] = "HCl";
        chemicalAgent[1] = "GaCl";
        chemicalAgent[2] = "NH3";
        chemicalAgent[3] = "AlCl3";
        chemicalAgent[4] = "H2";
        final double[] p = new double[5];
        final double[] d = new double[5];
        for (int i = 0; i < 5; i++) {
            p[i] = pressure.get(chemicalAgent[i]);
            d[i] = DataHolder.getD(chemicalAgent[i], T);
        }
        final double k9 = DataHolder.getK(T, 9);
        final double k10 = DataHolder.getK(T, 10);
        Function[] functions = new Function[6];
        functions[0] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[3] * x[2] - k9 * x[5] * x[0] * x[0] * x[0] * x[0];
            }
        };
        functions[1] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[1] * x[2] - k10 * (1 - x[5]) * x[0] * x[4];
            }
        };
        functions[2] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[0] * (p[0] - x[0]) + 2 * d[4] * (p[4] - x[4]) + 3 * d[2] * (p[2] - x[2]);
            }
        };
        functions[3] = new Function() {
            @Override
            public double calculate(double[] x) {
                return 3 * d[3] * (p[3] - x[3]) + d[1] * (p[1] - x[1]) + d[0] * (p[0] - x[0]);
            }
        };
        functions[4] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[3] * (p[3] - x[3]) + d[1] * (p[1] - x[1]) - d[2] * (p[2] - x[2]);
            }
        };
        functions[5] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[3] * (p[3] - x[3]) - x[5] * (d[1] * (p[1] - x[1]) + d[3] * (p[3] - x[3]));
            }
        };
        EquationSystem equationSystem = new EquationSystem(functions);
        boolean correct = false;
        double[] x = null;
        while (!correct) {
            x = equationSystem.universalMethod(1e-12, 1000000);
            correct = true;
            for (int i = 0; i < 5; i++) {
                correct &= (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.PRESSURE + ALLOWED_DISCREPANCY);
            }
            correct &= (x[5] >= 0 && x[5] <= 1);
        };
        System.out.println("Pg(AlCl3) = " + pressure.get("AlCl3"));
        out.println("Pg(AlCl3) = " + pressure.get("AlCl3"));
        for (int i = 0; i < 5; i++) {
            System.out.println("Pe(" + chemicalAgent[i] + ") = " + x[i]);
            out.println("Pe(" + chemicalAgent[i] + ") = " + x[i]);
        }
        System.out.println("x = " + x[5]);
        out.println("x = " + x[5]);
        double[] g = new double[5];
        for (int i = 0; i < 5; i++) {
            g[i] = d[i] * (p[i] - x[i]) / (8314 * T * delta);
            System.out.println("G(" + chemicalAgent[i] + ") = " + g[i]);
            out.println("G(" + chemicalAgent[i] + ") = " + g[i]);
        }
        double v = (g[3] * (DataHolder.getDouble("mu", "AlN") / DataHolder.getDouble("density", "AlN"))
            + g[1] * (DataHolder.getDouble("mu", "GaN") / DataHolder.getDouble("density", "GaN"))) * 1000000000;
        System.out.println("Vg(AlGaN) = " + v);
        out.println("Vg(AlGaN) = " + v);
    }

    public static void main(String[] args) throws FileNotFoundException {
        Map<String, Double> pressure = new HashMap<String, Double>();
        PrintWriter out;
        /*pressure.put("HCl", 10000d);
        pressure.put("N2", 90000d);
        pressure.put("AlCl", 0d);
        pressure.put("AlCl2", 0d);
        pressure.put("AlCl3", 0d);
        pressure.put("GaCl", 0d);
        pressure.put("GaCl2", 0d);
        pressure.put("GaCl3", 0d);
        pressure.put("H2", 0d);
        out = new PrintWriter("task1.out");
        for (int i = 35; i <= 65; i++) {
            double T = 10 * i + 273;
            solveAlClx(pressure, T, 0.01, out);
        }
        out.close();
        out = new PrintWriter("task2.out");
        for (int i = 65; i <= 95; i++) {
            double T = 10 * i + 273;
            solveGaClx(pressure, T, 0.01, out);
        }
        out.close(); */
        out = new PrintWriter("task3.out");
        pressure.put("NH3", 1500d);
        pressure.put("HCl", 0d);
        for (int q = 0; q < 2; q++) {
            if (q == 0) {
                System.out.println("Pure N2");
                out.println("Pure N2");
                pressure.put("N2", 98470d);
                pressure.put("H2", 0d);
            } else {
                System.out.println("H2/N2 = 1/9");
                out.println("H2/N2 = 1/9");
                pressure.put("N2", 88623d);
                pressure.put("H2", 9847d);
            }
            for (int i = 0; i <= 30; i++) {
                pressure.put("AlCl3", (double)i);
                pressure.put("GaCl", (double)(30 - i));
                solveAlGaN(pressure, 1100 + 273, 0.01, out);
            }
        }

        out.close();
    }
}
