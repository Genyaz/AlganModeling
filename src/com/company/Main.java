package com.company;

import java.io.*;
import java.util.*;

public class Main {

    // We allow answer pressures in [-ALLOWED_DISCREPANCY; ATMOSPHERIC_PRESSURE + ALLOWED_DISCREPANCY]
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
        // x[i] = Pe(chemicalAgent[i])
        Function[] functions = new Function[5];
        // 2 HCl + 2 Al = 2 AlCl + H2
        // Pe(HCl)^2 = K1 * Pe(AlCl)^2 * Pe(H2)
        functions[0] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] - k1 * x[1] * x[1] * x[4];
            }
        };
        // 2 HCl + Al = AlCl2 + H2
        // Pe(HCl) ^ 2 = K2 * Pe(AlCl2) * Pe(H2)
        functions[1] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] - k2 * x[2] * x[4];
            }
        };
        // 6 HCl + 2 Al = 2 AlCl3 + 3 H2
        // Pe(HCl)^6 = K3 * Pe(AlCl3)^2 * Pe(H2)^3
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
        // G(H) = G(HCl) + 2 * G(H2) = 0
        // D(HCl) * (Pg(HCl) - Pe(HCl)) + 2 * D(H2) * (Pg(H2) - Pe(H2)) = 0
        functions[3] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[0] * (p[0] - x[0]) + 2 * d[4] * (p[4] - x[4]);
            }
        };
        // G(Cl) = G(HCl) + G(AlCl) + 2 * G(AlCl2) + 3 * G(AlCl3) = 0
        // D(HCl) * (Pg(HCl) - Pe(HCl)) + D(AlCl) * (Pg(AlCl) - Pe(AlCl)) + 2 * D(AlCl2) * (Pg(AlCl2) - Pe(AlCl2))
        // + 3 * D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) = 0
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
                correct &= (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.ATMOSPHERIC_PRESSURE + ALLOWED_DISCREPANCY);
            }
        }
        System.out.println("T = " + T);
        out.println("T = " + T);
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
        // x[i] = Pe(chemicalAgent[i])
        Function[] functions = new Function[5];
        // 2 HCl + 2 Ga = 2 GaCl + H2
        // Pe(HCl)^2 = K4 * Pe(GaCl)^2 * Pe(H2)
        functions[0] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] - k4 * x[1] * x[1] * x[4];
            }
        };
        // 2 HCl + Ga = GaCl2 + H2
        // Pe(HCl)^2 = K5 * Pe(GaCl2) * Pe(H2)
        functions[1] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] - k5 * x[2] * x[4];
            }
        };
        // 6 HCl + 2 Ga = 2 GaCl3 + 3 H2
        // Pe(HCl)^6 = K6 * Pe(GaCl3)^2 * Pe(H2)^3
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
        // G(H) = G(HCl) + 2 * G(H2) = 0
        // D(HCl) * (Pg(HCl) - Pe(HCl)) + 2 * D(H2) * (Pg(H2) - Pe(H2)) = 0
        functions[3] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[0] * (p[0] - x[0]) + 2 * d[4] * (p[4] - x[4]);
            }
        };
        // G(Cl) = G(HCl) + G(GaCl) + 2 * G(GaCl2) + 3 * G(GaCl3) = 0
        // D(HCl) * (Pg(HCl) - Pe(HCl)) + D(GaCl) * (Pg(GaCl) - Pe(GaCl)) + 2 * D(GaCl2) * (Pg(GaCl2) - Pe(GaCl2))
        // + 3 * D(GaCl3) * (Pg(GaCl3) - Pe(GaCl3)) = 0
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
                correct &= (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.ATMOSPHERIC_PRESSURE + ALLOWED_DISCREPANCY);
            }
        }
        System.out.println("T = " + T);
        out.println("T = " + T);
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
        // x[i] = Pe(chemicalAgent[i]), i = 0..4
        // x[5] = x = G(AlCl3) / (G(AlCl3) + G(GaCl))
        Function[] functions = new Function[6];
        // AlCl3 + NH3 = AlN + 3 HCl
        // Pe(AlCl3) * Pe(NH3) = K9 * x * Pe(HCl)^3
        functions[0] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[3] * x[2] - k9 * x[5] * x[0] * x[0] * x[0];
            }
        };
        // GaCl + NH3 = GaN + HCl + H2
        // Pe(GaCl) * Pe(NH3) = K10 * (1 - x) * Pe(HCl) * Pe(H2)
        functions[1] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[1] * x[2] - k10 * (1 - x[5]) * x[0] * x[4];
            }
        };
        // G(H) = G(HCl) + 2 * G(H2) + 3 * G(NH3) = 0
        // D(HCl) * (Pg(HCl) - Pe(HCl)) + 2 * D(H2) * (Pg(H2) - Pe(H2)) + 3 * D(NH3) * (Pg(NH3) - Pe(NH3))
        functions[2] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[0] * (p[0] - x[0]) + 2 * d[4] * (p[4] - x[4]) + 3 * d[2] * (p[2] - x[2]);
            }
        };
        // G(Cl) = 3 * G(AlCl3) + G(GaCl) + G(HCl) = 0
        // 3 * D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) + D(GaCl) * (Pg(GaCl) - Pe(GaCl)) + D(HCl) * (Pg(HCl) - Pe(HCl)) = 0
        functions[3] = new Function() {
            @Override
            public double calculate(double[] x) {
                return 3 * d[3] * (p[3] - x[3]) + d[1] * (p[1] - x[1]) + d[0] * (p[0] - x[0]);
            }
        };
        // G(Al) + G(Ga) = G(AlCl3) + G(GaCl) = G(NH3) = G(N)
        // D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) + D(GaCl) * (Pg(GaCl) - Pe(GaCl)) = D(NH3) * (Pg(NH3) - Pe(NH3))
        functions[4] = new Function() {
            @Override
            public double calculate(double[] x) {
                return d[3] * (p[3] - x[3]) + d[1] * (p[1] - x[1]) - d[2] * (p[2] - x[2]);
            }
        };
        // G(AlCl3) = x * (G(AlCl3) + G(GaCl))
        // D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) = x * (D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) + D(GaCl) * (Pg(GaCl) - Pe(GaCl)))
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
                correct &= (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.ATMOSPHERIC_PRESSURE + ALLOWED_DISCREPANCY);
            }
            correct &= (x[5] >= 0 && x[5] <= 1);
        }
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

    public static Map<String, List<Double>> parseFile(String fileName) throws IOException {
        Map<String, List<Double>> result = new HashMap<String, List<Double>>();
        BufferedReader in = new BufferedReader(new FileReader(fileName));
        String s;
        while ((s = in.readLine()) != null) {
            String[] strings = s.split(" ");
            if (result.get(strings[0]) == null) {
                result.put(strings[0], new ArrayList<Double>());
            }
            result.get(strings[0]).add(Double.parseDouble(strings[2]));
        }
        in.close();
        return result;
    }

    public static void main(String[] args) throws IOException {
        PrintWriter out;
        Map<String, Double> pressure = new HashMap<String, Double>();
        pressure.put("HCl", 10000d);
        pressure.put("N2", 90000d);
        pressure.put("AlCl", 0d);
        pressure.put("AlCl2", 0d);
        pressure.put("AlCl3", 0d);
        pressure.put("GaCl", 0d);
        pressure.put("GaCl2", 0d);
        pressure.put("GaCl3", 0d);
        pressure.put("H2", 0d);

        // Task 1
        /*
        out = new PrintWriter("task1.out");
        for (int i = 35; i <= 65; i++) {
            double T = 10 * i + 273;
            solveAlClx(pressure, T, 0.01, out);
        }
        out.close();/**/

        // Task 2
        /*
        out = new PrintWriter("task2.out");
        for (int i = 65; i <= 95; i++) {
            double T = 10 * i + 273;
            solveGaClx(pressure, T, 0.01, out);
        }
        out.close();/**/

        // Task 3
        /*
        out = new PrintWriter("task3_pure_N2.out");
        pressure.put("NH3", 1500d);
        pressure.put("HCl", 0d);
        System.out.println("Pure N2");
        pressure.put("N2", 98470d);
        pressure.put("H2", 0d);
        for (int i = 0; i <= 30; i++) {
            pressure.put("AlCl3", (double)i);
            pressure.put("GaCl", (double)(30 - i));
            solveAlGaN(pressure, 1100 + 273, 0.01, out);
        }
        out.close();
        out = new PrintWriter("task3_N2_H2.out");
        System.out.println("H2/N2 = 1/9");
        pressure.put("N2", 88623d);
        pressure.put("H2", 9847d);
        for (int i = 0; i <= 30; i++) {
            pressure.put("AlCl3", (double)i);
            pressure.put("GaCl", (double)(30 - i));
            solveAlGaN(pressure, 1100 + 273, 0.01, out);
        }
        out.close();
        /**/
        // Parsing results
        Locale locale = new Locale("ru");
        Map<String, List<Double>> resMap = parseFile("task1.out");
        List<Double> t = resMap.get("T");
        String[] valueNames = new String[] {"G(AlCl)", "G(AlCl2)", "G(AlCl3)", "Ve(Al)"};
        for (String valueName: valueNames) {
            out = new PrintWriter("task1_" + valueName + ".out");
            List<Double> values = resMap.get(valueName);
            for (int i = 0; i < t.size(); i++) {
                out.printf(locale, "%f\t%f\n", (1 / t.get(i)), Math.log(Math.abs(values.get(i))));
            }
            out.close();
        }
        resMap = parseFile("task2.out");
        t = resMap.get("T");
        valueNames = new String[] {"G(GaCl)", "G(GaCl2)", "G(GaCl3)", "Ve(Ga)"};
        for (String valueName: valueNames) {
            out = new PrintWriter("task2_" + valueName + ".out");
            List<Double> values = resMap.get(valueName);
            for (int i = 0; i < t.size(); i++) {
                out.printf(locale, "%f\t%f\n", (1 / t.get(i)), Math.log(Math.abs(values.get(i))));
            }
            out.close();
        }
        resMap = parseFile("task3_pure_N2.out");
        List<Double> PgAlCl3 = resMap.get("Pg(AlCl3)");
        valueNames = new String[] {"G(AlCl3)", "G(GaCl)", "Vg(AlGaN)", "x"};
        for (String valueName: valueNames) {
            out = new PrintWriter("task3_" + valueName + "N2.out");
            List<Double> values = resMap.get(valueName);
            for (int i = 0; i < PgAlCl3.size(); i++) {
                out.printf(locale, "%f\t%e\n", (PgAlCl3.get(i) / 30), values.get(i));
            }
            out.close();
        }
        resMap = parseFile("task3_N2_H2.out");
        PgAlCl3 = resMap.get("Pg(AlCl3)");
        valueNames = new String[] {"G(AlCl3)", "G(GaCl)", "Vg(AlGaN)", "x"};
        for (String valueName: valueNames) {
            out = new PrintWriter("task3_" + valueName + "N2_H2.out");
            List<Double> values = resMap.get(valueName);
            for (int i = 0; i < PgAlCl3.size(); i++) {
                out.printf(locale, "%f\t%e\n", (PgAlCl3.get(i) / 30), values.get(i));
            }
            out.close();
        }
    }
}
