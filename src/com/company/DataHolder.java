package com.company;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class DataHolder {
    private final static double PRESSURE = 100000;
    private final static double R = 8.314;
    private final static Map<String, Map<String, String>> data = new HashMap<String, Map<String, String>>();
    private static boolean hasData = false;

    private static void readData() throws IOException {
        BufferedReader in = new BufferedReader(new FileReader("Bank_TD_Fragment.dat"));
        for (int i = 0; i < 6; i++) in.readLine();
        StringTokenizer tokens = new StringTokenizer(in.readLine());
        List<String> parameters = new ArrayList<String>();
        tokens.nextToken();
        while (tokens.hasMoreTokens()) {
            String parameter = tokens.nextToken().toUpperCase();
            parameters.add(parameter);
            data.put(parameter, new HashMap<String, String>());
        }
        in.readLine();
        tokens = new StringTokenizer(in.readLine());
        while (tokens.hasMoreTokens()) {
            String chemicalAgent = tokens.nextToken().toUpperCase();
            for (String parameter : parameters) {
                data.get(parameter).put(chemicalAgent, tokens.nextToken());
            }
            tokens = new StringTokenizer(in.readLine());
        }
        Map<String, String> densityMap = new HashMap<String, String>();
        densityMap.put("AL(S)", "2690");
        densityMap.put("GA(L)", "5900");
        densityMap.put("ALN(S)", "3200");
        densityMap.put("GAN(S)", "6150");
        data.put("DENSITY", densityMap);
    }

    public static String getData(String parameter, String chemicalAgent) {
        if (!hasData) {
            try {
                readData();
                hasData = true;
            } catch (IOException e) {
                e.printStackTrace();
                return null;
            }
        }
        Map<String, String> values = data.get(parameter.toUpperCase());
        return values.get(chemicalAgent.toUpperCase());
    }

    public static double getDouble(String parameter, String chemicalAgent) {
        return Double.parseDouble(getData(parameter, chemicalAgent));
    }

    public static double getD(String chemicalAgent, double T) {
        double sigmaIN2 = (getDouble("sigma", chemicalAgent) + getDouble("sigma", "N2")) / 2;
        double epsIN2 = Math.sqrt(getDouble("eps", chemicalAgent) * getDouble("eps", "N2"));
        double omega11 = 1.074 * Math.pow(T / epsIN2, -0.1604);
        double muIN2 = 2 * getDouble("mu", chemicalAgent) * getDouble("mu", "N2") / (getDouble("mu", chemicalAgent) + getDouble("mu", "N2"));
        return 0.02628 * Math.pow(T, 1.5) / (PRESSURE * sigmaIN2 * omega11 * Math.sqrt(muIN2));
    }

    public static double getG(String chemicalAgent, double T) {
        double x = T / 10000;
        return getDouble("H", chemicalAgent) - T * (getDouble("f1", chemicalAgent) + getDouble("f2", chemicalAgent) * Math.log(x)
            + getDouble("f3", chemicalAgent) / (x * x) + getDouble("f4", chemicalAgent) / x + getDouble("f5", chemicalAgent) * x
            + getDouble("f6", chemicalAgent) * x * x + getDouble("f7", chemicalAgent) * x * x * x);
    }

    public static double getK(double T, int number) {
        switch (number) {
            case 1:
                return Math.exp(-(2 * getG("HCl", T) + 2 * getG("Al", T) - 2 * getG("AlCl", T) - getG("H2", T)) / (R * T)) / PRESSURE;
            case 2:
                return Math.exp(-(2 * getG("HCl", T) + getG("Al", T) - getG("AlCl2", T) - getG("H2", T)) / (R * T));
            case 3:
                return Math.exp(-(6 * getG("HCl", T) + 2 * getG("Al", T) - 2 * getG("AlCl3", T) - 3 * getG("H2", T)) / (R * T)) * PRESSURE;
            case 4:
                return Math.exp(-(2 * getG("HCl", T) + 2 * getG("Ga", T) - 2 * getG("GaCl", T) - getG("H2", T)) / (R * T)) / PRESSURE;
            case 5:
                return Math.exp(-(2 * getG("HCl", T) + getG("Ga", T) - getG("GaCl2", T) - getG("H2", T)) / (R * T));
            case 6:
                return Math.exp(-(6 * getG("HCl", T) + 2 * getG("Ga", T) - 2 * getG("GaCl3", T) - 3 * getG("H2", T)) / (R * T)) * PRESSURE;
            case 7:
                return Math.exp(-(getG("AlCl", T) + getG("NH3", T) - getG("AlN", T) - getG("HCl", T) - getG("H2", T)) / (R * T));
            case 8:
                return Math.exp(-(2 * getG("AlCl", T) + 2 * getG("NH3", T) - 2 * getG("AlN", T) - 4 * getG("HCl", T) - getG("H2", T)) / (R * T)) / PRESSURE;
            case 9:
                return Math.exp(-(getG("AlCl3", T) + getG("NH3", T) - getG("AlN", T) - 3 * getG("HCl", T)) / (R * T)) / PRESSURE;
            case 10:
                return Math.exp(-(getG("GaCl", T) + getG("NH3", T) - getG("GaN", T) - getG("HCl", T) - getG("H2", T)) / (R * T));
            case 11:
                return Math.exp(-(2 * getG("GaCl", T) + 2 * getG("NH3", T) - 2 * getG("GaN", T) - 4 * getG("HCl", T) - getG("H2", T)) / (R * T)) / PRESSURE;
            case 12:
                return Math.exp(-(getG("GaCl3", T) + getG("NH3", T) - getG("GaN", T) - 3 * getG("HCl", T)) / (R * T)) / PRESSURE;
            default:
                return 0;
        }
    }
}
