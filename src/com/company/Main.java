package com.company;

import java.util.Arrays;

public class Main {



    public static void main(String[] args) {
        Function[] functions = new Function[3];
        functions[0] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 9;
            }
        };
        functions[1] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] * x[1] - 4;
            }
        };
        functions[2] = new Function() {
            @Override
            public double calculate(double[] x) {
                return x[0] + x[1] * x[2] - 4;
            }
        };
        EquationSystem equationSystem = new EquationSystem(functions);
        System.out.println("Newton method:");
        double[] answer = equationSystem.newtonMethod(1e-12, 1000000);
        System.out.println(Arrays.toString(answer));
        System.out.println("Universal method");
        answer = equationSystem.universalMethod(1e-12, 1000000);
        System.out.println(Arrays.toString(answer));
    }
}
