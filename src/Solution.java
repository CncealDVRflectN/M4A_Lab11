public class Solution {
    private static double calcF(double x) {
        return 1.5 * Math.pow(Math.E, x) - 0.5 * Math.cos(x);
    }

    private static double calcFDeriv(double x) {
        return 1.5 * Math.pow(Math.E, x) + 0.5 * Math.sin(x);
    }

    private static void printArr(double[] arr) {
        for (int i = 0; i < arr.length; i++) {
            System.out.print(arr[i] + " ");
        }
        System.out.println();
    }

    private static double[] systemSolve(double[][] left, double[] right, int n) {
        double[] alpha = new double[n - 1];
        double[] beta = new double[n];
        double[] result = new double[n];
        alpha[0] = -left[0][1] / left[0][0];
        for (int i = 1; i < n - 1; i++) {
            alpha[i] = -left[i][i + 1] / (left[i][i] - (-left[i][i - 1]) * alpha[i - 1]);
        }
        beta[0] = right[0] / left[0][0];
        for (int i = 1; i < n; i++) {
            beta[i] = right[i] + (-left[i][i - 1]) * beta[i - 1];
            beta[i] /= left[i][i] - (-left[i][i - 1]) * alpha[i - 1];
        }
        result[n - 1] = beta[n - 1];
        for (int i = n - 2; i > -1; i--) {
            result[i] = alpha[i] * result[i + 1] + beta[i];
        }
        return result;
    }

    private static double calcSpline(double[][] table, int N, double h, double x, int index) {
        double[][] moments = new double[N + 1][N + 1];
        double[] right = new double[N + 1];
        double[] momentsCalculated;
        moments[0][0] = 2;
        moments[0][1] = 1;
        right[0] = 6 * ((table[1][1] - table[1][0]) / (table[0][1] - table[0][0]) - calcFDeriv(table[0][0])) / (table[0][1] - table[0][0]);
        moments[N][N - 1] = 1;
        moments[N][N] = 2;
        right[N] = 6 * (calcFDeriv(table[0][N]) - (table[1][N] - table[1][N - 1]) / (table[0][N] - table[0][N - 1])) / (table[0][N] - table[0][N - 1]);
        for (int i = 1; i < N; i++) {
            moments[i][i] = 2 * h / 3;
            moments[i][i - 1] = h / 6;
            moments[i][i + 1] = h / 6;
            right[i] = (table[1][i + 1] - table[1][i]) / h - (table[1][i] - table[1][i - 1]) / h;
        }
        momentsCalculated = systemSolve(moments, right, N + 1);
        return momentsCalculated[index - 1] * Math.pow(table[0][index] - x, 3) / (6 * h) + momentsCalculated[index] *
                Math.pow(x - table[0][index - 1], 3) / (6 * h) + (table[1][index - 1] - momentsCalculated[index - 1] * Math.pow(h, 2) / 6) *
                (table[0][index] - x) / h + (table[1][index] - momentsCalculated[index] * Math.pow(h, 2) / 6) * (x - table[0][index - 1]) / h;
    }

    public static void main(String[] args) {
        double[][] table = new double[2][11];
        int N = 10;
        for (int i = 0; i <= N; i++) {
            table[0][i] = (double)i / N;
            table[1][i] = calcF(table[0][i]);
        }
        System.out.print("X: ");
        printArr(table[0]);
        System.out.print("Y: ");
        printArr(table[1]);
        System.out.println("Кубический сплайн:");
        System.out.println("r*: " + Math.abs(calcSpline(table, N, 0.1, 0.1 / 3, 1) - calcF(0.1 / 3)));
        System.out.println("r**: " + Math.abs(calcSpline(table, N, 0.1, 0.5 + 0.1 / 3, 5) - calcF(0.5 + 0.1 / 3)));
        System.out.println("r***: " + Math.abs(calcSpline(table, N, 0.1, 1 - 0.1 / 3, 10) - calcF(1 - 0.1 / 3)));
    }
}
