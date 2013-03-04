import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 */
public class Qap {

    int[][] dMatrix; // First matrix read from input file, typically distance matrix
    int[][] fMatrix; // Second matrix read from input file, typically flow matrix

    int[] solution; // Array containing the current solution. If dMatrix is the distance matrix and fMatrix the flow matrix, solution[i] is the item assigned to location i

    private int instanceSize;   // Instance size of the QAP instance

    long bestKnown;  // Best objective function value of QAP instance
    long bestFound;  // Best solution found so far with local search
    long veryBestFound;  // Best solution found so far with local search

    private boolean makeSymmetricFlag = false;  // Convert asymmetric instance into symmetric

    public static void main(String[] args) {

//        start_timers();             /* start timing routines */
//        init_program(argc, argv);   /* initialize all important data */

    }

    // Check whether the matrix is symmetric
    private boolean checkSymmetry(int[][] matrix, int size) {
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                if (matrix[i][j] != matrix[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }

    // Check whether the matrix has a zero diagonal
    private boolean checkNullDiagonal(int[][] matrix, int size) {
        for (int i = 0; i < size; i++) {
            if (matrix[i][i] != 0) {
                return false;
            }
        }
        return true;
    }

    // Makes an asymmetric matrix symmetric (calculates M = M + M-transpose)
    private int[][] makeMatrixSymmetric(int[][] matrix, int size) {
        int help;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < i; j++) {
                help = matrix[i][j] + matrix[j][i];
                matrix[i][j] = help;
                matrix[j][i] = help;
            }
        }
        return matrix;
    }

    // Prints solution
    public void printSolution(int[] solution) {
        System.out.println("Assignment:");
        for (int i : solution) {
            System.out.print(String.format("%d ", i));
        }
        System.out.println();
    }

    // Prints matrix
    public void printMatrix(int[][] matrix) {
        System.out.println();
        for (int i = 0; i < instanceSize; i++) {
            for (int j = 0; j < instanceSize; j++) {
                System.out.print(String.format("%d ", matrix[i][j]));
            }
            System.out.println();
        }
    }

    // Computes the objective function value for the QAP
    public long computeEvaluationFunction(int[] solution) {
        long objectiveFunctionValue = 0;
        for (int i = 0; i < instanceSize; i++) {
            for (int j = 0; j < instanceSize; j++) {
                objectiveFunctionValue += dMatrix[i][j] * fMatrix[solution[i]][solution[j]];
            }
        }
        if (makeSymmetricFlag) {
            /* Division by 2 has to be done if we have a n asymmetric instance which
            has been converted into a symmetric one (indicated by makeSymmetricFlag).
                    This is due to the particular way of doing this conversion. */
            objectiveFunctionValue /= 2;
        }
        System.out.println(String.format("Objective function value = %d", objectiveFunctionValue));
        System.out.println();
        return objectiveFunctionValue;
    }

    // Generates a random vector
    public int[] generateRandomVector(int size) {
        int[] vector = new int[size];
        int help, j;
        Random random = new Random();

        for(int i = 0; i < size; i++){
            vector[i] = i;
        }

        for (int i = 0; i < size - 1; i++) {
            j = random.nextInt(size - i);
            help = vector[i];
            vector[i] = vector[i + j];
            vector[i + j] = help;
        }

        System.out.println("Random vector:");
        for (int i = 0; i < size; i++) {
            System.out.print(String.format("%d ", vector[i]));
        }
        System.out.println();
        return vector;
    }

    // Swap items at positions pos11 and pos2
    public void swap(int pos1, int pos2, int[] vector){
        int help;
        help = vector[pos1];
        vector[pos1] = vector[pos2];
        vector[pos2] = help;
    }

    // First improvement 2-opt local search for asymmetric instances
    public void first2OptAssymetric(int[] vector){
        boolean improvement = true;
        int u, v;
        int tmp;
        int[] xVector = new int[instanceSize];  // Random vector to scan neighborhood in random order

        System.out.println("First imp, asymmetric case");

        bestFound  = computeEvaluationFunction(vector);
        xVector = generateRandomVector(instanceSize);
        while (improvement) {
            improvement = true;
            for (int i = 0 ; i < instanceSize; i++) {
                u = xVector[i];
                for (int j = 0 ; j < instanceSize; j++) {
                    v = xVector[j];
                    if (u == v){
                        continue;
                    }
                    tmp = 0;
                    for (int k = 0; k < instanceSize; k++ ) {
                        if ( (k != u) && (k != v) ) {
                            tmp += dMatrix[k][u] * ( fMatrix[vector[k]][vector[v]] - fMatrix[vector[k]][vector[u]] ) +
                                    dMatrix[k][v] * ( fMatrix[vector[k]][vector[u]] - fMatrix[vector[k]][vector[v]] ) +
                                    dMatrix[u][k] * ( fMatrix[vector[v]][vector[k]] - fMatrix[vector[u]][vector[k]] ) +
                                    dMatrix[v][k] * ( fMatrix[vector[u]][vector[k]] - fMatrix[vector[v]][vector[k]] );
                        }
                    }
                    tmp += dMatrix[u][u] * ( fMatrix[vector[v]][vector[v]] - fMatrix[vector[u]][vector[u]] )+
                            dMatrix[u][v] * ( fMatrix[vector[v]][vector[u]] - fMatrix[vector[u]][vector[v]] )+
                            dMatrix[v][u] * ( fMatrix[vector[u]][vector[v]] - fMatrix[vector[v]][vector[u]] )+
                            dMatrix[v][v] * ( fMatrix[vector[u]][vector[u]] - fMatrix[vector[v]][vector[v]] );
                    if (tmp < 0) {
                        improvement = true;
                        bestFound += tmp;
                        swap(u, v, vector);
                        System.out.println(String.format("Improvement %d, bestKnown %d", tmp, bestFound));
                    }
                }
            }
        }
        free ( xVector );
    }

    // First improvement 2-opt local search for symmetric instances
    public void first2OptSymmetric(int[] vector){
        boolean improvement = true;
        int u, v;
        int tmp;
        int[] xVector = new int[instanceSize];  // Scan neighborhood in random order
        int originalSymmetricFactor;  //2: original symmetric instance, 1: original asymmetric instance

        System.out.println("First imp, symmetric case");
        if (makeSymmetricFlag){
            originalSymmetricFactor = 1;  // Compensation because of not dividing matrix by 2
        } else {
            originalSymmetricFactor = 2;
        }
        bestFound  = computeEvaluationFunction(vector);
        improvement = true;
        xVector = generateRandomVector(instanceSize);
        while (improvement) {
            improvement = true;
            for (int i = 0; i < instanceSize; i++) {
                u = xVector[i];
                for (int j = 0; j < instanceSize; j++) {
                    v = xVector[j];
                    if (u == v){
                        continue;
                    }
                    tmp = 0;
                    for (int k = 0; k < instanceSize; k++) {
                        if ( (k != u) && (k != v) ) {
                            tmp += ( dMatrix[k][u] - dMatrix[k][v] ) * ( fMatrix[vector[k]][vector[v]] - fMatrix[vector[k]][vector[u]] );
                        }
                    }
                    tmp *= originalSymmetricFactor;
                    if (tmp < 0) {
                        improvement = true;
                        bestFound += tmp;
                        swap(u, v, vector);
                        System.out.println(String.format("Improvement %d, bestFound %d", tmp, bestFound));
                    }
                }
            }
        }
    }

    // Best improvement 2-opt local search for asymmetric instances
    public void best2OptAsymmetric(int[] vector){
        boolean improvement = true;
        int tmp;
        int originalSymmetricFactor;  //2: original symmetric instance, 1: original asymmetric instance

//        int[][] moveValues = new int[instanceSize][instanceSize]; // Matrix of move values in previous iteration allows for fast evaluation of neighbourhood

        boolean firstItFlag = true; // First iteration of local search: true
        long maxDecrease;    // Largest decrease found so far in neighbourhood scan
        int rchosen = instanceSize, schosen = instanceSize; // Memorize which is best move in current iteration
        int r = -1, s = -1;   // Memorize which is best move in previous iteration

        System.out.println("Best imp, asymmetric case");

        if ( makeSymmetricFlag ){
            originalSymmetricFactor = 1;
        }else{
            originalSymmetricFactor = 2;
        }
        bestFound  = computeEvaluationFunction(vector);

//        for (int k = 0; k < instanceSize; k++ ) {
//            moveValues[k] = (int[])(moveValues + instanceSize) + k * instanceSize;
//        }
        while ( improvement ) {
            improvement = false;
            maxDecrease = Long.MAX_VALUE;
        // In the first local search iteration the full neighborhood has to be evaluated
            if (firstItFlag) {
                firstItFlag = false;
                for (int u = 0 ; u < instanceSize - 1 ; u++) {
                    for (int v = u + 1 ; v < instanceSize; v++) {
                        tmp = 0;
                        for (int k = 0 ; k < instanceSize; k++ ) {
                            if ( (k != u) && (k != v) ) {
                                tmp += dMatrix[k][u] * ( fMatrix[vector[k]][vector[v]] - fMatrix[vector[k]][vector[u]] ) +
                                        dMatrix[k][v] * ( fMatrix[vector[k]][vector[u]] - fMatrix[vector[k]][vector[v]] ) +
                                        dMatrix[u][k] * ( fMatrix[vector[v]][vector[k]] - fMatrix[vector[u]][vector[k]] ) +
                                        dMatrix[v][k] * ( fMatrix[vector[u]][vector[k]] - fMatrix[vector[v]][vector[k]] );
                            }
                        }
                        tmp += dMatrix[u][u] * ( fMatrix[vector[v]][vector[v]] - fMatrix[vector[u]][vector[u]] )+
                                dMatrix[u][v] * ( fMatrix[vector[v]][vector[u]] - fMatrix[vector[u]][vector[v]] )+
                                dMatrix[v][u] * ( fMatrix[vector[u]][vector[v]] - fMatrix[vector[v]][vector[u]] )+
                                dMatrix[v][v] * ( fMatrix[vector[u]][vector[u]] - fMatrix[vector[v]][vector[v]] );
//                        moveValues[u][v] = tmp;
                        if (tmp < maxDecrease) {
                            maxDecrease = tmp;
                            rchosen = u;
                            schosen = v;
                        }
                    }
                }
            } else {
                for (int u = 0 ; u < instanceSize - 1 ; u++) {
                    for (int v = u+1 ; v < instanceSize ; v++) {
                        if (u == r || v == s || u == s || v == r) {
                            tmp = 0;
                            for (int k = 0 ; k < instanceSize ; k++ ) {
                                if ( (k != u) && (k != v) ) {
                                    tmp += dMatrix[k][u] * ( fMatrix[vector[k]][vector[v]] - fMatrix[vector[k]][vector[u]] ) +
                                            dMatrix[k][v] * ( fMatrix[vector[k]][vector[u]] - fMatrix[vector[k]][vector[v]] ) +
                                            dMatrix[u][k] * ( fMatrix[vector[v]][vector[k]] - fMatrix[vector[u]][vector[k]] ) +
                                            dMatrix[v][k] * ( fMatrix[vector[u]][vector[k]] - fMatrix[vector[v]][vector[k]] );
                                }
                            }
                            tmp += dMatrix[u][u] * ( fMatrix[vector[v]][vector[v]] - fMatrix[vector[u]][vector[u]] )+
                                    dMatrix[u][v] * ( fMatrix[vector[v]][vector[u]] - fMatrix[vector[u]][vector[v]] )+
                                    dMatrix[v][u] * ( fMatrix[vector[u]][vector[v]] - fMatrix[vector[v]][vector[u]] )+
                                    dMatrix[v][v] * ( fMatrix[vector[u]][vector[u]] - fMatrix[vector[v]][vector[v]] );
//                            moveValues[u][v] = tmp;
                            if (tmp < maxDecrease) {
                                maxDecrease = tmp;
                                rchosen = u;
                                schosen = v;
                            }
                        } else { /* Change derived from move_values */
                            tmp = ( dMatrix[r][u] - dMatrix[r][v] + dMatrix[s][v] - dMatrix[s][u] ) *
                                    ( fMatrix[vector[s]][vector[u]] - fMatrix[vector[s]][vector[v]] + fMatrix[vector[r]][vector[v]] - fMatrix[vector[r]][vector[u]] )
                                    + ( dMatrix[u][r] - dMatrix[v][r] + dMatrix[v][s] - dMatrix[u][s] ) *
                                    ( fMatrix[vector[u]][vector[s]] - fMatrix[vector[v]][vector[s]] + fMatrix[vector[v]][vector[r]] - fMatrix[vector[u]][vector[r]] );
//                            tmp += moveValues[u][v];
//                            moveValues[u][v] = tmp;
                        }
                        if (tmp < maxDecrease) {
                            maxDecrease = tmp;
                            rchosen = u;
                            schosen = v;
                        }
                    }
                }
            }
            if ( maxDecrease < 0 ) {      /* Objective function value can be improved */
                assert (rchosen < schosen);
                improvement = true;
                bestFound += maxDecrease;
                swap(rchosen,schosen,vector);
                r = rchosen;    // Memorize previously done move
                s = schosen;    // Memorize previously done move
                System.out.println(String.format("Improvement %d, bestFound %d, exchange %d and %d", maxDecrease, bestFound, rchosen, schosen));
            }
        }
    }

    public void best2OptSymmetric(int[] vector){

    }

    public void initProgram() {

    }

    private void readInstance() {
        /*
        try {
            FileInputStream fileInputStream = new FileInputStream(fileName);
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(fileInputStream));
            String strLine;

            Pattern intPattern = Pattern.compile("\\d+");
            Matcher matcher;

            // Read file line by line
            while ((strLine = bufferedReader.readLine()) != null) {
                // Skip comments
                if (!strLine.startsWith("c")) {
                    matcher = intPattern.matcher(strLine);
                    // Read graph information
                    if (strLine.startsWith("p")) {
                        // Number of nodes
                        if (matcher.find()) {
                            instanceSize = Integer.parseInt(matcher.group());
                            adj_mat = new boolean[instanceSize][instanceSize];
                        }
                        // Number of edges
                        if (matcher.find()) {
                            edgesNum = Integer.parseInt(matcher.group());
                        }
                        // Reading edges
                    } else if (strLine.startsWith("e")) {
                        int n1, n2;
                        if (matcher.find()) {
                            n1 = Integer.parseInt(matcher.group()) - 1;
                            if (matcher.find()) {
                                n2 = Integer.parseInt(matcher.group()) - 1;
                                adj_mat[n1][n2] = true;
                                adj_mat[n2][n1] = true;
                            }
                        }
                    }
                }
            }
        } catch (FileNotFoundException e) {
            System.out.println("ERROR: unable to open file");
            System.exit(1);
        } catch (IOException e) {
            System.out.println("ERROR: I/O exception");
            System.exit(1);
        }
        */
    }


}
