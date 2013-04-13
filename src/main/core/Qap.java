package main.core;

import main.core.utils.QapInstanceReader;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 */
public class Qap {

    int[][] dMatrix; // First matrix read from input file, typically distance matrix
    int[][] fMatrix; // Second matrix read from input file, typically flow matrix

    byte[] solution; // Array containing the current solution. If dMatrix is the distance matrix and fMatrix the flow matrix, solution[i] is the item assigned to location i

    private byte instanceSize;   // Instance size of the QAP instance

    boolean firstImproveFlag;   // First_improv version of local search

    long bestKnown;  // Best objective function value of QAP instance
    long bestFound;  // Best solution found so far with local search
    long veryBestFound;  // Best solution found so far with local search

    boolean dSymmetricFlag = false; // If first (d) matrix is symmetric: TRUE
    boolean fSymmetricFlag = false; // If second (f) matrix is symmetric: TRUE
    boolean nullDiagonalFlag = false;   // At least one matrix has zero diagonal: TRUE
    boolean makeSymmetricFlag = false;  // Convert asymmetric instance into symmetric instance: TRUE


    private List<byte[]> allSolutions;
    private List<SolutionType> solutionTypes;
    private EnumMap<SolutionType, Double> solutionTypeStatistics;
    private List<MacroState> collapsedSearchLandscape;
    private List<MacroStatePre> collapsedSearchLandscapePre;

    public Qap(byte instanceSize, int[][] dMatrix, int[][] fMatrix) {
        this.instanceSize = instanceSize;
        this.dMatrix = dMatrix;
        this.fMatrix = fMatrix;
        initProgram();
    }

    public static void main(String[] args) {

//        start_timers();             /* start timing routines */
//        init_program(argc, argv);   /* initialize all important data */

        Qap qap = QapInstanceReader.readQapInstance("aa.dat");

        qap.veryBestFound = Long.MAX_VALUE;

        for (int i = 0; i < 100; i++) {
            qap.solution = qap.generateRandomVector(qap.instanceSize);

            qap.bestFound = qap.computeEvaluationFunction(qap.solution);
            if (qap.dSymmetricFlag && qap.fSymmetricFlag && qap.nullDiagonalFlag) {
                qap.best2OptSymmetric(qap.solution);
            } else if (qap.makeSymmetricFlag) {
                qap.best2OptSymmetric(qap.solution);
            } else {
                qap.best2OptAsymmetric(qap.solution);
            }

            if (qap.computeEvaluationFunction(qap.solution) != qap.bestFound) {
                System.out.println("Some error must have occurred in local search routine,\n values do not match");
            }

            if (qap.bestFound < qap.veryBestFound) {
                qap.veryBestFound = qap.bestFound;
            }
        }

        System.out.println(String.format("best solution best improvement: %1$d", qap.veryBestFound));

        for (int i = 0; i < 100; i++) {
            qap.solution = qap.generateRandomVector(qap.instanceSize);

            qap.bestFound = qap.computeEvaluationFunction(qap.solution);
            if (qap.dSymmetricFlag && qap.fSymmetricFlag && qap.nullDiagonalFlag) {
                qap.first2OptSymmetric(qap.solution);
            } else if (qap.makeSymmetricFlag) {
                qap.first2OptSymmetric(qap.solution);
            } else {
                qap.first2OptAssymetric(qap.solution);
            }

            if (qap.computeEvaluationFunction(qap.solution) != qap.bestFound) {
                System.out.println("Some error must have occurred in local search routine,\n values do not match");
            }

            if (qap.bestFound < qap.veryBestFound) {
                qap.veryBestFound = qap.bestFound;
            }
        }

        System.out.println(String.format("best solution first improvement: %1$d", qap.veryBestFound));


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
    public void printSolution(byte[] solution) {
        System.out.println("Assignment: ");
        for (byte i : solution) {
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
    public long computeEvaluationFunction(byte[] solution) {
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
//        System.out.println(String.format("Objective function value = %d", objectiveFunctionValue));
//        System.out.println();
        return objectiveFunctionValue;
    }

    // Generates a random vector
    public byte[] generateRandomVector(int size) {
        byte[] vector = new byte[size];
        byte help, j;
        Random random = new Random();

        for (byte i = 0; i < size; i++) {
            vector[i] = i;
        }

        for (int i = 0; i < size - 1; i++) {
            j = (byte) random.nextInt(size - i);
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

    // Swap items at positions pos1 and pos2
    public void swap(int pos1, int pos2, byte[] vector) {
        byte help;
        help = vector[pos1];
        vector[pos1] = vector[pos2];
        vector[pos2] = help;
    }

    // First improvement 2-opt local search for asymmetric instances
    public void first2OptAssymetric(byte[] vector) {
        boolean improvement = true;
        int u, v;
        int tmp;

        System.out.println("First imp, asymmetric case");

        bestFound = computeEvaluationFunction(vector);
        byte[] xVector = generateRandomVector(instanceSize); // Random vector to scan neighborhood in random order
        while (improvement) {
            improvement = true;
            for (int i = 0; i < instanceSize; i++) {
                u = xVector[i];
                for (int j = 0; j < instanceSize; j++) {
                    v = xVector[j];
                    if (u == v) {
                        continue;
                    }
                    tmp = 0;
                    for (int k = 0; k < instanceSize; k++) {
                        if ((k != u) && (k != v)) {
                            tmp += dMatrix[k][u] * (fMatrix[vector[k]][vector[v]] - fMatrix[vector[k]][vector[u]]) +
                                    dMatrix[k][v] * (fMatrix[vector[k]][vector[u]] - fMatrix[vector[k]][vector[v]]) +
                                    dMatrix[u][k] * (fMatrix[vector[v]][vector[k]] - fMatrix[vector[u]][vector[k]]) +
                                    dMatrix[v][k] * (fMatrix[vector[u]][vector[k]] - fMatrix[vector[v]][vector[k]]);
                        }
                    }
                    tmp += dMatrix[u][u] * (fMatrix[vector[v]][vector[v]] - fMatrix[vector[u]][vector[u]]) +
                            dMatrix[u][v] * (fMatrix[vector[v]][vector[u]] - fMatrix[vector[u]][vector[v]]) +
                            dMatrix[v][u] * (fMatrix[vector[u]][vector[v]] - fMatrix[vector[v]][vector[u]]) +
                            dMatrix[v][v] * (fMatrix[vector[u]][vector[u]] - fMatrix[vector[v]][vector[v]]);
                    if (tmp < 0) {
                        improvement = true;
                        bestFound += tmp;
                        swap(u, v, vector);
                        System.out.println(String.format("Improvement %d, bestKnown %d", tmp, bestFound));
                    }
                }
            }
        }
    }

    // First improvement 2-opt local search for symmetric instances
    public void first2OptSymmetric(byte[] vector) {
        boolean improvement = true;
        int u, v;
        int tmp;
        byte[] xVector = new byte[instanceSize];  // Scan neighborhood in random order
        int originalSymmetricFactor;  //2: original symmetric instance, 1: original asymmetric instance

        System.out.println("First imp, symmetric case");
        if (makeSymmetricFlag) {
            originalSymmetricFactor = 1;  // Compensation because of not dividing matrix by 2
        } else {
            originalSymmetricFactor = 2;
        }
        bestFound = computeEvaluationFunction(vector);
        improvement = true;
        xVector = generateRandomVector(instanceSize);
        while (improvement) {
            improvement = true;
            for (int i = 0; i < instanceSize; i++) {
                u = xVector[i];
                for (int j = 0; j < instanceSize; j++) {
                    v = xVector[j];
                    if (u == v) {
                        continue;
                    }
                    tmp = 0;
                    for (int k = 0; k < instanceSize; k++) {
                        if ((k != u) && (k != v)) {
                            tmp += (dMatrix[k][u] - dMatrix[k][v]) * (fMatrix[vector[k]][vector[v]] - fMatrix[vector[k]][vector[u]]);
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
    public void best2OptAsymmetric(byte[] vector) {
        boolean improvement = true;
        int tmp;
        int originalSymmetricFactor;  //2: original symmetric instance, 1: original asymmetric instance

//        int[][] moveValues = new int[instanceSize][instanceSize]; // Matrix of move values in previous iteration allows for fast evaluation of neighbourhood

        boolean firstItFlag = true; // First iteration of local search: true
        long maxDecrease;    // Largest decrease found so far in neighbourhood scan
        int rchosen = instanceSize, schosen = instanceSize; // Memorize which is best move in current iteration
        int r = -1, s = -1;   // Memorize which is best move in previous iteration

        System.out.println("Best imp, asymmetric case");

        if (makeSymmetricFlag) {
            originalSymmetricFactor = 1;
        } else {
            originalSymmetricFactor = 2;
        }
        bestFound = computeEvaluationFunction(vector);

//        for (int k = 0; k < instanceSize; k++ ) {
//            moveValues[k] = (int[])(moveValues + instanceSize) + k * instanceSize;
//        }
        while (improvement) {
            improvement = false;
            maxDecrease = Long.MAX_VALUE;
            // In the first local search iteration the full neighborhood has to be evaluated
            if (firstItFlag) {
                firstItFlag = false;
                for (int u = 0; u < instanceSize - 1; u++) {
                    for (int v = u + 1; v < instanceSize; v++) {
                        tmp = 0;
                        for (int k = 0; k < instanceSize; k++) {
                            if ((k != u) && (k != v)) {
                                tmp += dMatrix[k][u] * (fMatrix[vector[k]][vector[v]] - fMatrix[vector[k]][vector[u]]) +
                                        dMatrix[k][v] * (fMatrix[vector[k]][vector[u]] - fMatrix[vector[k]][vector[v]]) +
                                        dMatrix[u][k] * (fMatrix[vector[v]][vector[k]] - fMatrix[vector[u]][vector[k]]) +
                                        dMatrix[v][k] * (fMatrix[vector[u]][vector[k]] - fMatrix[vector[v]][vector[k]]);
                            }
                        }
                        tmp += dMatrix[u][u] * (fMatrix[vector[v]][vector[v]] - fMatrix[vector[u]][vector[u]]) +
                                dMatrix[u][v] * (fMatrix[vector[v]][vector[u]] - fMatrix[vector[u]][vector[v]]) +
                                dMatrix[v][u] * (fMatrix[vector[u]][vector[v]] - fMatrix[vector[v]][vector[u]]) +
                                dMatrix[v][v] * (fMatrix[vector[u]][vector[u]] - fMatrix[vector[v]][vector[v]]);
//                        moveValues[u][v] = tmp;
                        if (tmp < maxDecrease) {
                            maxDecrease = tmp;
                            rchosen = u;
                            schosen = v;
                        }
                    }
                }
            } else {
                for (int u = 0; u < instanceSize - 1; u++) {
                    for (int v = u + 1; v < instanceSize; v++) {
                        if (u == r || v == s || u == s || v == r) {
                            tmp = 0;
                            for (int k = 0; k < instanceSize; k++) {
                                if ((k != u) && (k != v)) {
                                    tmp += dMatrix[k][u] * (fMatrix[vector[k]][vector[v]] - fMatrix[vector[k]][vector[u]]) +
                                            dMatrix[k][v] * (fMatrix[vector[k]][vector[u]] - fMatrix[vector[k]][vector[v]]) +
                                            dMatrix[u][k] * (fMatrix[vector[v]][vector[k]] - fMatrix[vector[u]][vector[k]]) +
                                            dMatrix[v][k] * (fMatrix[vector[u]][vector[k]] - fMatrix[vector[v]][vector[k]]);
                                }
                            }
                            tmp += dMatrix[u][u] * (fMatrix[vector[v]][vector[v]] - fMatrix[vector[u]][vector[u]]) +
                                    dMatrix[u][v] * (fMatrix[vector[v]][vector[u]] - fMatrix[vector[u]][vector[v]]) +
                                    dMatrix[v][u] * (fMatrix[vector[u]][vector[v]] - fMatrix[vector[v]][vector[u]]) +
                                    dMatrix[v][v] * (fMatrix[vector[u]][vector[u]] - fMatrix[vector[v]][vector[v]]);
//                            moveValues[u][v] = tmp;
                            if (tmp < maxDecrease) {
                                maxDecrease = tmp;
                                rchosen = u;
                                schosen = v;
                            }
                        } else { /* Change derived from move_values */
                            tmp = (dMatrix[r][u] - dMatrix[r][v] + dMatrix[s][v] - dMatrix[s][u]) *
                                    (fMatrix[vector[s]][vector[u]] - fMatrix[vector[s]][vector[v]] + fMatrix[vector[r]][vector[v]] - fMatrix[vector[r]][vector[u]])
                                    + (dMatrix[u][r] - dMatrix[v][r] + dMatrix[v][s] - dMatrix[u][s]) *
                                    (fMatrix[vector[u]][vector[s]] - fMatrix[vector[v]][vector[s]] + fMatrix[vector[v]][vector[r]] - fMatrix[vector[u]][vector[r]]);
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
            if (maxDecrease < 0) {      /* Objective function value can be improved */
                assert (rchosen < schosen);
                improvement = true;
                bestFound += maxDecrease;
                swap(rchosen, schosen, vector);
                r = rchosen;    // Memorize previously done move
                s = schosen;    // Memorize previously done move
                System.out.println(String.format("Improvement %d, bestFound %d, exchange %d and %d", maxDecrease, bestFound, rchosen, schosen));
            }
        }
    }

    public void best2OptSymmetric(byte[] vector) {

    }

    public void initProgram() {
        dSymmetricFlag = checkSymmetry(dMatrix, instanceSize);
        fSymmetricFlag = checkSymmetry(fMatrix, instanceSize);
        nullDiagonalFlag = checkNullDiagonal(dMatrix, instanceSize);
        if (!nullDiagonalFlag) {
            nullDiagonalFlag = checkNullDiagonal(fMatrix, instanceSize);
        }

        makeSymmetricFlag = (dSymmetricFlag != fSymmetricFlag);
        if (makeSymmetricFlag && nullDiagonalFlag) {
            if (!dSymmetricFlag)
                makeMatrixSymmetric(dMatrix, instanceSize);
            else if (!fSymmetricFlag)
                makeMatrixSymmetric(fMatrix, instanceSize);
            else {
                System.out.println("One matrix should have been symmetric");
                System.exit(1);
            }
        }
    }

    public double ran01() {
        //TODO
        return 0;
    }


    public boolean areNeighbours(byte[] vector1, byte[] vector2) {
        int changes = 0;

        if (vector1.length != vector2.length) {
            return false;
        }

        for (int i = 0; i < vector1.length; i++) {
            if (vector1[i] != vector2[i]) {
                changes++;
                if (changes > 2) {
                    return false;
                }
            }
        }

        return changes == 2;
    }

    public void computeAllSolutions() {
        if (allSolutions == null) {
            allSolutions = getVectorPermutations(null);
        }
    }

    public List<byte[]> getAllSolutions() {
        return allSolutions;
    }

    private List<byte[]> getVectorPermutations(List<byte[]> previousList) {
        if (previousList == null) {
            previousList = new ArrayList<byte[]>();
        }
        List<byte[]> newList = new ArrayList<byte[]>();
        byte previousLength = 0;

        if (previousList.size() == 0) {
            newList.add(new byte[]{0});
        } else {
            previousLength = (byte) previousList.get(0).length;
            System.out.println(previousLength);
            for (byte[] vector : previousList) {
                for (int j = 0; j < previousLength + 1; j++) {
                    byte[] newVector = Arrays.copyOf(vector, previousLength + 1);

                    for (int k = previousLength; k > j; k--) {
                        newVector[k] = newVector[k - 1];
                    }

                    newVector[j] = previousLength;

                    newList.add(newVector);
                }
            }
        }

        if (previousLength + 1 == instanceSize) {
            return newList;
        } else {
            return getVectorPermutations(newList);
        }
    }

    public static List<byte[]> getAllNeighboursOfSolution(byte[] vector) {
        List<byte[]> neighboursList = new LinkedList<byte[]>();

        for (int i = 0; i < vector.length; i++) {
            for (int j = i + 1; j < vector.length; j++) {
                byte[] newVector = Arrays.copyOf(vector, vector.length);
                byte help = newVector[i];
                newVector[i] = newVector[j];
                newVector[j] = help;
                neighboursList.add(newVector);
            }
        }

        return neighboursList;
    }

    public void computeAllSolutionTypes() {
        if (solutionTypes == null) {
            computeAllSolutions();
            solutionTypes = new LinkedList<SolutionType>();

            collapsedSearchLandscape = new ArrayList<MacroState>();
            collapsedSearchLandscapePre = new LinkedList<MacroStatePre>();

            for (int i = 0; i < allSolutions.size(); i++) {
                int down = 0;
                int up = 0;
                int side = 0;

                byte[] currentSolution = allSolutions.get(i);
                long currentSolutionEvaluation = computeEvaluationFunction(currentSolution);

                List<byte[]> currentSolutionNeighbours = getAllNeighboursOfSolution(currentSolution);

                for (byte[] currentSolutionNeighbour : currentSolutionNeighbours) {
                    long newSolutionEvaluation = computeEvaluationFunction(currentSolutionNeighbour);
                    if (newSolutionEvaluation < currentSolutionEvaluation) {
                        down++;
                    } else if (newSolutionEvaluation == currentSolutionEvaluation) {
                        side++;
                        collapsedSearchLandscapePre.add(new MacroStatePre(currentSolutionEvaluation, i, currentSolutionNeighbour));
                    } else if (newSolutionEvaluation > currentSolutionEvaluation) {
                        up++;
                    }
                }

                if (down == 0 && side == 0) {
                    solutionTypes.add(SolutionType.SLMIN);
                } else if (down == 0 && side > 0 && up > 0) {
                    solutionTypes.add(SolutionType.LMIN);
                } else if (down == 0 && up == 0) {
                    solutionTypes.add(SolutionType.IPLAT);
                    System.out.println(SolutionType.IPLAT);
                } else if (down > 0 && side > 0 && up > 0) {
                    solutionTypes.add(SolutionType.LEDGE);
                } else if (down > 0 && side == 0 && up > 0) {
                    solutionTypes.add(SolutionType.SLOPE);
                } else if (down > 0 && side > 0 && up == 0) {
                    solutionTypes.add(SolutionType.LMAX);
                } else if (side == 0 && up == 0) {
                    solutionTypes.add(SolutionType.SLMAX);
                }
            }
        }

        for (MacroStatePre macroStatePre : collapsedSearchLandscapePre) {
            boolean isNeighbourFind = false;
            for (MacroState macroState : collapsedSearchLandscape) {
                if (macroState.getSolutionsEvaluation() == macroStatePre.getSolutionsEvaluation()) {
                    for (int i = 0; i < macroState.getSolutions().size(); i++) {
                        if (areNeighbours(allSolutions.get(macroState.getSolutions().get(i)), allSolutions.get(macroStatePre.getSolutionNum())) && !macroState.getSolutions().contains(macroStatePre.getSolutionNum())) {
                            isNeighbourFind = true;
                            macroState.addSolutionToMacroState(macroStatePre.getSolutionNum());
                            break;
                        }
                    }
                }
                if(macroState.getSolutions().contains(macroStatePre.getSolutionNum())){
                    isNeighbourFind = true;
                    break;
                }
            }
            if(!isNeighbourFind){
                MacroState macroState = new MacroState(macroStatePre.getSolutionsEvaluation());
                macroState.addSolutionToMacroState(macroStatePre.getSolutionNum());
                collapsedSearchLandscape.add(macroState);
            }
        }

        for(int i = 0; i < collapsedSearchLandscape.size(); i++){
            for(int j = i + 1; j < collapsedSearchLandscape.size(); j++){
                if(collapsedSearchLandscape.get(i).getSolutionsEvaluation() == collapsedSearchLandscape.get(j).getSolutionsEvaluation()){
                    boolean areNeighbours = false;

                    for(int sol1Num : collapsedSearchLandscape.get(i).getSolutions()){
                        for(int sol2Num : collapsedSearchLandscape.get(j).getSolutions()){
                            if(areNeighbours(allSolutions.get(sol1Num), allSolutions.get(sol2Num))){
                                areNeighbours = true;
                            }
                        }
                    }

                    if(areNeighbours){
                        for(int solNum : collapsedSearchLandscape.get(j).getSolutions()){
                            collapsedSearchLandscape.get(i).addSolutionToMacroState(solNum);
                        }
                        collapsedSearchLandscape.remove(j);
                        j--;
                    }
                }
            }
        }
    }

    public void computeSolutionTypeStatistics() {
        if (solutionTypeStatistics == null) {
            computeAllSolutionTypes();
            solutionTypeStatistics = new EnumMap<SolutionType, Double>(SolutionType.class);
            for (SolutionType solutionType : SolutionType.values()) {
                solutionTypeStatistics.put(solutionType, 0.0);
            }
            for (SolutionType solutionType : solutionTypes) {
                solutionTypeStatistics.put(solutionType, solutionTypeStatistics.get(solutionType) + 1);
            }
            for (SolutionType solutionType : SolutionType.values()) {
                solutionTypeStatistics.put(solutionType, solutionTypeStatistics.get(solutionType) / solutionTypes.size() * 100);
            }
        }
    }

    public EnumMap<SolutionType, Double> getSolutionTypeStatistics() {
        computeSolutionTypeStatistics();
        return solutionTypeStatistics;
    }

    public List<MacroState> getCollapsedSearchLandscape() {
        return collapsedSearchLandscape;
    }

    public List<MacroStatePre> getCollapsedSearchLandscapePre() {
        return collapsedSearchLandscapePre;
    }
}
