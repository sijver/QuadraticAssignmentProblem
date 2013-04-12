package main.core;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class NeighbourhoodMatrix {

    private boolean[][] neighbourhoodMatrix;

    public NeighbourhoodMatrix(int solutionsNum) {
        neighbourhoodMatrix = new boolean[solutionsNum][solutionsNum];
    }

    public void print() {
        for (boolean[] row : neighbourhoodMatrix) {
            for (boolean b : row) {
                System.out.print(b + " ");
            }
            System.out.println();
        }
    }

    public void setCell(int i, int j, boolean value) {
        neighbourhoodMatrix[i][j] = value;
        neighbourhoodMatrix[j][i] = value;
    }

    public boolean getCell(int i, int j) {
        return neighbourhoodMatrix[i][j];
    }

    public List<Short> getAllNeighboursOfSolution(int solutionNumber) {
        List<Short> allNeighbours = new LinkedList<Short>();

        for (short i = 0; i < neighbourhoodMatrix.length; i++) {
            if (neighbourhoodMatrix[solutionNumber][i]) {
                allNeighbours.add(i);
            }
        }

        return allNeighbours;
    }


}
