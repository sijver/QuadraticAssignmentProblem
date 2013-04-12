package main.core;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class MacroState {

    private List<Short> solutions;

    private int solutionsEvaluation;

    public MacroState(int solutionsEvaluation) {
        this.solutionsEvaluation = solutionsEvaluation;
        solutions = new LinkedList<Short>();
    }

    public void addSolutionToMacroState(Short solutionNumber){
        solutions.add(solutionNumber);
    }

    public List<Short> getSolutions() {
        return solutions;
    }

    public int getSolutionsEvaluation() {
        return solutionsEvaluation;
    }
}
