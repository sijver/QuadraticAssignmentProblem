package main.core;

/**
 * Created with IntelliJ IDEA.
 */
public class MacroStatePre {

    private long solutionsEvaluation;

    private int solutionNum;

    private byte[] appropriateNeighbour;

    public MacroStatePre(long solutionsEvaluation, int solutionNum, byte[] appropriateNeighbour) {
        this.solutionsEvaluation = solutionsEvaluation;
        this.solutionNum = solutionNum;
        this.appropriateNeighbour = appropriateNeighbour;
    }

    public long getSolutionsEvaluation() {
        return solutionsEvaluation;
    }

    public int getSolutionNum() {
        return solutionNum;
    }

    public byte[] getAppropriateNeighbour() {
        return appropriateNeighbour;
    }
}
