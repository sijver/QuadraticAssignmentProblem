package main.run;

import main.core.MacroState;
import main.core.MacroStatePre;
import main.core.Qap;
import main.core.SolutionType;
import main.core.utils.QapInstanceReader;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 */
public class Start {

    public static void main(String[] args) {
        Qap qap = QapInstanceReader.readQapInstance("tai11a.dat");

        EnumMap<SolutionType, Double> solutionTypeStatistics = qap.getSolutionTypeStatistics();

        System.out.println(qap.getAllSolutions().size());
        for (SolutionType solutionType : SolutionType.values()) {
            System.out.println(String.format("%1$s %2$.4f%%", solutionType, solutionTypeStatistics.get(solutionType)));
        }

        for(MacroStatePre macroStatePre : qap.getCollapsedSearchLandscapePre()){
            System.out.println(macroStatePre.getSolutionsEvaluation());
            System.out.println("  " + macroStatePre.getSolutionNum() + "   " + Arrays.toString(macroStatePre.getAppropriateNeighbour()));
        }
        System.out.println();
        System.out.println();
        for(MacroState macroState : qap.getCollapsedSearchLandscape()){
            if(macroState.getSolutions().size() == 1 || macroState.getSolutions().size() == 2){
            System.out.println(macroState.getSolutionsEvaluation());
            System.out.println("  " + Arrays.toString(macroState.getSolutions().toArray()));       }
        }
    }


}
