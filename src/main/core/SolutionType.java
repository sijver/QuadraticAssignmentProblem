package main.core;

/**
 * Created with IntelliJ IDEA.
 */
public enum SolutionType {
    SLMIN("strict local minima"),
    LMIN("local minima"),
    IPLAT("interior plateau"),
    LEDGE("ledge"),
    SLOPE("slope"),
    LMAX("local maxima"),
    SLMAX("strict local maxima");

    private SolutionType(String fullName) {
        this.fullName = fullName;
    }

    private String fullName;

    public String getFullName() {
        return fullName;
    }
}
