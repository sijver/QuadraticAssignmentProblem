package main.core.utils;

import main.core.Qap;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 */
public class QapInstanceReader {

    public static Qap readQapInstance(String filePath) {
        byte instanceSize = 0;
        int[][] dMatrix = null;
        int[][] fMatrix = null;

        try {
            FileInputStream fileInputStream = new FileInputStream(filePath);
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(fileInputStream));
            String strLine;

            Pattern intPattern = Pattern.compile("\\d+");
            Matcher matcher;

            boolean firstLine = true;
            boolean isDMatrixReading = true;
            int rowCounter = 0;

            // Read file line by line
            while ((strLine = bufferedReader.readLine()) != null) {
                matcher = intPattern.matcher(strLine);
                // Read instance size information
                if (firstLine) {
                    if (matcher.find()) {
                        instanceSize = (byte) Integer.parseInt(matcher.group());
                        dMatrix = new int[instanceSize][instanceSize];
                        fMatrix = new int[instanceSize][instanceSize];
                    }
                    firstLine = false;
                    // Reading matrices
                } else if (strLine.length() > 0) {
                    for (int i = 0; i < instanceSize; i++) {
                        if (matcher.find()) {
                            if(isDMatrixReading){
                                dMatrix[rowCounter][i] = Integer.parseInt(matcher.group());
                            } else {
                                fMatrix[rowCounter][i] = Integer.parseInt(matcher.group());
                            }
                        }
                    }
                    rowCounter++;
                    if (rowCounter == instanceSize) {
                        rowCounter = 0;
                        isDMatrixReading = false;
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

        return new Qap(instanceSize, dMatrix, fMatrix);
    }

    public static void readBestKnown(String filePath){
        //TODO
    }

}
