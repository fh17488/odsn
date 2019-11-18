import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

public class CsvParser {

    int[][] network;
    int tot_agents = 0;

    public static void main (String []args){
        String directory = "../Results/networks/";
        //String fileName = "ring-lattice_network_0.025.csv";
        String fileName = args[0];

        CsvParser csvParser = new CsvParser(directory, fileName);
        csvParser.printNetwork();
        return;
    }

    public CsvParser(String directory, String fileName) {        
        File file= new File(directory + fileName);

        // this gives you a 2-dimensional array of strings
        List<List<String>> lines = new ArrayList<>();
        Scanner inputStream;

        try{
            inputStream = new Scanner(file);

            while(inputStream.hasNext()){
                String line= inputStream.next();
                String[] values = line.split(",");
                // this adds the currently parsed line to the 2-dimensional string array
                lines.add(Arrays.asList(values));
            }

            inputStream.close();
        }catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        this.tot_agents = lines.size() - 1;
        this.network = new int[this.tot_agents][this.tot_agents];
        for(int i = 0; i < this.tot_agents; i++)
            for(int j = 0; j < this.tot_agents; j++)
                this.network[i][j] = -1;

        // the following code lets you iterate through the 2-dimensional array
        int lineNo = 1;
        for(List<String> line: lines) {
            int columnNo = 1;
            if(lineNo > 1){
                for (String value: line) {                    
                    if(columnNo > 1){
                        int rowIdx = lineNo - 2;
                        int colIdx = columnNo - 2;
                        this.network[rowIdx][colIdx] = Integer.valueOf(value);
                        //System.out.println("Line " + lineNo + " Column " + columnNo + ": " + value);
                    }
                    columnNo++;
                }
            }
            lineNo++;
        }
    }

    public void printNetwork()
    {
        for(int i = 0; i < this.tot_agents; i++){
            for(int j = 0; j < this.tot_agents; j++){
                System.out.print(this.network[i][j]);
                System.out.print('\t');
            }
            System.out.println();
        }
        return;
    }

    int getTotAgents(){
        return this.tot_agents;
    }

    int[][] getNetwork(){
        return this.network;
    }
}
