import java.util.Random;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class RAwNetworkModel {

    public static final String networksDirectory = "../Results/networks/";
    public static final String opinionDynamicsDirectory = "../Results/opinion_dynamics/";
    public static final boolean logOpinionDynamicsOfFirstSample = true;
    public static final int SAMPLES = 20; // Number of repetitions at each parameter setting.
    public static final int INTERACTIONS = 4000000;

    public static final double OPINION_MIN = -1.0;
    public static final double OPINION_MAX = 1.0;

    public static final double EXTREMIST_OPINION_BOUNDARY = 0.8;
    public static final double EXTREMIST_UNCERTAINTY = 0.1;

    public static final double MU = 0.2; // How much agents affect each other.

    public static final double Pe_MAX = 0.3; // Maximum proportion of Extremists.
    public static final double U_MAX = 2.0; // Maximum uncertainty.
    //public static final int Pe_STEPS = 50; // Increments of Pe to test.
    //public static final int U_STEPS = 10; // Increments of uncertainty to test.

    public static final int Nm = 1000; // Total number of moderate agents
    //private static final double delta = 0.00; //The code effectively uses this value for delta.
    private String networkTopology = null;
    int[][] network = null;
    private List<int[]> neighbourhoods = null;
    public int N = 0; // Total number of agents
    private double Pe = 0;
    private int tot_extremists = 0;
    private int Nx1 = 0;
    private int Nx2 = 0;
    private double U = 0; 
    private boolean initializationSuccessful = false;
    private double avgY = 0;
    private double sdY = 0;
    private List<String> messageLog = null;
    private String messagesLogFile = null;
    private String metricYAvgLogFile = null;
    private String metricYStdDevLogFile = null;
    private String modelParamsLogFile = null;
    private String opinionDynamicsLogFile = null;
    
    //Action functions
    public static void main (String [] args) {      
        RAwNetworkModel model = new RAwNetworkModel(String.valueOf(args[0]), Double.valueOf(args[1]), Double.valueOf(args[2]));        
        if(model.isInitializationSuccessful() == true){
            model.logModelParamsToFile();
            long startTime = System.nanoTime();
            model.simulateScenario();
            long endTime = System.nanoTime();
            long duration = (endTime - startTime)/1000000000;  //divided by 1000000000 to get seconds.
            model.logAvgY();
            model.logStdDevY();
            model.logMessage("Runtime:\t"+duration+"sec");
        }
        model.exportMessageLogToFile();
        return;
    }
    public RAwNetworkModel(String networkTopology, double Pe, double U){        
        this.networkTopology = networkTopology;
        this.Pe = Pe;
        this.U = U;
        this.messageLog = new ArrayList<>();            
        if(areParametersValid() == true){
            this.metricYAvgLogFile = this.networkTopology+"_network_"+Double.valueOf(this.Pe)+"_"+Double.valueOf(this.U)+"_MetricYAvgLog.txt";
            this.metricYStdDevLogFile = this.networkTopology+"_network_"+Double.valueOf(this.Pe)+"_"+Double.valueOf(this.U)+"_MetricYStdDevLog.txt";
            this.modelParamsLogFile = this.networkTopology+"_network_"+Double.valueOf(this.Pe)+"_"+Double.valueOf(this.U)+"_ModelParamsLog.txt";
            this.opinionDynamicsLogFile = this.networkTopology+"_network_"+Double.valueOf(this.Pe)+"_"+Double.valueOf(this.U)+"_OpinionDynamicsLog";//csv ext. added in function.
            this.messagesLogFile = this.networkTopology+"_network_"+Double.valueOf(this.Pe)+"_"+Double.valueOf(this.U)+"_MessageLog.txt";

            String networkCSVFile = this.networkTopology+"_network_"+Double.valueOf(this.Pe)+".csv";
            CsvParser csvParser = new CsvParser(RAwNetworkModel.networksDirectory, networkCSVFile);
            this.network = csvParser.getNetwork();
            if(isNetworkValid() == false)
                return;
            this.N = this.network.length;
            this.tot_extremists = this.N - RAwNetworkModel.Nm;
            this.Nx1 = this.tot_extremists/2;
            this.Nx2 = this.tot_extremists/2;   

            this.neighbourhoods = new ArrayList<>();
            for(int i=0; i < this.N; i++){
                int [] neighbours_unfiltered = new int[this.N];
                int neighbour_count = 0;
                for(int j=0; j < this.N; j++){
                    if(this.network[i][j] == 1){
                        neighbours_unfiltered[neighbour_count] = j;
                        neighbour_count++;
                    }
                }
                int [] neighbours = new int[neighbour_count];
                for(int ii=0; ii<neighbour_count; ii++){
                    neighbours[ii] = neighbours_unfiltered[ii];
                }
                this.neighbourhoods.add(neighbours);
            }
            this.logNeighbourhoods();
            
            this.initializationSuccessful = true;    
        }
        return;
    }//TBD: ensure the neighbourhoods datastructure is valid.
    public void simulateScenario(){
        Random random = new Random();
        double [] ySample = new double[RAwNetworkModel.SAMPLES];
        for (int sample = 0; sample < RAwNetworkModel.SAMPLES; sample++) {
            int totalInteractions;            
            double yMetric = -1;
            double yMetricNewValue = -1;
            Agent [] agent = new Agent[this.N];
            agent = this.initializeAgents(agent);
            // Agents interact.
            for (totalInteractions = 0; totalInteractions < RAwNetworkModel.INTERACTIONS; totalInteractions++) {
                // Choose 2 random agents.
                int agenti = -1, agentj = -1;
                int[] neighbours = null;
                int neighbour_count = 0;
                while(neighbour_count == 0){
                    agenti = random.nextInt(this.N);
                    neighbours = neighbourhoods.get(agenti);
                    neighbour_count = neighbours.length;
                }
                agentj = random.nextInt(neighbour_count);
                agent = interact(agent, agenti, agentj);                
            } // interaction for loop end.

            ySample[sample] = calculateY(agent);
            if(sample==0 && logOpinionDynamicsOfFirstSample){      
                this.logMessage("Y value for which plot was made is:\t"+ySample[sample]);
                this.logOpinionDynamicsToFile(agent, totalInteractions, sample);
            }            
        } // sample for loop end.

        // Calculate and save average y for those parameters.
        Statistics calculator = new Statistics(ySample);
        this.avgY = calculator.getMean();
        this.sdY = calculator.getStdDev();

        return;
    }//TBD: Check function
    private Agent [] initializeAgents(Agent [] agent){
        // Initialise agents (moderates then extremists).
        Random random = new Random();
        double randomOpinion;
        for (int i = 0; i < RAwNetworkModel.Nm; i++) {
            randomOpinion = (random.nextDouble() * (RAwNetworkModel.EXTREMIST_OPINION_BOUNDARY * 2)) - RAwNetworkModel.EXTREMIST_OPINION_BOUNDARY;
            agent[i] = new Agent(randomOpinion, this.U, RAwNetworkModel.logOpinionDynamicsOfFirstSample);
        }
        for(int i = RAwNetworkModel.Nm; i < RAwNetworkModel.Nm + this.Nx1; i++){
            randomOpinion = (random.nextDouble() * (RAwNetworkModel.OPINION_MAX - RAwNetworkModel.EXTREMIST_OPINION_BOUNDARY)) + RAwNetworkModel.EXTREMIST_OPINION_BOUNDARY;
            agent[i] = new Agent(randomOpinion, RAwNetworkModel.EXTREMIST_UNCERTAINTY, RAwNetworkModel.logOpinionDynamicsOfFirstSample);
        }
        for (int i = RAwNetworkModel.Nm + this.Nx1; i < this.N; i++) {
            randomOpinion = (random.nextDouble() * (RAwNetworkModel.OPINION_MAX - RAwNetworkModel.EXTREMIST_OPINION_BOUNDARY)) + RAwNetworkModel.EXTREMIST_OPINION_BOUNDARY;
            agent[i] = new Agent(-randomOpinion, RAwNetworkModel.EXTREMIST_UNCERTAINTY, RAwNetworkModel.logOpinionDynamicsOfFirstSample);
        }
        return agent;
    }
    private Agent [] interact(Agent [] agent, int agenti, int agentj){
        // Interact!
        double Xi, Xj, Ui, Uj;
        Xi = agent[agenti].getOpinion();
        Xj = agent[agentj].getOpinion();
        Ui = agent[agenti].getUncertainty();
        Uj = agent[agentj].getUncertainty();

        double Hji = Math.min(Xj + Uj, Xi + Ui) - Math.max(Xj - Uj, Xi - Ui);
        double Hij = Math.min(Xi + Ui, Xj + Uj) - Math.max(Xi - Ui, Xj - Uj);

        double RAji = (Hji / Uj) - 1.0;
        double RAij = (Hij / Ui) - 1.0;

        // Update
        if (Hji > Uj) {
            agent[agenti].setUncertainty( Ui + (RAwNetworkModel.MU * RAji * (Uj - Ui)) );
            agent[agenti].setOpinion( Xi + (RAwNetworkModel.MU * RAji * (Xj - Xi)) );
        } else {
            agent[agenti].setUncertainty(Ui); // For history.
            agent[agenti].setOpinion(Xi); // For history.
        }
        if (Hij > Ui) {
            agent[agentj].setUncertainty( Uj + (RAwNetworkModel.MU * RAij * (Ui - Uj)) );
            agent[agentj].setOpinion( Xj + (RAwNetworkModel.MU * RAij * (Xi - Xj)) );
        } else {
            agent[agentj].setUncertainty(Uj); // For history.
            agent[agentj].setOpinion(Xj); // For history.
        }
        return agent;
    }

    //Helper functions
    public boolean isInitializationSuccessful(){
        return this.initializationSuccessful;
    }
    public void logMessage(String message){
        this.messageLog.add(message);
        return;
    }
    private boolean areParametersValid(){
        boolean validStatus = true;
        if(this.networkTopology.equals("ring-lattice")==false && this.networkTopology.equals("small-world")==false && this.networkTopology.equals("random")== false){
            this.messageLog.add("Abort:\tThis is an invalid network topology.");
            validStatus = false;
        }
        if(this.Pe<0 || this.Pe>RAwNetworkModel.Pe_MAX){
            this.messageLog.add(String.format("Abort:\tThe value of Pe=%.3f is not in the range [0,%.3f].",this.Pe,RAwNetworkModel.Pe_MAX));
            validStatus = false;
        }
        if(this.U<0 || this.U>2){
            this.messageLog.add(String.format("Abort:\tThe value of U=%.3f is not in the range [0,2].",this.U));
            validStatus = false;
        }
        return validStatus;
    }
    private boolean isNetworkValid(){
        boolean validStatus = true;
        int N = RAwNetworkModel.Nm + (int)Math.round(this.Pe*RAwNetworkModel.Nm);
        if(N%2!=0)
            N++;
        int rows = this.network.length;
        if(rows != N){//check rows=n
            this.messageLog.add("Abort:\tThe number of rows="+rows+" in the CSV do not equal Nm+(Pe*Nm)="+N);
            validStatus = false;
        }
        for(int i=0; i<N; i++){//check all values = {0,1}
            for(int j=0; j<N; j++){
                if(this.network[i][j] != 0 && this.network[i][j] != 1){
                    this.messageLog.add("Abort:\tThe value at index (row,col)=("+i+","+j+") is "+network[i][j]+" and not in {0,1}.");
                    validStatus = false;
                }
            }
        } 
        return validStatus;
    }
    private double calculateY(Agent[] agent){
        double yMetricValue = 0;
        int newPositiveExtremists = 0;
        int newNegativeExtremists = 0;
        for (int i = 0; i < RAwNetworkModel.Nm; i++) { // Don't count original extremists.
            if (agent[i].getOpinion() > RAwNetworkModel.EXTREMIST_OPINION_BOUNDARY-0.1)
                newPositiveExtremists++;
            if (agent[i].getOpinion() < -RAwNetworkModel.EXTREMIST_OPINION_BOUNDARY+0.1)
                newNegativeExtremists++;
        }
        double pplus = (double)newPositiveExtremists / (double)(RAwNetworkModel.Nm);
        double pminus = (double)newNegativeExtremists / (double)(RAwNetworkModel.Nm);
        yMetricValue = (pplus * pplus) + (pminus * pminus);
        return yMetricValue;
    }
    /*private boolean continueSimulation(int totalInteractions, double yMetric, double yMetricNewValue){
        boolean continueSimulation = true;
        if(Math.abs(yMetricNewValue - yMetric) > 0.005)
            this.logMessage(String.format("The new value of metric y:\t%.5f", yMetricNewValue));
        else if(totalInteractions >= (RAwNetworkModel.INTERACTIONS_PER_AGENT*this.N)/2){
            this.logMessage(String.format("Simulation Ended:\tThe number of interactions has reached a limit that is:\t%d", ((RAwNetworkModel.INTERACTIONS_PER_AGENT*this.N)/2)));
            continueSimulation = false;
        }
        else{
            this.logMessage(String.format("The final value of metric y:\t%.5f", yMetricNewValue));
            continueSimulation = false;
        }
        return continueSimulation;
    }*/

    //Log functions
    public void logModelParamsToFile(){
        FileWriter fileWriter = null;
        PrintWriter printWriter = null;
        try{
           fileWriter = new FileWriter(RAwNetworkModel.opinionDynamicsDirectory+this.modelParamsLogFile, true);
           printWriter = new PrintWriter(fileWriter);            
           printWriter.printf("N=%d\t", this.N);
           printWriter.printf("Nm=%d\t", RAwNetworkModel.Nm);
           printWriter.printf("pe=%.4f\t", this.Pe);
           printWriter.printf("tot_extremists=%d\t", this.tot_extremists);
           printWriter.printf("Nx1=%d\t", this.Nx1);
           printWriter.printf("Nx2=%d\t", this.Nx2);
           printWriter.printf("ue=%.2f\t", RAwNetworkModel.EXTREMIST_UNCERTAINTY);
           printWriter.printf("U=%.2f\t", this.U);
           printWriter.printf("MU=%.2f\t", RAwNetworkModel.MU);
           printWriter.printf("topology=%s\t", this.networkTopology);
           printWriter.close();
           fileWriter.close();
        }
        catch(IOException excep){
           excep.printStackTrace();
        }
        return;
    }
    public void logOpinionDynamicsToFile(Agent[] agent, int totalInteractions, int sampleID){
        //Export detailed logs of the opinion dynamics process for a single scenario.
        try {
            FileWriter file = new FileWriter(RAwNetworkModel.opinionDynamicsDirectory+this.opinionDynamicsLogFile+"_"+sampleID+".csv");
            BufferedWriter buffer = new BufferedWriter(file);

            // Format and output data.
            // agent[0].opinion[0] \t agent[0].uncertanity[0] \t agent[i].opinion[0] \t agent[i].uncertanity[0] \t agent[N-1].opinion[0] \t agent[N-1].uncertanity[0]
            // agent[0].opinion[i] \t agent[0].uncertanity[i] \t agent[i].opinion[i] \t agent[i].uncertanity[i] \t agent[N-1].opinion[i] \t agent[N-1].uncertanity[i]
            // agent[0].opinion[(INTERACTIONS/N)] \t agent[0].uncertanity[(INTERACTIONS/N)] \t agent[i].opinion[(INTERACTIONS/N)] \t agent[i].uncertanity[(INTERACTIONS/N)] \t agent[N-1].opinion[(INTERACTIONS/N)] \t agent[N-1].uncertanity[(INTERACTIONS/N)]

            //System.out.println("Total number of interactions:\t"+totalInteractions);
	        //int totalInteractionsPerAgent = totalInteractions/this.N;
            //for (int iteration_idx = 0; iteration_idx < totalInteractionsPerAgent; iteration_idx++) {
                for (int agent_idx = 0; agent_idx < this.N; agent_idx++){
                    buffer.write("\t" + String.format("%.5f",agent[agent_idx].getOpinion()));
                    buffer.write("\t" + String.format("%.5f",agent[agent_idx].getUncertainty()));
                }
                //buffer.newLine();
            //}
            buffer.close();
	        file.close();
        } catch (IOException e) {
            System.out.println("Error:\tFile writing failed");
        }
        return;
    }//TBD: Check function
    public void logAvgY(){
        this.logMetricToFile(this.metricYAvgLogFile, this.avgY);
        return;
    }
    public void logStdDevY(){
        this.logMetricToFile(this.metricYStdDevLogFile, this.sdY);
        return;
    }
    public void logMetricToFile(String fileName, double metricValue){
       FileWriter fileWriter = null;
       PrintWriter printWriter = null;
       try{
          fileWriter = new FileWriter(RAwNetworkModel.opinionDynamicsDirectory+fileName);
          printWriter = new PrintWriter(fileWriter); 
          printWriter.printf("%.3f", metricValue);
          printWriter.close();
          fileWriter.close();
       }
       catch(IOException excep){
          excep.printStackTrace();
       }
       return;
    }
    private void logNeighbourhoods(){
        List<Integer> moderates = new ArrayList<>();
        List<Integer> positiveExtremists = new ArrayList<>();
        List<Integer> negativeExtremists = new ArrayList<>();
        int agentIdx = 0;
        for(int[] neighbourhood: neighbourhoods){
            int numNeighbours = neighbourhood.length;
            boolean foundCount = false;
            if(agentIdx < RAwNetworkModel.Nm){//moderate
                for(int i=0; i < moderates.size(); i++){
                    if(moderates.get(i).intValue()==numNeighbours)
                        foundCount = true;
                }
                if(foundCount == false){                        
                    moderates.add(Integer.valueOf(numNeighbours));
                }
            }
            else if(agentIdx < RAwNetworkModel.Nm+this.Nx1){//positive extremist
                for(int i=0; i < positiveExtremists.size(); i++){
                    if(positiveExtremists.get(i).intValue()==numNeighbours)
                        foundCount = true;
                }
                if(foundCount == false){                        
                    positiveExtremists.add(Integer.valueOf(numNeighbours));
                }
            }
            else{//negative extremist
                for(int i=0; i < negativeExtremists.size(); i++){
                    if(negativeExtremists.get(i).intValue()==numNeighbours)
                        foundCount = true;
                }
                if(foundCount == false){                        
                    negativeExtremists.add(Integer.valueOf(numNeighbours));
                }
            }
            agentIdx++;
        }
        FileWriter fileWriter = null;
        BufferedWriter buffer = null;
        try{
           fileWriter = new FileWriter(RAwNetworkModel.opinionDynamicsDirectory+this.modelParamsLogFile);
           buffer = new BufferedWriter(fileWriter);            
           buffer.write("Distinct neighbour count of moderates:");
           for(Integer neighbourCount: moderates){
                buffer.write("\t" + neighbourCount.intValue());
           }
           buffer.newLine();
           buffer.write("Distinct neighbour count of positive extremists:");
           for(Integer neighbourCount: positiveExtremists){
                buffer.write("\t" + neighbourCount.intValue());
           }
           buffer.newLine();
           buffer.write("Distinct neighbour count of negative extremists:");
           for(Integer neighbourCount: negativeExtremists){
                buffer.write("\t" + neighbourCount.intValue());
           }
           buffer.newLine();
           buffer.close();
           fileWriter.close();
        }
        catch(IOException excep){
           excep.printStackTrace();
        }
        return;
    }
    public void exportMessageLogToFile(){
        FileWriter fileWriter = null;
        PrintWriter printWriter = null;
        if(this.messageLog.size()>0){
            try{
                fileWriter = new FileWriter(RAwNetworkModel.opinionDynamicsDirectory+this.messagesLogFile);
                printWriter = new PrintWriter(fileWriter);
                int messageCount = 1;
                for(String message: this.messageLog) {
                        printWriter.printf("%d:\t%s\n", messageCount++, message);
                }
                printWriter.close();
                fileWriter.close();
            }
            catch(IOException excep){
                excep.printStackTrace();
            }
        }
        return;
    }
} // RAwNetworkModel

