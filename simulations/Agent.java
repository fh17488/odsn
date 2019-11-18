 public class Agent {
	
	private double uncertainty;
	// Each agent is expected to interact (INTERACTIONS / N) * 2 times.  (Each interaction uses 2 agents.)
	// We double the expected length of history for safety.
	private double [] opinion;
	private double uncertanity;
	private double [] uncertanity_log;
	private int opinionCounter;
	private static final boolean logHistory = false;
	
	public Agent(double o, double u, boolean logOpinionDynamicsOfFirstSample) {
		if(Agent.logHistory==true){
			this.opinion = new double[(RAwNetworkModel.INTERACTIONS/RAwNetworkModel.Nm)*100];
			this.uncertanity_log = new double[(RAwNetworkModel.INTERACTIONS/RAwNetworkModel.Nm)*100];
		}
		else{
			this.opinion = new double[1];
		}		
		this.opinionCounter = 0;
		this.setUncertainty(u);
		this.setOpinion(o);
	} // Agent
	
	public double getUncertainty() {
		return this.uncertainty;
	} // getUncertainty
	public double getUncertainty(int i) {
		if(Agent.logHistory == false)
			return this.uncertainty;
		else if (i < opinionCounter)
			return this.uncertanity_log[i];
		else
			return this.uncertanity_log[opinionCounter-1];
	} // getUncertainty
	public void setUncertainty(double u) {
		this.uncertainty = u;
		if(Agent.logHistory==true)
			this.uncertanity_log[opinionCounter] = u;
	} // setUncertainty
	//call setUncertanity before calling setOpinion since setOpinion increments the value of the counter after it sets the opinion. Therefore, if setUncertanity is called before setOpinion it will use the same value of index counter as setOpinion.

	public double getOpinion() {
		if(Agent.logHistory == false)
			return this.opinion[0];
		else
			return this.opinion[opinionCounter-1];
	} // getOpinion	
	public double getOpinion(int i) {
		if(Agent.logHistory == false)
			return this.opinion[0];
		else if (i < opinionCounter)
			return this.opinion[i];
		else
			return this.opinion[opinionCounter-1];
	} // getOpinion
	
	public void setOpinion(double o) {
		if(Agent.logHistory == false)
			this.opinion[0] = o;
		else if (this.opinionCounter < this.opinion.length) {
			this.opinion[opinionCounter] = o;
			this.opinionCounter++;
		} else {
			System.out.println("Error:\tOpinion array not long enough.");
		}
	} // setOpinion

	public int getOpinionCounter(){
		if(Agent.logHistory==false)
			return 0;
		else
			return this.opinionCounter;
	}

} // Agent

