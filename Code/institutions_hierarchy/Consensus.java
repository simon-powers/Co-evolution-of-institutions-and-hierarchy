package institutions_hierarchy;

public class Consensus {
	private double h;
	private int timeToConsensus;
	
	
	
	public Consensus(double h, int timeToConsensus) {
		this.h = h;
		this.timeToConsensus = timeToConsensus;
	}
	public double getH() {
		return h;
	}
	public void setH(double h) {
		this.h = h;
	}
	public int getTimeToConsensus() {
		return timeToConsensus;
	}
	public void setTimeToConsensus(int timeToConsensus) {
		this.timeToConsensus = timeToConsensus;
	}
	
	
}
