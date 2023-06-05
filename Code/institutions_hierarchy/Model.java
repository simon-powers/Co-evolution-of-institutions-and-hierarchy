package institutions_hierarchy;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.util.ArrayList;
import cern.jet.random.Poisson;
import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister;
import java.util.Date;

/**
 *
 * @author spowers
 */
public class Model {
    
    private int NUM_DEMES;
    private int T;
    private int kBase;
    private double instCostScale;
    private double MUT_RATE;
    private double MIG_RATE;
    private double H_MUT_RATE;
    private double beta;
    private double grad;
    private double coopB;
    private double coopC;
    
    //consensus formation parameters
    private int numListeners;
    private double kScale;
    private double xThr; 
    
    private double rBase=2;
        
    public static final int C=0;
    public static final int S=1;
    public static final int L=2;
    public static final int NUM_TYPES = 3;
        
    private final double MUT_SIZE=0.1;
    private final double a_gl=0.05;
    private final double a_lg=0.05;
        
    //private double[] grandCoopMean = new double[T];
   // private double[] grandSelfishMean = new double[T];
    //private double[] grandLonerMean = new double[T];
    private double[] meanH;
    private double[][] globalProps;
    private double[] kGroupMean;
    private int[] groupMeanConsensusTime;
    private double[] alphaMean;
    private double[] instCostMean;
    private int[] numLeaders;
    private double[] alphaSkewness;
    private double[][] groupHRecord; //record of agreed H on each patch at each generation
    private double[][] skewnessRecord;
    private ArrayList<Individual>[] demes;
    private ArrayList<Individual>[] migrantPool;
    private double[] propLeaders;
    private Date seed;
    private MersenneTwister generator;
    private Normal normal;
    
    private ArrayList<String> indivRecord;
    private int printEvery;
    private boolean evolveAlpha;

    public Model(int NUM_DEMES,int T, int kBase, double instCostScale,
                        double MUT_RATE, double MIG_RATE,
                        double beta, double grad, double coopB, double coopC,
                        int numListeners, double kScale, double xThr, double H_MUT_RATE,
                        int printEvery, boolean evolveAlpha){
        this.NUM_DEMES=NUM_DEMES;
        this.T=T;
        this.kBase=kBase;
        this.instCostScale=instCostScale;
        this.MUT_RATE=MUT_RATE;
        this.MIG_RATE=MIG_RATE;
        this.H_MUT_RATE = H_MUT_RATE;
        this.beta=beta;
        this.grad=grad;
        this.coopB=coopB;
        this.coopC=coopC;
        
        this.numListeners = numListeners;
        this.kScale = kScale;
        this.xThr = xThr;
        
        this.demes = new ArrayList[NUM_DEMES];
        this.migrantPool = new ArrayList[NUM_DEMES];
        this.meanH = new double[T];
        this.kGroupMean = new double[T];
        this.globalProps = new double[T][NUM_TYPES];
        this.groupHRecord = new double[T][NUM_DEMES];
        this.skewnessRecord = new double[T][NUM_DEMES];
        this.groupMeanConsensusTime = new int[T];
        this.alphaMean = new double[T];
        this.instCostMean = new double[T];
        this.numLeaders = new int[T];
        this.propLeaders = new double[T];
        this.alphaSkewness = new double[T];
        this.indivRecord = new ArrayList<>();
        this.printEvery = printEvery;
        this.evolveAlpha = evolveAlpha;
        //initialise random number generator
        this.seed = new Date();
        this.generator = new MersenneTwister(seed);
        //this.generator = new MersenneTwister(1);
        double normMean=0.0;
        this.normal=new Normal(normMean, MUT_SIZE, generator);
                
        
    }
    
    private void initDemes(){
        //initialise demes full of loners with random hValues, using default 
        //constructor
    	 double alphaLeader = 0.99;
         double alphaFollower = 0.99;
        for(int i=0; i<NUM_DEMES; i++){
            ArrayList<Individual> a = new ArrayList<Individual>(kBase);
            for(int j=0; j<kBase; j++){
            	//adding a leader
            	/*if(j == 0){
            		a.add(new Individual(generator, normal, L, generator.nextDouble(), alphaLeader));
            	}
            	else{
            		a.add(new Individual(generator, normal, L, generator.nextDouble(), alphaFollower));
            	}*/
            	if(evolveAlpha){
            		a.add(new Individual(generator, normal, L, generator.nextDouble(), generator.nextDouble()));
            
            	}
            	else{
            		a.add(new Individual(generator, normal, L, generator.nextDouble(), 0.5));
            	}
            	
            	
            }
            demes[i]=a;
        }
        
    }
    
    private void initDemesAllSocial(){
        //initialise demes full of loners with random hValues, using default 
        //constructor
        for(int i=0; i<NUM_DEMES; i++){
            ArrayList<Individual> a = new ArrayList<Individual>(kBase);
            double alphaLeader = 0.9;
            double alphaFollower = 0.1;
            for(int j=0; j<kBase; j++){
            	//adding a leader
            	if(j == 0){
            		a.add(new Individual(generator, normal, S, generator.nextDouble(), alphaLeader));
            	}
            	else{
            		a.add(new Individual(generator, normal, S, generator.nextDouble(), alphaFollower));
            	}
                
            }
            demes[i]=a;
        }
    }
    
    private double[] localProps(ArrayList<Individual> deme){
        //return local proportions within deme in a 1D array
        double[] props = new double[this.NUM_TYPES];
        int numDemeCoop=0;
        int numDemeSelfish=0;
        int numDemeLoner=0;
        int demeSize=0;
        Individual indiv;
        for(int i=0; i<deme.size(); i++){
            indiv = deme.get(i);
            if(indiv.getType()==this.C){
                numDemeCoop++;
            } else if (indiv.getType()==this.S){
                numDemeSelfish++;
            } else {
                numDemeLoner++;
            }
           
        }
         demeSize=numDemeCoop+numDemeSelfish+numDemeLoner;
            props[this.C]= (double) numDemeCoop/demeSize;
            //have to cast to double to avoid integer division.
            props[this.S] = (double) numDemeSelfish/demeSize;
            props[this.L] = (double) numDemeLoner/demeSize;
            return props;
        
                
    }
    
    private int[] localTypeNums(ArrayList<Individual> deme){
        //return local proportions within deme in a 1D array
        int[] typeNums = new int[this.NUM_TYPES];
        int numDemeCoop=0;
        int numDemeSelfish=0;
        int numDemeLoner=0;
        //int demeSize=0;
        Individual indiv;
        for(int i=0; i<deme.size(); i++){
            indiv = deme.get(i);
            if(indiv.getType()==Model.C){
                numDemeCoop++;
            } else if (indiv.getType()==Model.S){
                numDemeSelfish++;
            } else {
                numDemeLoner++;
            }
           
        }
         //demeSize=numDemeCoop+numDemeSelfish+numDemeLoner;
            typeNums[this.C]= numDemeCoop;
            //have to cast to double to avoid integer division.
            typeNums[this.S] = numDemeSelfish;
            typeNums[this.L] = numDemeLoner;
            return typeNums;
        
                
    }
    
    private double[] globalProps(){
        double[] props = new double[this.NUM_TYPES];
        int numCoop = 0;
        int numSelfish = 0;
        int numLoner = 0;
        //int popSize=0;
        ArrayList<Individual> deme;
        Individual indiv;
        for (int i = 0; i < NUM_DEMES; i++) {
            deme = demes[i];
            for (int j = 0; j < deme.size(); j++) {
                indiv = deme.get(j);
                if (indiv.getType() == this.C) {
                    numCoop++;
                } else if (indiv.getType() == this.S) {
                    numSelfish++;
                } else {
                    numLoner++;
                }
            }
        }
        int popSize = numCoop + numSelfish + numLoner;
        props[this.C] = (double) numCoop / popSize;
        //have to cast to double to avoid integer division.
        props[this.S] = (double) numSelfish / popSize;
        props[this.L] = (double) numLoner / popSize;
        return props;

    }
    
    private double meanH(ArrayList<Individual> deme){
        //return mean Hvalue in deme (not of loners).
        double hTotal = 0;
        Individual indiv;
        int numNotLoner=0;
        for(int i=0; i<deme.size(); i++){
            indiv=deme.get(i);
            if(indiv.getType()!=this.L){
                hTotal += indiv.getHValue();
                numNotLoner++;
            }
        }
        if(numNotLoner > 0){
        	return hTotal/numNotLoner;
        }
        else{
        	return 0;
        }
        
    }
    
    private Consensus setHByConsensus(ArrayList<Individual> demeAll)
    {
    	
    	for(Individual indiv : demeAll){
    		indiv.setHBehaviourTime(indiv.getHValue());
    	}
    	
    	//Make deme without asocials
    	ArrayList<Individual> deme = new ArrayList<>();
    	for(Individual indiv : demeAll){
    		if(indiv.getType() != L){
    			deme.add(indiv); //shallow copy!
    		}
    		
    	}
    	//System.out.println("Deme size " + deme.size());
    	if(deme.size() == 0){
    		return new Consensus(0,0);
    	}
    	else if(deme.size() == 1){
    		return new Consensus(deme.get(0).getHValue(), 0);
    	}
    	
    	double sdX; double sdXLead; int nEvent; int speaker; int numListeners2 = this.numListeners; String resTemp ="";
    	double xMax = 1;
    	double alphaSum; 
    	double diffAlpha; double newFSpeaker;                                                       //CHANGE
    	int[] listenerList;

    	if(deme.size() <= this.numListeners){numListeners2 = deme.size()-1;}  

    	sdX = Utility.sdPref(deme);                                    //Initial variance of f
    	nEvent = 0;

    	alphaSum = 0.0d;                                                 // Sum of the alpha to calculate the probability of being chosen as a speaker
    	double[] probSpeaker = new double[deme.size()];                     //Array of the probability of each ConsensusIndividual
    	for(int i=0; i<deme.size(); i++){                       
    		alphaSum += Math.pow(deme.get(i).getAlpha(),kScale);
    	}
    	double probTemp = 0;                
    	for(int i=0; i<deme.size(); i++){                                   //We calculate the weigthed probabilities
    		probSpeaker[i]= probTemp + (Math.pow(deme.get(i).getAlpha(),kScale)/alphaSum);
    		probTemp = probSpeaker[i];
    	}                                                   
    	//Simulation of consensus decision making
    	while(sdX > xThr){

    		speaker = Utility.probSample(probSpeaker, Math.random());         //Sample a speaker as a function of alpha

    		//Sample a random speaker
    		listenerList = Utility.randomSampleOtherList(deme.size(), numListeners2, speaker);              //Create the list of listener (random ConsensusIndividuals except speaker)
    		for(int i=0; i<numListeners2; i++){
    			//ConsensusMod1 Only comparison
    			diffAlpha =  deme.get(speaker).getAlpha() - deme.get(listenerList[i]).getAlpha();
    			if(diffAlpha <= 0.01){diffAlpha = 0.01;} //can set 0.05 here for faster convergence
    			deme.get(listenerList[i]).setHBehaviourTime(deme.get(listenerList[i]).getHBehaviourTime() + diffAlpha * (deme.get(speaker).getHBehaviourTime()-deme.get(listenerList[i]).getHBehaviourTime()));

    		}
    		sdX = Utility.sdPref(deme);

    		//Counter of time of consensus
    		nEvent++;

    	}
    	//now take mean of hBehaviourTime
    	double hTotal = 0;
    	for(Individual indiv : deme){
    		hTotal += indiv.getHBehaviourTime();
    	} 
    	double hConsensus = hTotal / deme.size();

    	return new Consensus(hConsensus, nEvent);

    }
    
    
    private double meanHAllIndivs(){
    	double hTotal = 0;
    	int numIndivs = 0;
    	for(ArrayList<Individual> deme : demes){
    		
    		for(Individual indiv : deme){
                hTotal += indiv.getHValue();
                numIndivs++;
    		} 
        }
    	return hTotal / numIndivs;
    }
    
    private double meanAlpha(){
    	double alphaTotal = 0;
    	int numIndivs = 0;
    	for(ArrayList<Individual> deme : demes){
    		
    		for(Individual indiv : deme){
                alphaTotal += indiv.getAlpha();
                numIndivs++;
    		} 
        }
    	return alphaTotal / numIndivs;
    }
    
    private double skewAlpha(){
    	ArrayList<Individual> population = new ArrayList<>();
    	for(ArrayList<Individual> deme : demes){
    		population.addAll(deme);
    	}
    	return Utility.skewnessAlpha(population);
    }

    private int numLeaders(){
    	int numLeaders = 0;
    	for(ArrayList<Individual> deme : demes){
    		for(Individual indiv : deme){
    			if(indiv.getAlpha() > 0.9){
    				numLeaders++;
    			}
    		} 
    	}
    	return numLeaders;
    }
    
    private int popSize(){
    	int numIndivs = 0;
    	for(ArrayList<Individual> deme : demes){
    		numIndivs += deme.size();
    	}
    	return numIndivs;
    }
    
    private ArrayList<Individual> generateOffspring(ArrayList<Individual> deme,
            double lambda_coop, double lambda_selfish, double lambda_loner ){
        
        Poisson coopPoisson = new Poisson(lambda_coop, generator);
        Poisson selfishPoisson = new Poisson(lambda_selfish, generator);
        Poisson lonerPoisson = new Poisson(lambda_loner, generator);
        int numAdults=deme.size();
        int[] typeNums = localTypeNums(deme);
        ArrayList<Individual> offspring =new ArrayList<Individual>(numAdults*2);
        Individual parent;
        Individual child;
        int parentType;
        int numIndivOffspring;
        //generate offspring from poisson distribution 
        for(int i=0; i<numAdults; i++){
            parent = deme.get(i);
            parentType = parent.getType();
            if(parentType == this.C){
                numIndivOffspring = coopPoisson.nextInt();
                for(int j=0; j<numIndivOffspring; j++){
                    child = new Individual(generator, normal, 
                            parentType, parent.getHValue(), parent.getAlpha());
                    if(generator.nextDouble()<H_MUT_RATE)
                        child.mutateH();
                    if(generator.nextDouble()<MUT_RATE)
                        child.mutateType();
                    if(generator.nextDouble()<MUT_RATE && evolveAlpha)
                        child.mutateAlpha();
                    offspring.add(child);
                }
            }
            else if(parentType == this.S){
                numIndivOffspring = selfishPoisson.nextInt();
                for(int j=0; j<numIndivOffspring; j++){
                    child = new Individual(generator, normal, 
                            parentType, parent.getHValue(), parent.getAlpha());
                    if(generator.nextDouble()<H_MUT_RATE)
                        child.mutateH();
                    if(generator.nextDouble()<MUT_RATE)
                        child.mutateType();
                    if(generator.nextDouble()<MUT_RATE && evolveAlpha)
                        child.mutateAlpha();
                    offspring.add(child);
                }
            }
            else{
                numIndivOffspring = lonerPoisson.nextInt();
                for(int j=0; j<numIndivOffspring; j++){
                    child = new Individual(generator, normal,
                            parentType, parent.getHValue(), parent.getAlpha());
                    if(generator.nextDouble()<H_MUT_RATE)
                        child.mutateH();
                    if(generator.nextDouble()<MUT_RATE)
                        child.mutateType();
                    if(generator.nextDouble()<MUT_RATE && evolveAlpha)
                        child.mutateAlpha();
                    offspring.add(child);
                }
            }
        }
        return offspring;
        
    }
    
    private void doMigration(){
        //initialise migant pool
        for(int md=0; md<NUM_DEMES; md++){
            migrantPool[md] = new ArrayList<Individual>();
        }
        ArrayList<Individual> natalDeme;
        //ArrayList<Individual> immigrantDeme;
        Individual migrant;
        int r;
        for(int d=0; d<NUM_DEMES; d++){
            natalDeme = demes[d];
            for(int i=natalDeme.size()-1; i>=0; i--){
                if(generator.nextDouble()<MIG_RATE){
                    migrant = natalDeme.get(i);
                    //migrant.setAlpha(0.5); //MIGRANTS CANNOT BE LEADERS
                    //choose random deme to migrate to, excluding natal deme
                    r = (int)(Math.round((NUM_DEMES-1)*generator.nextDouble()));
                    while(r == d){
                        r = (int)(Math.round((NUM_DEMES-1)*
                                generator.nextDouble()));
                    }
                    migrantPool[r].add(migrant);
                    //natalDeme.remove(i);
                    demes[d].remove(i);
                }
            }
        }
        //add individuals from migrant pool into their new demes
        for(int d=0; d<NUM_DEMES; d++){
            for(int i=0; i<migrantPool[d].size(); i++){
                demes[d].add(migrantPool[d].get(i));
            }
        }
    }
    
    
    public Results runFunc(){
        
       
        initDemes(); //was init_demes
        
        int demeNumCoop=0;
        int demeNumSelfish=0;
        int demeNumLoner=0;
        int[] demeNums;
        double groupH;
        int consensusTime;
        double publicGood;
        double helpingGood;
        double punishmentGood;
        double pcPun;
        double kLoner;
        double kGroup;
        double rCoop;
        double rSelfish;
        double rLoner;
        double lambdaCoop;
        double lambdaSelfish;
        double lambdaLoner;
        double kGroupTotal;
        double groupHTotal;
        double inst_cost;
        int consensusTimeTotal = 0;
        double instCostTotal;
        
        double socialHTotal = 0;
        double[] gp;
        double [][] demeCoop = new double[T][NUM_DEMES];
        double [][] demeSelfish = new double[T][NUM_DEMES];
        double[][] demeHDiff = new double[T][NUM_DEMES];  
        ArrayList<Individual> deme;
        for(int t=0; t<T; t++){
        	socialHTotal = 0;
            if(t % 100 == 0){
            	System.out.println(t);
            }
            kGroupTotal=0;
            groupHTotal=0;
            consensusTimeTotal = 0;
            instCostTotal = 0;
            for (int d = 0; d < NUM_DEMES; d++) {
                deme = demes[d];
                socialHTotal += meanH(deme);
             
                demeNums = localTypeNums(deme);
                demeNumCoop = demeNums[C];
                demeCoop[t][d] = (double)demeNumCoop / deme.size();
                demeSelfish[t][d]= (double)demeNumSelfish / deme.size();
                
                demeNumSelfish = demeNums[S];
                demeNumLoner = demeNums[L];
                //political game
                
                Consensus consensus = setHByConsensus(deme);
                groupH = consensus.getH();
          
                demeHDiff[t][d] = groupH - meanH(deme);
                consensusTime = consensus.getTimeToConsensus();
                groupHRecord[t][d] = groupH;
                skewnessRecord[t][d] = Utility.skewnessAlpha(deme);
               
                groupHTotal += groupH;
                consensusTimeTotal += consensusTime;
                inst_cost = instCostScale * consensusTime;
                instCostTotal += inst_cost;
                publicGood = demeNumCoop * coopB;
                helpingGood = groupH * publicGood;
                helpingGood = beta*(1-Math.exp(-grad*helpingGood));
                punishmentGood = (1 - groupH) * publicGood;
                pcPun = punishmentGood / demeNumSelfish;
                //print format: gen, deme, coop, p, alpha, p_consensus, consensus_time, i_cost, demeSkew
                if(t % printEvery == 0){
                	for(Individual indiv : deme){
                		String indivType;
                		if(indiv.getType() == Model.C){
                			indivType = "C";
                		}
                		else if(indiv.getType() == Model.S){
                			indivType = "D";
                		}
                		else{
                			indivType = "A";
                		}
                		indivRecord.add(t + "," + d + "," + indivType + "," + (1-indiv.getHValue()) +
                				"," + indiv.getAlpha() + "," + (1-groupH) + "," + consensusTime + "," +
                				inst_cost + "," + skewnessRecord[t][d]);
                	}
                }
                kLoner = kBase;
                kGroup = kBase + helpingGood;
                kGroupTotal += kGroup;
                //payoff
                rLoner = rBase;
                rCoop = rBase - coopC - inst_cost;
                if (pcPun > 0) {
                    rSelfish = rBase - pcPun - inst_cost;
                } else {
                    rSelfish = rBase - inst_cost;
                }

                if (rCoop < 1) {
                    rCoop = 1;
                }

                if (rSelfish < 1) {
                    rSelfish = 1;
                }
                
                lambdaCoop = rCoop / (1 + (1 / kGroup)
                        * (demeNumCoop + demeNumSelfish) + 
                        a_gl * demeNumLoner);
                lambdaSelfish = rSelfish / (1 + (1 / kGroup)
                        * (demeNumCoop + demeNumSelfish)
                        + a_gl * demeNumLoner);
                lambdaLoner = rLoner / (1 + (1 / kLoner) * demeNumLoner
                        + a_lg * (demeNumCoop + demeNumSelfish));
                
                //non-overlapping generations.
                demes[d]=generateOffspring(deme, lambdaCoop, lambdaSelfish,
                        lambdaLoner);

            }
            
            //synchronous migration
            doMigration();
            //record stats
            kGroupMean[t] = kGroupTotal/NUM_DEMES;
            //group_meanH[t] = groupH_total/NUM_DEMES;
            meanH[t] = meanHAllIndivs();
            //meanH[t] = socialHTotal / NUM_DEMES;
            alphaMean[t] = meanAlpha();
            groupMeanConsensusTime[t] = consensusTimeTotal / NUM_DEMES;
            instCostMean[t] = instCostTotal / NUM_DEMES;
            numLeaders[t] = numLeaders();
            propLeaders[t] = (double)numLeaders[t]/popSize();
            alphaSkewness[t] = skewAlpha();
            gp=globalProps();
            globalProps[t][Model.C]=gp[Model.C];
            globalProps[t][Model.S]=gp[Model.S];
            globalProps[t][Model.L]=gp[Model.L];
            
            
        }
        return new Results
        		(globalProps, kGroupMean, meanH, this.groupMeanConsensusTime, alphaMean, 
        				instCostMean, alphaSkewness, numLeaders, propLeaders,
        				groupHRecord, skewnessRecord, demeCoop, demeSelfish,
        				demeHDiff, indivRecord);
        
        
    }
    
}
