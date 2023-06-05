package institutions_hierarchy;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.io.IOException;

/**
 *
 * @author spowers
 */
public class Results {
    
    private double[][] globalProps;
    private double[] kGroupMean;
    private double[] hMean;
    private int[] groupMeanConsensusTime;
    private double[] alphaMean;
    private double[] instCostMean;
    private double[] alphaSkewness;
    private int[] numLeaders;
    private double[] propLeaders;
    
    private double[][] groupHRecord; //record of agreed H on each patch at each generation
    private double[][] skewnessRecord;
    private double[][] demeCoop; 
    private double[][] demeSelfish;
    
    private double[][] demeHDiff;
    private ArrayList<String> indivRecord;
    
    public Results(double[][] globalProps, double[] kGroupMean,
                    double[] groupHMean, int[] groupMeanConsensusTime,
                    double[] alphaMean, double[] instCostMean,
                    double[] alphaSkewness, int[] numLeaders, double[] propLeaders,
                    double[][] groupHRecord, double[][] skewnessRecord,
                    double[][] demeCoop, double[][] demeSelfish,
                    double[][] demeHDiff,
                    ArrayList<String> indivRecord){
        this.globalProps = globalProps;
        this.kGroupMean = kGroupMean;
        this.hMean=groupHMean;
        this.groupMeanConsensusTime = groupMeanConsensusTime;
        this.alphaMean = alphaMean;
        this.instCostMean = instCostMean;
        this.alphaSkewness = alphaSkewness;
        this.numLeaders = numLeaders;
        this.propLeaders = propLeaders;
        this.groupHRecord = groupHRecord;
        this.skewnessRecord = skewnessRecord;
        this.demeCoop = demeCoop;
        this.demeSelfish = demeSelfish;
        this.demeHDiff = demeHDiff;
        this.indivRecord = indivRecord;
    }
    
    public double meanFreqCoop(){
    	double sumFreq = 0;
    	for(int row=0; row<globalProps.length; row++){
            sumFreq += globalProps[row][Model.C];
    	}
    	
    	return sumFreq / globalProps.length;
    }
    
    public double meanFreqSelfish(){
    	double sumFreq = 0;
    	for(int row=0; row<globalProps.length; row++){
            sumFreq += globalProps[row][Model.S];
    	}
    	
    	return sumFreq / globalProps.length;
    }
    
    public double meanFreqLoner(){
    	double sumFreq = 0;
    	for(int row=0; row<globalProps.length; row++){
            sumFreq += globalProps[row][Model.L];
    	}
    	
    	return sumFreq / globalProps.length;
    }
    
    public double meanAlphaSkew(){
    	double sumAlphaSkew = 0;
    	for(double alphaSkew : alphaSkewness){
    		sumAlphaSkew += alphaSkew;
    	}
    	return sumAlphaSkew / alphaSkewness.length;
    }
    
    /**
     * @return the mean H trait of all individuals
     */
    public double meanH(){
    	double sumH = 0;
    	for(double h : hMean){
    		sumH += h;
    	}
    	return sumH / hMean.length;
    }
    
    public double meanGroupK(){
    	double kTotal = 0;
    	for(double k : this.kGroupMean){
    		kTotal += k;
    	}
    	return kTotal / kGroupMean.length;
    }
    
    public double meanHByConsensus(){
    	double demesTotalH = 0;
    	for(int t=0; t<groupHRecord.length; t++){
    		double demeTotalH = 0;
    		for(int d=0; d<groupHRecord[0].length; d++){
    			demeTotalH += groupHRecord[t][d];
    		}
    		demesTotalH += (demeTotalH / groupHRecord[0].length);
    	}
    	return demesTotalH / groupHRecord.length;
    }
    
    public double meanAlpha(){
    	double sumAlpha = 0;
    	for(double alpha : alphaMean){
    		sumAlpha += alpha;
    	}
    	return sumAlpha / alphaMean.length;
    }
    
    public double meanInstCost(){
    	double sumInstCost = 0;
    	for(double i : instCostMean){
    		sumInstCost += i;
    	}
    	return sumInstCost / instCostMean.length;
    }
    
   
    
    public double[][] getGroupHRecord() {
		return groupHRecord;
	}

	public double[][] getSkewnessRecord() {
		return skewnessRecord;
	}

	public void printGlobalProps(){
        for(int row=0; row<globalProps.length; row++){
            for(int col=0; col<globalProps[0].length-1; col++){
                //System.out.print(global_props[row][col] + ", ");
                System.out.printf("%1.3f, ", globalProps[row][col]);
            }
            //System.out.println(global_props[row][global_props[0].length-1]);
            System.out.printf("%1.3f %n",
                    globalProps[row][globalProps[0].length-1]);
        }
        
    }
    
    public void printK_group_mean(){
        for(int i=0; i<kGroupMean.length; i++){
            System.out.println(kGroupMean[i]);
        }
    }
    
    public void printGroup_hMean(){
        for(int i=0; i<hMean.length; i++){
            System.out.printf("%1.3f %n",hMean[i]);
        }
    }
    
    public void printAll(int every_t){
        //print stats every_t generations.
        System.out.println("Generation, propCoop, propSelfish, propLoner, " + 
                "meanGroupK, meanGroupH, meanConsensusTime, meanAlpha, meanInstCost,"
                + " alpha_skewness, numLeaders, propLeaders");
        //always print initial stats on first generation
        //int maxGens = global_props.length;
        //double numGenDigits = Math.ceil(maxGens/10);
        System.out.printf("%6d, ", 1);
        for(int col=0; col<globalProps[0].length; col++){
            System.out.printf("%1.4f, ", globalProps[0][col]);
        }
        System.out.printf("%4.1f, ", kGroupMean[0]);
        System.out.printf("%1.3f, ", hMean[0]);
        System.out.printf("%6d, ", this.groupMeanConsensusTime[0]);
        System.out.printf("%1.3f, ", alphaMean[0]);
        System.out.printf("%1.3f, ", instCostMean[0]);
        System.out.printf("%1.3f, ", alphaSkewness[0]);
        System.out.printf("%3d, ", numLeaders[0]);
        System.out.printf("%1.3f %n", propLeaders[0]);
        for (int row = 1; row < globalProps.length; row++) {
            if ((row % every_t) == 0) {
                //print stats this generation.
                System.out.printf("%6d, ", row+1);
                for (int col = 0; col < globalProps[0].length; col++) {
                    System.out.printf("%1.4f, ", globalProps[row][col]);
                }
                System.out.printf("%4.1f, ", kGroupMean[row]);
                System.out.printf("%1.3f, ", hMean[row]);
                System.out.printf("%6d, ", this.groupMeanConsensusTime[row]);
                System.out.printf("%1.3f, ", alphaMean[row]);
                System.out.printf("%1.3f, ", instCostMean[row]);
                System.out.printf("%1.3f, ", alphaSkewness[row]);
                System.out.printf("%3d, ", numLeaders[row]);
                System.out.printf("%1.3f %n", propLeaders[row]);
            }
        }
        //System.out.println()
        
    }
    
    public void writeToFile(int every_t, String file_path) throws IOException{
        FileWriter write = new FileWriter(file_path);
        PrintWriter printer = new PrintWriter(write);
        FileWriter writeP = new FileWriter("p.csv"); //P is the new H
        PrintWriter printP = new PrintWriter(writeP);
        FileWriter writeSkew = new FileWriter("skew.csv");
        PrintWriter printSkew = new PrintWriter(writeSkew);
        FileWriter writeCoop = new FileWriter("coop.csv");
        PrintWriter printCoop = new PrintWriter(writeCoop);
        FileWriter writeSelfish = new FileWriter("selfish.csv");
        PrintWriter printSelfish = new PrintWriter(writeSelfish);
        FileWriter writePDiff = new FileWriter("Pdiff.csv");
        PrintWriter printPDiff = new PrintWriter(writePDiff);
        //print stats every_t generations.
        printer.println("Generation, propCoop, propSelfish, propLoner, " + 
                "meanGroupK, meanGroupP");
        //always print initial stats on first generation
        //int maxGens = global_props.length;
        //double numGenDigits = Math.ceil(maxGens/10);
        printer.printf("%6d, ", 1);
        for(int col=0; col<globalProps[0].length; col++){
            printer.printf("%1.4f, ", globalProps[0][col]);
        }
        printer.printf("%4.1f, ", kGroupMean[0]);
        printer.printf("%1.4f %n", 1 - hMean[0]);
        for (int row = 1; row < globalProps.length; row++) {
            if ((row % every_t) == 0) {
                //print stats this generation.
                printer.printf("%6d, ", row);
                for (int col = 0; col < globalProps[0].length; col++) {
                    printer.printf("%1.4f, ", globalProps[row][col]);
                }
                printer.printf("%4.1f, ", kGroupMean[row]);
                printer.printf("%1.4f %n", 1 - hMean[row]);
            }
        }
        printer.close();
        
        for(int t=0; t<groupHRecord.length; t++){
        	for(int d=0; d<groupHRecord[0].length; d++){
        		printP.printf("%1.4f,", 1 - groupHRecord[t][d]);
        		printSkew.printf("%1.4f,", skewnessRecord[t][d]);
        		printCoop.printf("%1.4f,", demeCoop[t][d]);
        		printSelfish.printf("%1.4f,", demeSelfish[t][d]);
        		printPDiff.printf("%1.4f,", demeHDiff[t][d]);
        	}
        	printP.println();
        	printSkew.println();
        	printCoop.println();
        	printSelfish.println();
        	printPDiff.println();
        }
        printP.close();
        printSkew.close();
        printCoop.close();
        printSelfish.close();
        printPDiff.close();
        
    }
    
    public void printIndivs(String filename) throws IOException{
    	 FileWriter write = new FileWriter(filename);
         PrintWriter printer = new PrintWriter(write);
         //print format: gen, deme, coop, p, alpha, p_consensus, consensus_time, i_cost, demeSkew
         printer.println("generation,deme,socialType,p,alpha,p_consensus,consensus_time,instCost,demeSkew");
         for(String s : indivRecord){
        	 printer.println(s); 
         }
         printer.close();
    }
    
}
