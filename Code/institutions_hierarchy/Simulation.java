package institutions_hierarchy;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.util.ArrayList;
import java.io.IOException;

/**
 *
 * @author spowers
 */
public class Simulation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
    	int NUM_DEMES=50;
    	int T=10000;
    	int k_base=20;
    	double inst_cost=0.01; //was 0.01
    	double MUT_RATE=0.01; // was 0.0001
    	double H_MUT_RATE = 0.02;
    	double MIG_RATE=0.1;
    	double beta=300;  //WAS 300
    	double grad=0.0075; 
    	double coopB=0.9;
    	double coopC=0.1;
    	int numListeners = 30; //was 30
    	double kScale = 4;
        double xThr = 0.03; //was 0.05, didn't work if too small at 0.01,
        					//as migrants have too small an effect on 
        					//consensus
        boolean evolveAlpha = true;
        //good params: H_MUT_RATE 0.02, xThr 0.025, inst cost scalar 0.01
        //best is with xThr 0.035, get 0.95 coop and 0.985 skew
        //with xThr 0.045 get 0.95 coop, 0.02 selfish, and 0.78 skew
        // with xthr 0.55 Mean coop: 0.9717026501853732, Mean selfish: 0.016909343597032892, Mean Loner: 0.011388006217571795, Mean P0.22531944100919832, Mean skew: 0.6131010217351806,
        //with xThr 0.065 Mean coop: 0.9404093868899905, Mean selfish: 0.017682015954004836, Mean Loner: 0.04190859715598152, Mean P0.22649855373142413, Mean skew: 0.3124628884082521, Mean Alpha0.45484417676915373,  Mean inst cost: 0.15921028800001796
        //with xThr 0.015 Mean coop: 0.4740309594499337, Mean selfish: 0.49194991459900156, Mean Loner: 0.03401912595107591, Mean P0.27358802397594206, Mean skew: 1.0763970325940142, Mean Alpha0.3049508179517963,  Mean inst cost: 0.14897940399999562
        //with xThr 1, Mean coop: 0.9808363801364096, Mean selfish: 0.014043081156724226, Mean Loner: 0.005120538706869239, Mean P0.21051973701431836, Mean skew: -0.028106489226438408, Mean Alpha0.5089108576984863,  Mean inst cost: 0.0
        
        //with Beta raised to 600 get group sizes of 430, Mean coop: 0.9447099550233594, Mean selfish: 0.039542882284047426, Mean Loner: 0.015747162692578413, Mean P0.32183209831158366, Mean skew: 1.0305884489064687, Mean Alpha0.3186812767530237,  Mean inst cost: 0.41794731199999546
    	/*Model(int NUM_DEMES,int T, int k_base, double inst_cost,
                        double MUT_RATE, double MIG_RATE,
                        double beta, double grad, double coopB, double coopC){*/
    	if(NUM_DEMES==1)
    		MIG_RATE = 0;
    	int every_t=10; //print stats every t generations
    	Model model = new Model(NUM_DEMES, T, k_base,
    			inst_cost, MUT_RATE,
    			MIG_RATE, beta, grad, coopB, coopC,
    			numListeners, kScale, xThr, H_MUT_RATE, every_t, evolveAlpha);
    	Results results = model.runFunc();

        
        //results.printAll(every_t);
        System.out.println("Mean coop: " + results.meanFreqCoop() + ", " 
        + "Mean selfish: " + results.meanFreqSelfish() + ", " +
        		"Mean Loner: " + results.meanFreqLoner() + ", " + 
        "Mean P" + (1- results.meanH()) + ", " + "Mean skew: " + results.meanAlphaSkew() + 
       ", " + "Mean Alpha" + results.meanAlpha() + ", " + " Mean inst cost: " + results.meanInstCost());
        String filename = "runStats,T" + T + ",k_base" + k_base + ",mig_rate"
                + MIG_RATE + ",inst_cost" + inst_cost + "mut_rate" + 
                MUT_RATE +  "beta" + beta + ".csv";
        
        /*String indivFilename = "NumGens" + T + "xThr" + xThr +"inst_cost" + inst_cost +
        		"numListeners" + numListeners + "PMut" + H_MUT_RATE +
        		"MigRate" + MIG_RATE + "beta" + beta + ".csv";*/
        String indivFilename = "Params" + "xThr" + xThr + "Beta" + beta + ".csv";
        
       
        try{
            results.writeToFile(every_t, filename);
            results.printAll(every_t);
            results.printIndivs(indivFilename);
            
        }
        catch(IOException e){
            System.out.println("Error writing to file!");
        }
    }
    
    
            
}
