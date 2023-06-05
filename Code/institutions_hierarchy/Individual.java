package institutions_hierarchy;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.math.*;
import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister;
/**
 *
 * @author spowers
 */
public class Individual {
    
    private int type;
    private double hValue; //inherited h (x) value;
    private double alpha;
    private double hBehaviourTime; //stores the h being altered during consensus
    private MersenneTwister generator;
    private Normal normal;
    
    public Individual(MersenneTwister generator, Normal normal){
        //default constructor
        this.generator=generator;
        this.normal=normal;
        this.type=Model.L;
        //this.hValue=Math.random();
        this.hValue=generator.nextDouble();
        this.hBehaviourTime=hValue;
        this.alpha = 0.5;
        
    }
    
    public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public double getHBehaviourTime() {
		return hBehaviourTime;
	}

	public void setHBehaviourTime(double hBehaviourTime) {
		this.hBehaviourTime = hBehaviourTime;
	}

	public Individual(MersenneTwister generator, Normal normal, int type){
        this.generator=generator;
        this.normal=normal;
        this.type=type;
        //this.hValue=Math.random();
        this.hValue=generator.nextDouble();
    }
    
    public Individual(MersenneTwister generator, Normal normal, 
            int type, double hValue){
        this.generator=generator;
        this.normal=normal;
        this.type=type;
        this.hValue=hValue;
    }
    
    public Individual(MersenneTwister generator, Normal normal, 
            int type, double hValue, double alpha){
        this.generator=generator;
        this.normal=normal;
        this.type=type;
        this.hValue=hValue;
        this.hBehaviourTime=hValue;
        this.alpha = alpha;
    }
    
   
    
   
    
    public int getType(){
        return this.type;
    }
    
    public double getHValue(){
        return this.hValue;
    }
    
    public void setType(int newType){
        //call to change type by mutation
        this.type=newType;
    }
    
    public void setHValue(double newHValue){
        this.hValue=newHValue;
    }
    
    public void mutateH(){
    	double newH;
    	newH = hValue + normal.nextDouble();
        if(newH<0){
            newH=0.00; //was 0
        }else if(newH>1){
            newH=1.0; //was 1
        }
        this.hValue=newH;
    }
    
    public void mutateType(){
    	int r = (int)Math.round((Model.NUM_TYPES-1)*generator.nextDouble());
        while(r==this.type){
            r = (int)Math.round((Model.NUM_TYPES-1)*generator.nextDouble());
        }
        this.type=r;
    }
    
    public void mutateAlpha(){
    	double newAlpha;
        newAlpha = alpha + normal.nextDouble();
        if(newAlpha < 0){
            newAlpha=0;
        }else if(newAlpha > 1){
            newAlpha=1;
        }
        this.alpha = newAlpha;
    }
    
    public void mutate(){
        double newH;
        if(generator.nextDouble() < 0.33){
            //mutate the Hvalue locus
            newH = hValue + normal.nextDouble();
            if(newH<0){
                newH=0;
            }else if(newH>1){
                newH=1;
            }
            this.hValue=newH;
        }
        else if(generator.nextDouble() < 0.66){
            //mutate the alpha locus
        	double newAlpha;
            newAlpha = alpha + normal.nextDouble();
            if(newAlpha < 0){
                newAlpha=0;
            }else if(newAlpha > 1){
                newAlpha=1;
            }
            this.alpha = newAlpha;
        }
        else{
            int r = (int)Math.round((Model.NUM_TYPES-1)*generator.nextDouble());
            while(r==this.type){
                r = (int)Math.round((Model.NUM_TYPES-1)*generator.nextDouble());
            }
            this.type=r;
        }
    }
    
}
