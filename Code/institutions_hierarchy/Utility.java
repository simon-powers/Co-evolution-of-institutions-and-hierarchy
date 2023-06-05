package institutions_hierarchy;


import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;
public class Utility {
   private static Skewness skewnessCalculator = new Skewness();
    // Probability test
    public static int testProb(double pMu){
        if(Math.random() > pMu){
            return 0;
        }
        else{
           return 1;
        }
    }
    //Writer
    public static List writeFile(List pList){
        
        return(pList);
    }
    
  //Sample an individual except one
    public static int randomSampleOther(int pArrayL, int pIndex){
        //if(pLArray == 1){System.out.println("BUG randomSampleOther List too short");}
        ArrayList<Integer> resList = new ArrayList<>();
        for(int i=0; i<pArrayL; i++){resList.add(i);}
        resList.remove(pIndex);
        int index = (int)(Math.floor(Math.random()*(pArrayL-1)));
        return resList.get(index);
        }
    
    //Sample a number of individuals except one
    public static int[] randomSampleOtherList(int pArrayL, int pSampleSize, int pIndex){
        //if(pLArray == 1){System.out.println("BUG randomSampleOther List too short");}
        ArrayList<Integer> resList = new ArrayList<>();
        for(int i=0; i<pArrayL; i++){resList.add(i);}
        resList.remove(pIndex);
        int indexTemp;
        int[] indexL = new int[pSampleSize];
        for(int i=0; i<pSampleSize; i++){
            indexTemp = (int)(Math.floor(Math.random()*(resList.size())));
            indexL[i]= resList.get(indexTemp);
            resList.remove(indexTemp);
        }
        return indexL;
        }

    public static int probSample(double[] pArrayProb, double pKey){             //Binary search algorithm for sampling probability
                            int lower = 0;
                            int upper = pArrayProb.length-1;
                            int mid;
                            while (lower < upper){ 
                                mid = (int)Math.floor((lower + upper )/2);      
                                if((pArrayProb[mid] - pKey) > 0){
                                    upper = mid;
                                }
                                else{
                                    lower = mid + 1;
                                }
                            }
                            return lower;
                        }
    
    //Calcul of standart deviation
    public static double sdPref(List<Individual> pList){
        double M= 0;
        double S = 0;
        double oldM;
        for(int i = 0; i < pList.size(); i++){
            oldM = M; 
            M += (pList.get(i).getHBehaviourTime()-M)/(i+1);
            S += (pList.get(i).getHBehaviourTime()-M) * (pList.get(i).getHBehaviourTime()-oldM);
        }
        return (Math.sqrt(S/(pList.size()-1)));
    }
    //Calcul of mean of the preferences
    public static double meanPref(List<Individual> pList){
        double sumF = 0;
        for(int i = 0; i < pList.size(); i++){
            sumF += pList.get(i).getHBehaviourTime();
        }
        return(sumF/(pList.size()));
    }
    
    public static double skewnessAlpha(List<Individual> pList) {
    	double[] alphaList = new double[pList.size()];
    	for(int i =0; i<pList.size(); i++) {alphaList[i]=pList.get(i).getAlpha();}
    	double res = skewnessCalculator.evaluate(alphaList);
    	return res;
    }
}
