
/**
 * This code finds all the genes in a DNA string and then store 
 * them using the StorageResource class.
 * 
 * @author Nirlan Neckir Zamprogno de Souza
 * @version 18/03/2016
 */
import edu.duke.*;
import java.io.*;

public class StringGenes {
    public int findStopIndex(String dna, int index) {
        int stop1 = dna.indexOf("TAG", index +3);
        int stop2 = dna.indexOf("TGA", index +3);
        int stop3 = dna.indexOf("TAA", index +3);
        if (stop1 == -1 || (stop1 - index) % 3 != 0) {
            stop1 = dna.length();
        }
        if (stop2 == -1 || (stop2 - index) % 3 != 0) {
            stop2 = dna.length();
        }
        if (stop3 == -1 || (stop3 - index) % 3 != 0) {
            stop3 = dna.length();
        }
        return Math.min(stop1, Math.min(stop2, stop3));
    }
    public float cgRatio(String dna) {
        StorageResource storex = new StorageResource();
        StorageResource storey = new StorageResource();
        dna = dna.toUpperCase();
        int idx = 0;
        int idy = 0;
        while (true) {
            if (idx+1 > dna.length()) {
                break;
            }
            int c = dna.indexOf("C", idx);
            if (c == -1) {
                break;
            }
            String cs = dna.substring(c);            
            storex.add(cs);            
            idx = c +1;            
        }
        while (true) {
            if (idy+1 > dna.length()) {
                break;
            }
            int g = dna.indexOf("G", idy);
            if (g == -1) {
                break;
            }
            String gs = dna.substring(g);
            storey.add(gs);
            idy = g +1;
        } 
        int size = storex.size() + storey.size();
        //float cgRatio = size / dna.length();
        //return cgRatio;
        return (float)size/dna.length();
    }
     public StorageResource storeAll(String dna) {
        StorageResource store = new StorageResource();
        dna = dna.toUpperCase();
        int start = 0;    
        while (true) {
            if ((start+1) >= dna.length()) {
                break;
            }
            int index = dna.indexOf("ATG", start);
            if (index == -1) {
                break;                  
            }
            int stop = findStopIndex(dna, index);
            //if stop returns dna.legth() is because the no codon was found this turn
            while (stop == dna.length()) {
                if ((index+4) >= dna.length()) {
                    break;
                }
                index = dna.indexOf("ATG", index +3);
                if (index == -1) {
                    break;
                }
                stop = findStopIndex(dna, index);
            }
            if (index == -1 || (stop+4) >= dna.length() || (index+1) >= dna.length()) {
                break;                  
            }
            String found = dna.substring(index, stop +3);
            store.add(found);
            start = stop +3;
        }
        return store;
    }
    public void testStorageFinder() {
        FileResource fr = new FileResource();
        String dna = fr.asString();
        StorageResource s1 = storeAll(dna);
        System.out.println("Genes found: " + s1.size());
        System.out.println("");
        printGenes(s1);
        //codonCTG(dna);
    }
    public void printGenes(StorageResource sr) {
        StorageResource strlen1 = new StorageResource();
        StorageResource strlen2 = new StorageResource();
        
        System.out.println("Strings with more than 60 characters: ");
        System.out.println("");
        for (String str : sr.data()) {
            if (str.length() > 60) {                               
                System.out.println(str);
                strlen1.add(str);
            }
        }
        System.out.println("");
        System.out.println("Number of strings: " + strlen1.size());
        System.out.println("");
        System.out.println("");
        System.out.println("");
        System.out.println("");
        System.out.println("C-G ratio higher than 0.35: ");
        System.out.println("");
        for (String str : sr.data()) {
            if (cgRatio(str) > 0.35) {
                System.out.println(str);
                strlen2.add(str);
            }
        }
        System.out.println("");
        System.out.println("Number of strings: " + strlen2.size());
    }
    public StorageResource storeCtg(String dna) {
        StorageResource sr = new StorageResource();                      
        int ctgidx = 0;
        //int numCtg = 0;
        while (true) {
            int index = dna.indexOf("CTG", ctgidx);
            if (index == -1) {
                break;
            }
            //numCtg += 1;
            ctgidx = index +1; 
            if (index +1 > dna.length()) {
                break;
            }
            String found = dna.substring(index, index +3);
            sr.add(found);
        }
            return sr;
    }
    public void ctgCount() {
        FileResource fr = new FileResource();
        String dna = fr.asString();
        dna = dna.toUpperCase();
        StorageResource sr1 = storeCtg(dna);
        int size = sr1.size();
        System.out.println("Number of CTG in the DNA string: " + size);
    }
    public void longestGeneSize() {
        FileResource fr = new FileResource();
        String dna = fr.asString();
        StorageResource s1 = storeAll(dna);
        String longest = findLongestGene(s1);
        System.out.println("Size of the longest gene: " + longest.length() + " characters.");
    }
    public String findLongestGene(StorageResource sr) {               
        int maxLength = 0;
        String longestString = null;
        for (String s : sr.data()) {        
        if (s.length() > maxLength) {
            maxLength = s.length();
            longestString = s;      
        }        
        }
        return longestString;
    }
}
    