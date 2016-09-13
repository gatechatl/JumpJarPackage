package METABOLOMICS.AIM;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

import org.openscience.cdk.formula.IsotopeContainer;
import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.formula.IsotopePatternManipulator;
//import org.openscience.cdk.formula.IsotopePatternSimilarity;


import CDKFunction.IsotopePatternSimilarity;
import Interface.JumpInterface;

/**
 * Generates the isotope table
 * @author tshaw
 */
public class AIMPuritySampler implements JumpInterface {

	public static double A = -1;
	public static double A_freq = -1;
	public static double B = -1;
	public static double B_freq = -1;
	public static double C = -1;
	public static double C_freq = -1;
	
	public static int C_num = 0;
	public static int H_num = 0;
	public static int N_num = 0;
	public static int S_num = 0;
	public static int O_num = 0;
	public static int P_num = 0;
	
	public static double _ppm = 10;
	public static BIGNUM cutoff = new BIGNUM(1, -50);
	
	public static void main(String[] args) {
		String formula = "C15H24N3O4";
		String charge = "1";
		String ppm = "50";
		String C12_File = "C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Drew_MISSILE_Exp1_DataFitting\\Supplementary\\" + formula + "\\" + formula + "_C12.txt";
		String C13_File = "C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Drew_MISSILE_Exp1_DataFitting\\Supplementary\\" + formula + "\\" + formula + "_C13.txt";
		String N15_File = "C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Drew_MISSILE_Exp1_DataFitting\\Supplementary\\" + formula + "\\" + formula + "_N15.txt";
		String sampled_element = "C12";
		String[] parameter = {"C15H24N3O4", C12_File, C13_File, N15_File, sampled_element, charge, ppm};
		
		AIMPuritySampler sampler = new AIMPuritySampler();
		sampler.execute(parameter);
	}
	

	/**
	 * Generating the isotope table
	 * @param args
	 */
	public void execute(String[] args) {
		
		try {
			
			String formula = args[0];
			String C12_File = args[1];
			String C13_File = args[2];
			String N15_File = args[3];
			String sampled_element = args[4];
			int charge = new Integer(args[5]);
			double ppm = new Double(args[6]);
			boolean single_check = false;
			double ppm_val = ppm;
			//for (double ppm_val = 10; ppm_val <= 100; ppm_val = ppm_val + 10) {
			_ppm = ppm_val;
			
			
			String H2_OutFile = "C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Drew_MISSILE_Exp1_DataFitting\\Supplementary\\" + formula + "\\" + formula + "_Sampling_H2.txt";
			String C12_OutFile = "C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Drew_MISSILE_Exp1_DataFitting\\Supplementary\\" + formula + "\\" + formula + "_Sampling_C12.txt";
			String C13_OutFile = "C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Drew_MISSILE_Exp1_DataFitting\\Supplementary\\" + formula + "\\" + formula + "_Sampling_C13.txt";
			String N15_OutFile = "C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Drew_MISSILE_Exp1_DataFitting\\Supplementary\\" + formula + "\\" + formula + "_Sampling_N15.txt";
			
			double big_step = 0.05;
			
			int _charge = charge;
			
			//String formula = "C100N29O36H147"; // this looks the best

			//String formula = "C100N29O28H274";
			//String formula = "C100N29O26S1H274";
			//String formula = "C100N29O24S2H274"; // pretty close too
			//String formula = "C100N29O34SH147"; //
			//String formula = "O9H290C100S9N29"; // another candidate
			//String formula = "C100N29O33S5H36"; //
			//String formula = "C100N29O21S7H163";
			
			//

			//printTime();
			//IsotopePattern reference_pattern = getPeakInfo("C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\DataFitting\\C24H49NO7P_C12.peak");
			//IsotopePattern reference_pattern = getPeakInfo("C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\DataFitting\\unknown_C12.peak");
			IsotopePattern reference_pattern = getPeakInfo(C12_File);
			
			reference_pattern = IsotopePatternManipulator.normalize(reference_pattern);
			
			double C12 = 0.9893;
			double N14 = 0.99632;
			double H1 = 0.999885;
			
			double simulated_C12 = -1;
			double simulated_N15 = -1;
			double simulated_C13 = -1;
			double simulated_H2 = -1;
			IsotopePattern query_pattern = null;
			
			if (single_check) {
				query_pattern = calculate_pattern(formula, C12, N14, H1, _charge);
				
				for (int i = 0; i < query_pattern.getNumberOfIsotopes(); i++) {
					IsotopeContainer container = query_pattern.getIsotope(i);
					System.out.println(container.getMass() + "\t" + container.getIntensity());
				}
				System.exit(0);
			}
			
			boolean check_C12 = false;
			boolean check_C13 = false;
			boolean check_H2 = false;
			boolean check_N15 = false;
			if (sampled_element.equals("C12")) {
				check_C12 = true;
			}
			if (sampled_element.equals("C13")) {
				check_C13 = true;
			}
			if (sampled_element.equals("N15")) {
				check_N15 = true;
			}
			double step = 0.001;
			
			
			
			if (check_H2) {
					
				boolean doOnce = true;
				double smallest = 0;
				double largest = 1;
				double prevIndex = 0;
				for (double H1_sample = 0; H1_sample <= 1; H1_sample += big_step) {
					
					query_pattern = calculate_pattern(formula, C12, N14, H1_sample, _charge);
					
					//System.out.println("Number of isotopes: " + query_pattern.getNumberOfIsotopes());
					IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
					double similarity_score = similarity.compare(reference_pattern, query_pattern);
					if (similarity_score > 0) {
						if (doOnce) {
							smallest = prevIndex;
							doOnce = false;
						}
						largest = H1_sample;
					}
					prevIndex = H1_sample;
				}
				largest = largest + big_step;
				if (largest <= (1 - big_step)) {
					largest += big_step;
					if (largest >= 1) {
						largest = 1.0;
					}
				}
				
				double max_score = 0;
				for (double H1_sample = 0; H1_sample < smallest; H1_sample += step) {
					System.out.println(H1_sample + "\t" + 0.0);
				}
				for (double H1_sample = smallest; H1_sample <= largest; H1_sample += step) {
				
					query_pattern = calculate_pattern(formula, C12, N14, H1_sample, _charge);
					
					//System.out.println("Number of isotopes: " + query_pattern.getNumberOfIsotopes());
					IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
					double similarity_score = similarity.compare(reference_pattern, query_pattern);
					if (max_score < similarity_score) {
						max_score = similarity_score;
						simulated_H2 = H1_sample;
					}
					//System.out.println(similarity_score);
					System.out.println(H1_sample + "\t" + similarity_score);
					
				}
				for (double H1_sample = largest + step; H1_sample <= 1.0; H1_sample += step) {
					System.out.println(H1_sample + "\t" + 0.0);
				}
			}
			
			if (check_C13) {
				reference_pattern = getPeakInfo(C13_File);
				reference_pattern = IsotopePatternManipulator.normalize(reference_pattern);
				
				boolean doOnce = true;
				double smallest = 0;
				double largest = 1;
				double prevIndex = 0;
				for (double C12_sample = 0; C12_sample <= 1; C12_sample += big_step) {
					
					query_pattern = calculate_pattern(formula, C12_sample, N14, H1, _charge);					
					//System.out.println("Number of isotopes: " + query_pattern.getNumberOfIsotopes());
					IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
					double similarity_score = similarity.compare(reference_pattern, query_pattern);
					if (similarity_score > 0) {
						if (doOnce) {
							smallest = prevIndex;
							doOnce = false;
						}
						largest = C12_sample;
					}
					prevIndex = C12_sample;
				}
				largest = largest + big_step;
				if (largest <= 1 - big_step) {
					largest += big_step;
					if (largest >= 1) {
						largest = 1.0;
					}
				}
				
				double max_score = 0;
				for (double C12_sample = 0; C12_sample < smallest; C12_sample += step) {
					System.out.println(C12_sample + "\t" + 0.0);
				}
				for (double C12_sample = smallest; C12_sample <= largest; C12_sample += step) {
				
					query_pattern = calculate_pattern(formula, C12_sample, N14, H1, _charge);
					
					IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
					double similarity_score = similarity.compare(reference_pattern, query_pattern);
					if (max_score < similarity_score) {
						max_score = similarity_score;
						simulated_C13 = C12_sample;
					}
					
					System.out.println((1 - C12_sample) + "\t" + similarity_score);
				}
				for (double C12_sample = largest + step; C12_sample <= 1.0; C12_sample += step) {
					System.out.println(C12_sample + "\t" + 0.0);
				}
				
			}
			
			if (check_N15) {
				reference_pattern = getPeakInfo(N15_File);
				reference_pattern = IsotopePatternManipulator.normalize(reference_pattern);
				
				
				boolean doOnce = true;
				double smallest = 0;
				double largest = 1;
				double prevIndex = 0;
				for (double N14_sample = 0; N14_sample < smallest; N14_sample += step) {
					System.out.println(N14_sample + "\t" + 0.0);
				}
				for (double N14_sample = 0; N14_sample <= 1; N14_sample += big_step) {
					
					query_pattern = calculate_pattern(formula, C12, N14_sample, H1, _charge);
					
					//System.out.println("Number of isotopes: " + query_pattern.getNumberOfIsotopes());
					IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
					double similarity_score = similarity.compare(reference_pattern, query_pattern);
					if (similarity_score > 0) {
						if (doOnce) {
							smallest = prevIndex;
							doOnce = false;
						}
						largest = N14_sample;
					}
					prevIndex = N14_sample;
				}
				largest = largest + big_step;
				if (largest <= 1 - big_step) {
					largest += big_step;
					if (largest >= 1) {
						largest = 1.0;
					}
				}
				
				double max_score = 0;
				for (double N14_sample = smallest; N14_sample <= largest; N14_sample += step) {
				//double N14_sample = 0.01;
					query_pattern = calculate_pattern(formula, C12, N14_sample, H1, _charge);
					
					//System.out.println("Number of isotopes: " + query_pattern.getNumberOfIsotopes());
					IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
					double similarity_score = similarity.compare(reference_pattern, query_pattern);
					if (max_score < similarity_score) {
						max_score = similarity_score;
						simulated_N15 = N14_sample;
					}
					//System.out.println(similarity_score);
					
					System.out.println(N14_sample + "\t" + similarity_score);
				}
				for (double N14_sample = largest + step; N14_sample <= 1.0; N14_sample += step) {
					System.out.println(N14_sample + "\t" + 0.0);
				}
			}

			if (check_C12) {
				reference_pattern = getPeakInfo(C12_File);
				reference_pattern = IsotopePatternManipulator.normalize(reference_pattern);
				
				
				boolean doOnce = true;
				double smallest = 0;
				double largest = 1;
				double prevIndex = 0;
				for (double C12_sample = 0; C12_sample <= 1; C12_sample += big_step) {
					
					query_pattern = calculate_pattern(formula, C12_sample, N14, H1, _charge);
					
					//System.out.println("Number of isotopes: " + query_pattern.getNumberOfIsotopes());
					IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
					double similarity_score = similarity.compare(reference_pattern, query_pattern);
					if (similarity_score > 0) {
						if (doOnce) {
							smallest = prevIndex;
							doOnce = false;
						}
						largest = C12_sample;
					}
					prevIndex = C12_sample;
				}
				largest = largest + big_step;
				if (largest <= 1 - big_step) {
					largest += big_step;
					if (largest >= 1) {
						largest = 1.0;
					}
				}
				
				
				double max_score = 0;
				for (double C12_sample = 0; C12_sample < smallest; C12_sample += step) {
					System.out.println(C12_sample + "\t" + 0.0);
				}
				for (double C12_sample = smallest; C12_sample <= largest; C12_sample += step) {
				
					query_pattern = calculate_pattern(formula, C12_sample, N14, H1, _charge);
					
					//System.out.println("Number of isotopes: " + query_pattern.getNumberOfIsotopes());
					IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
					double similarity_score = similarity.compare(reference_pattern, query_pattern);
					if (max_score < similarity_score) {
						max_score = similarity_score;
						simulated_C12 = C12_sample;
					}
					//System.out.println(similarity_score);
					//out.write(C12_sample + "\t" + similarity_score + "\n");
					System.out.println(C12_sample + "\t" + similarity_score);
				}
				for (double C12_sample = largest + step; C12_sample <= 1.0; C12_sample += step) {
					System.out.println(C12_sample + "\t" + 0.0);
				}
			}
			if (sampled_element.equals("C12")) {
				System.out.println("Simulated Result: " + formula + "\tcharge:" + _charge + "\tC12:" + simulated_C12);
			} else if (sampled_element.equals("C13")) {
				double result = (1 - simulated_C13);
				if (result < 0) {
					result = 0;
				}
				System.out.println("Simulated Result: " + formula + "\tcharge:" + _charge + "\tC13:" + result);
			} else if (sampled_element.equals("N15")) {
				double result = (1 - simulated_N15);
				if (result < 0) {
					result = 0;
				}
				System.out.println("Simulated Result: " + formula + "\tcharge:" + _charge + "\tN15:" + result);
			}
		
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculate the pattern for a formula
	 * @param C12
	 * @param N14
	 * @param H1
	 * @return
	 */
	public static IsotopePattern calculate_pattern(String formula, double C12, double N14, double H1, int charge) {
		
		C_num = retrieve_num_element(formula, "C");
		H_num = retrieve_num_element(formula, "H");
		N_num = retrieve_num_element(formula, "N");
		O_num = retrieve_num_element(formula, "O");
		S_num = retrieve_num_element(formula, "S");
		P_num = retrieve_num_element(formula, "P");
		
		
		LinkedList ini_list = create_ini_list(C12, N14, H1);
		
		LinkedList listC = preload_ini_file(ini_list, "C");
		LinkedList listH = preload_ini_file(ini_list, "H");
		LinkedList listN = preload_ini_file(ini_list, "N");
		LinkedList listS = preload_ini_file(ini_list, "S");
		LinkedList listO = preload_ini_file(ini_list, "O");
		LinkedList listP = preload_ini_file(ini_list, "P");
		// possibly the following line might be causing trouble
		
		boolean tolerance_type = true;
		double ppm_val = _ppm;
		
		listC = merge_mass_tolerance(listC, tolerance_type, ppm_val);
		listH = merge_mass_tolerance(listH, tolerance_type, ppm_val);
		listN = merge_mass_tolerance(listN, tolerance_type, ppm_val);
		listS = merge_mass_tolerance(listS, tolerance_type, ppm_val);
		listO = merge_mass_tolerance(listO, tolerance_type, ppm_val);
		listP = merge_mass_tolerance(listP, tolerance_type, ppm_val);
		
		boolean ppm = true;
		
		
		
		//System.out.println(C_num + "\t" + H_num + "\t" + N_num + "\t" + O_num + "\t" + S_num + "\t" + P_num);
		//write_example(list, path, 1);
		
		
		LinkedList newList = new LinkedList();
		for (int i = 0; i < C_num; i++) {
			if (newList.size() == 0) {
				newList = listC;
			} else {
				
				newList = MergeList(newList, listC);					
				newList = merge_mass_tolerance(newList, ppm, ppm_val);
			}
		}
		for (int i = 0; i < H_num; i++) {
			if (newList.size() == 0) {
				newList = listH;
			} else {
				//LinkedList list = new LinkedList();
				//list = MergeList(list, listH);
				//print_list(list);
				newList = MergeList(newList, listH);
				newList = merge_mass_tolerance(newList, ppm, ppm_val);
			}
		}
		
		for (int i = 0; i < N_num; i++) {
			if (newList.size() == 0) {
				newList = listN;
			} else {
				newList = MergeList(newList, listN);
				newList = merge_mass_tolerance(newList, ppm, ppm_val);
			}
		}
		
		for (int i = 0; i < S_num; i++) {
			if (newList.size() == 0) {
				newList = listS;
			} else {
				newList = MergeList(newList, listS);
				newList = merge_mass_tolerance(newList, ppm, ppm_val);
			}
		}
		
		for (int i = 0; i < O_num; i++) {
			if (newList.size() == 0) {
				newList = listO;
			} else {
				newList = MergeList(newList, listO);
				newList = merge_mass_tolerance(newList, ppm, ppm_val);
			}
		}
		
		for (int i = 0; i < P_num; i++) {
			if (newList.size() == 0) {
				newList = listP;
			} else {
				newList = MergeList(newList, listP);
				newList = merge_mass_tolerance(newList, ppm, ppm_val);
			}
		}
		
		//print_list(newList);
		IsotopePattern query_pattern = convert_Compounds2IsotopePattern(newList, -1, charge);
		return query_pattern;
	}
	/**
	 * Convert a linkedlist of compound and output Isotopic Pattern
	 * @param list
	 */
	public static IsotopePattern convert_Compounds2IsotopePattern(LinkedList list, int top, int charge) {

		
		if (top <= 0) {
			top = list.size();
		}
		int count = 0;
		
	
		double[] intensities = new double[list.size()];
		Iterator itr = list.iterator();
		while (itr.hasNext()) {
			COMPOUND compound = (COMPOUND)itr.next();
				
			double mass = compound.MASS.toDouble();
			double intensity = compound.TOTALPROB.toDouble();
			intensities[count] = intensity;
			//System.out.println(mass + "\t" + intensity);
			count++;
		}
		/*
		Arrays.sort(intensities);
		
		for (double intensity: intensities) {
			System.out.println(intensity);
		}
		*/
		count = 0;
		IsotopePattern pattern = new IsotopePattern();
		itr = list.iterator();
		while (itr.hasNext()) {
			COMPOUND compound = (COMPOUND)itr.next();
			if (count < top) {
				
				double mass = compound.MASS.toDouble();
				double intensity = compound.TOTALPROB.toDouble();
				if (charge > 1) {
					//mass = (mass / charge + (charge - 1) * 1.007825);
					// mass here is neutral monoisotopic mass 
					//mass = (mass + (charge) * 1.007276466812) / charge ;
					mass = mass / charge;
				
					//mass = (mass - 1.007825) / charge + 1.007825;
				}
				pattern.addIsotope(new IsotopeContainer(mass, intensity));
				//System.out.println(mass + "\t" + intensity);
			}
			count++;
		}
		return pattern;
	}
	/**
	 * based on the input file, obtain the estimated mass
	 * @param fileName
	 * @return
	 */
	public static IsotopePattern getPeakInfo(String fileName) {
		
		try {

			IsotopePattern result = new IsotopePattern();
			
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				if (!str.trim().equals("")) {
					String[] split = str.split("\t");
					double mass = new Double(split[0]);
					double intensity = new Double(split[1]);
					result.addIsotope(new IsotopeContainer(mass, intensity));										
					//System.out.println`(mass + "\t" + intensity);
				}
			}
			in.close();
			
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	public static LinkedList create_ini_list(double C12, double N14, double H1) {
		LinkedList list = new LinkedList();

		list.add("C: 12," + C12);		
		list.add("C: 13.0033548378," + (1 - C12));
		list.add("H: 1.0078250321," + H1);
		list.add("H: 2.0141017780," + (1 - H1));
		list.add("N: 14.0030740052," + N14);
		list.add("N: 15.0001088984," + (1 - N14));
		list.add("O: 15.9949146221,0.99757");
		list.add("O: 16.99913150,0.00038");
		list.add("O: 17.9991604,0.00205");
		list.add("S: 31.97207069,0.9493");
		list.add("S: 32.97145850,0.0076");
		list.add("S: 33.96786683,0.0429");
		list.add("P: 30.97376151,1.0");
		return list;
	}
	public static void write_example(LinkedList list, String path, int num) {
		
		File f = new File(path);
		if (!f.exists()) {
			f.mkdir();
		}
		
		int minnum = 0;
		for (int i = 0; i <= 100000; i = i + 1000) {
			if (num > i) {
				minnum = i;
			}
		}
		
		f = new File(path + "/" + minnum);
		if (!f.exists()) {
			f.mkdir();
		}
		try {
			FileWriter fwriter = new FileWriter(path + "/" + minnum + "/" + num + ".txt");
			BufferedWriter out = new BufferedWriter(fwriter);
			
			Iterator itr = list.iterator();
			
			double total = 0;
			
			while (itr.hasNext()) {
				
				COMPOUND compound = (COMPOUND)itr.next();
				BIGNUM prob = (BIGNUM)compound.PROB; 
				BIGNUM freq = (BIGNUM)compound.FREQ;
				BIGNUM mass = (BIGNUM)compound.MASS;
				BIGNUM total_prob = (BIGNUM)compound.TOTALPROB;				
				//int isotope = (Integer)compound.NUM_ISOTOPE;
		
				
				out.write(mass.getString() + "\t" + total_prob.getString() + "\n");
				//System.out.println(n + "\t" + weight + "\t" + val.NUM + "E" + new Double(val.E).intValue() + "\t" + count.NUM + "E" + new Double(count.E).intValue() + "\t" + mult.NUM + "E" + new Double(mult.E).intValue());
				//System.out.println(num + "\t" + val + "\t" + count + "\t" + val * count);
				//total += val * count;
			
			}
			//System.out.println(total);
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	
	public static LinkedList make_longer(LinkedList one, LinkedList list) {
		LinkedList newList = MergeList(list, one);
		return newList;
		//LinkedList shortList = condense_filer(newList);
		//return shortList;
	}
	
	/**
	 * 
	 * @param fileName
	 * @param chem
	 * @return
	 */
	public static LinkedList preload_ini_file(LinkedList list, String chem) {
		LinkedList isotope_list = new LinkedList();
		try {
			
			int C_ID = 0;
			int H_ID = 0;
			int N_ID = 0;
			int S_ID = 0;
			int O_ID = 0;
			int P_ID = 0;
			Iterator itr = list.iterator();
			while (itr.hasNext()) {				
				String str = (String)itr.next();
				if (str.startsWith("C:") && str.split(":")[0].equals(chem)) {
					str = str.replaceAll("C:",  "").trim();
					String[] split = str.split(",");
					double mass = new Double(split[0]);
					double prob = new Double(split[1]);					
					COMPOUND compound = new COMPOUND();
					compound.TYPE = "C";
					compound.SETTYPE(C_ID);
					//compound.TYPE2 = C_ID;
					compound.MASS = new BIGNUM(mass);
					compound.SIZE = 1;
					compound.PROB = new BIGNUM(prob);
					compound.FREQ = new BIGNUM(1.0);
					compound.TOTALPROB = new BIGNUM(prob);
					C_ID++;
					isotope_list.add(compound);
				}
				if (str.startsWith("H:") && str.split(":")[0].equals(chem)) {
					str = str.replaceAll("H:",  "").trim();
					String[] split = str.split(",");
					double mass = new Double(split[0]);
					double prob = new Double(split[1]);					
					COMPOUND compound = new COMPOUND();
					compound.TYPE = "H";
					compound.SETTYPE(H_ID);
					//compound.TYPE2 = H_ID;
					compound.MASS = new BIGNUM(mass);
					compound.SIZE = 1;
					compound.PROB = new BIGNUM(prob);
					compound.FREQ = new BIGNUM(1.0);
					compound.TOTALPROB = new BIGNUM(prob);
					H_ID++;
					isotope_list.add(compound);
				}
				if (str.startsWith("N:") && str.split(":")[0].equals(chem)) {
					str = str.replaceAll("N:",  "").trim();
					String[] split = str.split(",");
					double mass = new Double(split[0]);
					double prob = new Double(split[1]);					
					COMPOUND compound = new COMPOUND();
					compound.TYPE = "N";
					compound.SETTYPE(N_ID);
					//compound.TYPE2 = N_ID;
					compound.MASS = new BIGNUM(mass);
					compound.SIZE = 1;
					compound.PROB = new BIGNUM(prob);
					compound.FREQ = new BIGNUM(1.0);
					compound.TOTALPROB = new BIGNUM(prob);
					N_ID++;
					isotope_list.add(compound);
				}
				if (str.startsWith("O:") && str.split(":")[0].equals(chem)) {
					str = str.replaceAll("O:",  "").trim();
					String[] split = str.split(",");
					double mass = new Double(split[0]);
					double prob = new Double(split[1]);					
					COMPOUND compound = new COMPOUND();
					compound.TYPE = "O";
					compound.SETTYPE(O_ID);
					//compound.TYPE2 = O_ID;
					compound.MASS = new BIGNUM(mass);
					compound.SIZE = 1;
					compound.PROB = new BIGNUM(prob);
					compound.FREQ = new BIGNUM(1.0);
					compound.TOTALPROB = new BIGNUM(prob);
					O_ID++;
					isotope_list.add(compound);
				}
				if (str.startsWith("S:") && str.split(":")[0].equals(chem)) {
					str = str.replaceAll("S:",  "").trim();
					String[] split = str.split(",");
					double mass = new Double(split[0]);
					double prob = new Double(split[1]);					
					COMPOUND compound = new COMPOUND();
					compound.TYPE = "S";
					compound.SETTYPE(S_ID);
					//compound.TYPE2 = S_ID;
					compound.MASS = new BIGNUM(mass);
					compound.SIZE = 1;
					compound.PROB = new BIGNUM(prob);
					compound.FREQ = new BIGNUM(1.0);
					compound.TOTALPROB = new BIGNUM(prob);
					S_ID++;
					isotope_list.add(compound);
				}
				if (str.startsWith("P:") && str.split(":")[0].equals(chem)) {
					str = str.replaceAll("P:",  "").trim();
					String[] split = str.split(",");
					double mass = new Double(split[0]);
					double prob = new Double(split[1]);					
					COMPOUND compound = new COMPOUND();
					compound.TYPE = "P";
					compound.SETTYPE(S_ID);
					//compound.TYPE2 = S_ID;
					compound.MASS = new BIGNUM(mass);
					compound.SIZE = 1;
					compound.PROB = new BIGNUM(prob);
					compound.FREQ = new BIGNUM(1.0);
					compound.TOTALPROB = new BIGNUM(prob);
					P_ID++;
					isotope_list.add(compound);
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return isotope_list;
	}
	
	
	/**
	 * Merge the compound based on the mass tolerance
	 * @param list
	 * @param ppm_cutoff
	 * @return
	 */
	public static LinkedList merge_mass_tolerance(LinkedList list, boolean ppm, double ppm_val) {
		HashMap final_map = new HashMap();
		LinkedList mass_list = new LinkedList();
		LinkedList final_list = new LinkedList();
		Iterator itr = list.iterator();
		while (itr.hasNext()) {
			COMPOUND compound1 = (COMPOUND)itr.next();
			double mass = compound1.MASS.NUM * Math.pow(10, compound1.MASS.E);
			if (final_map.containsKey(mass)) {
				COMPOUND old = (COMPOUND)final_map.get(mass);
				//System.out.println("Merge1: " + compound1.TOTALPROB + "\t" + old.TOTALPROB);
				compound1.TOTALPROB = add(compound1.TOTALPROB, old.TOTALPROB);				
				final_map.put(mass, compound1);
			} else {
				final_map.put(mass, compound1);
				mass_list.add(mass);
			}
		}
		double[] double_list = new double[mass_list.size()];
		int i = 0;
		itr = mass_list.iterator();
		while (itr.hasNext()) {
			double mass = (Double)itr.next();
			double_list[i] = mass;
			i++;
		}
		Arrays.sort(double_list);
		//System.out.println(double_list.length);
		//final_map.clear();
		
		double mass1 = double_list[0];;
		for (i = 1 ;i < double_list.length; i++) {
			
			double mass2 = double_list[i];
			COMPOUND compound1 = (COMPOUND)final_map.get(mass1);
			COMPOUND compound2 = (COMPOUND)final_map.get(mass2);
			boolean update = false;
			if (ppm) {
				//System.out.println(mass1 + "\t" + mass2);
				//if (check_within_ppm(compound1.MASS, compound2.MASS, new BIGNUM(10))) {
				if (check_within_ppm(compound1.MASS, compound2.MASS, new BIGNUM(ppm_val))) {
					update = true;
					
					if (isGreater(compound1.TOTALPROB, compound2.TOTALPROB)) {
						//System.out.println("Merge2: " + compound1.TOTALPROB.getString() + "\t" + compound2.TOTALPROB.getString());
						compound1.TOTALPROB = add(compound1.TOTALPROB, compound2.TOTALPROB);
						final_map.remove(mass2);
						final_map.put(mass1,  compound1);		
						mass1 = mass1;
					} else {
						compound2.TOTALPROB = add(compound1.TOTALPROB, compound2.TOTALPROB);
						final_map.remove(mass1);
						final_map.put(mass2,  compound2);
						mass1 = mass2;
					}
				} else {
					//System.out.println(mass1 + "\t" + mass2);
				}
			} else {
				//if (within_distance(compound1.MASS, compound2.MASS, new BIGNUM(0.1))) {
				if (within_distance(compound1.MASS, compound2.MASS, new BIGNUM(ppm_val))) {
					update = true;
					if (isGreater(compound1.TOTALPROB, compound2.TOTALPROB)) {
						//System.out.println("Merge2: " + compound1.TOTALPROB.getString() + "\t" + compound2.TOTALPROB.getString());
						compound1.TOTALPROB = add(compound1.TOTALPROB, compound2.TOTALPROB);
						final_map.remove(mass2);
						final_map.put(mass1,  compound1);
						mass1 = mass1;
					} else {
						compound2.TOTALPROB = add(compound1.TOTALPROB, compound2.TOTALPROB);
						final_map.remove(mass1);
						final_map.put(mass2,  compound2);
						mass1 = mass2;
					}
					
					//System.out.println("Merge2: " + compound1.TOTALPROB.getString() + "\t" + compound2.TOTALPROB.getString());
					
				} else {
					//System.out.println(mass1 + "\t" + mass2);
				}
			}
			if (!update) {
				mass1 = mass2;
			}
		}
		//System.out.println(final_map.size());
		i = 0;
		double_list = new double[final_map.size()];
		itr = final_map.keySet().iterator();
		while (itr.hasNext()) {
			double mass = (Double)itr.next();
			double_list[i] = mass;
			i++;
		}
		Arrays.sort(double_list);
		for (double mass: double_list) {
			final_list.add(final_map.get(mass));
		}

		return final_list;
	}	
	public static LinkedList condense_filer(LinkedList list) {
		LinkedList short_list = new LinkedList();
		HashMap map = new HashMap();
		Iterator itr = list.iterator();
		while (itr.hasNext()) {
			COMPOUND compound = (COMPOUND)itr.next();
			//if (isGreater(compound.PROB, cutoff)) {
			String key = compound.TYPE1 + "\t" + compound.TYPE2 + "\t" + compound.TYPE3 + "\t" + compound.TYPE4; //compound.NUM_ISOTOPE + "_" + compound.PROB.E + "E" + compound.PROB.NUM;
			if (map.containsKey(key)) {
				COMPOUND old = (COMPOUND)map.get(key);
				compound.FREQ = add(old.FREQ, compound.FREQ);
				//compound.TOTALPROB = multiply(compound.FREQ, compound.PROB);
				compound.TOTALPROB = multiply(compound.TOTALPROB, compound.TOTALPROB);
				map.put(key,  compound);
			} else {
				map.put(key, compound);
				//short_list.add(compound);
			} 
		}
		Iterator itr2 = map.keySet().iterator();
		while (itr2.hasNext()) {
			String key = (String)itr2.next();
			short_list.add(map.get(key));
		}
		return short_list;
	}

	public static void printTime() {
		Calendar cal = Calendar.getInstance();
    	cal.getTime();
    	SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss");
    	System.out.println( sdf.format(cal.getTime()) );
	}
	
	public static void print_list(LinkedList list) {
		Iterator itr = list.iterator();
		while (itr.hasNext()) {
			COMPOUND compound = (COMPOUND)itr.next();
			
			//System.out.println(compound.SIZE + "\t" + compound.NUM_ISOTOPE + "\t" + compound.MASS.getString() + "\t" + compound.PROB.getString() + "\t" + compound.FREQ.getString() + "\t" + multiply(compound.FREQ, compound.PROB).getString());
			System.out.println(compound.SIZE + "\t" + compound.MASS.getString() + "\t" + compound.TOTALPROB.getString());
		}
	}
	
	/**
	 * Subtract two bignum
	 * subtract bn2 from bn1
	 */
	public static BIGNUM subtract(BIGNUM bn1, BIGNUM bn2) {
		BIGNUM newbn = new BIGNUM(-bn2.NUM, bn2.E);
		return add(bn1, newbn);
	}
	
	/**
	 * Merge Resolution 
	 * @param mass1 other mass
	 * @param mass2 theoretical mass
	 * @param ppm_cutoff
	 * @return
	 */
	public static boolean within_distance(BIGNUM mass1, BIGNUM mass2, BIGNUM ppm_cutoff) {
		BIGNUM million = new BIGNUM(1e6);
		BIGNUM difference;
		if (isGreater(mass1, mass2)) {
			difference = subtract(mass1, mass2);
		} else {
			difference = subtract(mass2, mass1);	
		}
		
		
		//System.out.println(mass1.getString() + "\t" + mass2.getString() + "\t" + ppm.getString());
		if (isGreater(ppm_cutoff, difference)) {
			//System.out.println("Is greater: " + mass1.getString() + "\t" + mass2.getString() + "\t" + ppm.getString() + "\t" + ppm_cutoff.getString());
			
			return true;
		}
		//System.out.println("Is smaller: " + mass1.getString() + "\t" + mass2.getString() + "\t" + ppm.getString() + "\t" + ppm_cutoff.getString());
		return false;
	}
	/**
	 * PPM = 1E6 * (Mass1 - Mass2) / Mass2 
	 * @param mass1 other mass
	 * @param mass2 theoretical mass
	 * @param ppm_cutoff
	 * @return
	 */
	public static boolean check_within_ppm(BIGNUM mass1, BIGNUM mass2, BIGNUM ppm_cutoff) {
		BIGNUM million = new BIGNUM(1e6);
		BIGNUM difference;
		if (isGreater(mass1, mass2)) {
			difference = subtract(mass1, mass2);
		} else {
			difference = subtract(mass2, mass1);	
		}
		BIGNUM ppm = divide(multiply(million, difference), divide(add(mass2, mass1), new BIGNUM(2)));
		
		//System.out.println(mass1.getString() + "\t" + mass2.getString() + "\t" + ppm.getString());
		if (isGreater(ppm_cutoff, ppm)) {
			//System.out.println("Is greater: " + mass1.getString() + "\t" + mass2.getString() + "\t" + ppm.getString() + "\t" + ppm_cutoff.getString());
			
			return true;
		}
		//System.out.println("Is smaller: " + mass1.getString() + "\t" + mass2.getString() + "\t" + ppm.getString() + "\t" + ppm_cutoff.getString());
		return false;
	}
	/*public static COMPOUND merge_compound(COMPOUND compound_A, COMPOUND compound_B) {
		COMPOUND new_compound = new COMPOUND();
		//System.out.println(compound_A.TYPE + "\t" + compound_B.TYPE);
		if (compound_A.TYPE.equals(compound_B.TYPE)) {
			new_compound.PROB = multiply(compound_A.PROB, compound_B.PROB);
			new_compound.SIZE = compound_A.SIZE + compound_B.SIZE;
			//new_compound.NUM_ISOTOPE = compound_A.NUM_ISOTOPE + compound_B.NUM_ISOTOPE;
			new_compound.TYPE1 = compound_A.TYPE1 + compound_B.TYPE1;
			new_compound.TYPE2 = compound_A.TYPE2 + compound_B.TYPE2;
			new_compound.TYPE3 = compound_A.TYPE3 + compound_B.TYPE3;
			new_compound.FREQ = multiply(compound_A.FREQ, compound_B.FREQ); //calc_N_choose_M(new_compound.SIZE, new_compound.NUM_ISOTOPE);
			new_compound.TYPE = compound_A.TYPE;
			new_compound.MASS = add(compound_A.MASS, compound_B.MASS);
			//new_compound.TOTALPROB = multiply(new_compound.FREQ, new_compound.PROB);
			new_compound.TOTALPROB = multiply(compound_A.TOTALPROB, compound_B.TOTALPROB);
		}
		return new_compound;
	}*/
	
	public static COMPOUND merge_compound(COMPOUND compound_A, COMPOUND compound_B) {
		COMPOUND new_compound = new COMPOUND();
		
		new_compound.NUM_ISOTOPE = compound_A.NUM_ISOTOPE + compound_B.NUM_ISOTOPE;
		new_compound.SIZE = compound_A.SIZE + compound_B.SIZE;
		new_compound.MASS = add(compound_A.MASS, compound_B.MASS);
		new_compound.TOTALPROB = multiply(compound_A.TOTALPROB, compound_B.TOTALPROB);
		return new_compound;
	}
	/**
	 * Merging two different elements and return a merged list
	 * @param element_A
	 * @param element_B
	 * @return
	 */
	public static LinkedList MergeList(LinkedList element_A, LinkedList element_B) {
		LinkedList merged = new LinkedList();
		if (element_B.size() == 0) {
			return element_A;
		}
		if (element_A.size() == 0) {
			return element_B;
		}
		Iterator itr = element_A.iterator();
		while (itr.hasNext()) {
			COMPOUND compound_A = (COMPOUND)itr.next();
			Iterator itr2 = element_B.iterator();
			while (itr2.hasNext()) {
				COMPOUND compound_B = (COMPOUND)itr2.next();				
				COMPOUND new_compound = merge_compound(compound_A, compound_B);
				
				// add whether to save this value
				if (isGreater(new_compound.TOTALPROB, cutoff)) {
					merged.add(new_compound);
				}
			}
		}
		return merged;
	}
	
	/**
	 * Grab Table
	 * @param fileName
	 * @param size
	 * @param type
	 * @return
	 */
	public static LinkedList grabTable(String fileName, int size, String type) {
		LinkedList list = new LinkedList();
		try {

			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				int num = new Integer(split[0]);
				COMPOUND compound = new COMPOUND();
				//compound.NUM_ISOTOPE = num - 1;
				compound.MASS = new BIGNUM(split[1]);
				compound.PROB = new BIGNUM(split[2]);
				compound.FREQ = new BIGNUM(split[3]);
				compound.SIZE = size;
				compound.TYPE = type;
				list.add(compound);
			}
			in.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return list;
	}

	public static BIGNUM total_val(int N, int M) {
		BIGNUM val = new BIGNUM(0.0);
		for (int i = 0; i <= M; i++) {
			val = add(val, calc_N_choose_M(N, i));
		}
		return val;
	}
	
	public static BIGNUM calc_N_choose_M(int N, int M) {
		BIGNUM top = new BIGNUM(new Double(N));
		BIGNUM bottom = new BIGNUM(1.0);
		if (M == 0) {
			return new BIGNUM(1.0);
		}
		for (int j = 1; j <= M; j++) {
			bottom = multiply(bottom, new BIGNUM(new Double(j)));
		}
		for (int j = 1; j < M; j++) {
			top = multiply(top, new BIGNUM(new Double(N - j)));
		}		
		return divide(top, bottom);
	}
	/**
	 * Test if first number is larger than second number
	 * @param num
	 * @param num2
	 * @return
	 */
	public static boolean isGreater(BIGNUM num, BIGNUM num2) {
		if (num.E > num2.E) {
			return true;
		} else if (num.E == num2.E) {
			if (num.NUM > num2.NUM) {
				return true;
			}
		}
		return false;
	}
	/**
	 * Dividing two big num
	 * @param bn1
	 * @param bn2
	 * @return
	 */
	public static BIGNUM divide(BIGNUM bn1, BIGNUM bn2) {
		double divide = bn1.NUM / bn2.NUM;
		double newE = bn1.E - bn2.E;
		//System.out.println("divide function: " + multiply + "\t" + newE);
		BIGNUM result = new BIGNUM(divide, newE);
		return result;
	}	
	/**
	 * Multiplying two big num
	 * @param bn1
	 * @param bn2
	 * @return
	 */
	public static BIGNUM multiply(BIGNUM bn1, BIGNUM bn2) {
		double multiply = bn1.NUM * bn2.NUM;
		double newE = bn1.E + bn2.E;
		//System.out.println("multiply function: " + multiply + "\t" + newE);
		BIGNUM result = new BIGNUM(multiply, newE);
		return result;
	}
	
	/**
	 * Calculate the base factor 10 difference between two number
	 * num1 should be greater than num2
	 * @param num1
	 * @param num2
	 * @return
	 */
	public static double DIFF(double num1, double num2) {
		double DIFF = 0;
		double value = num1 / num2;
		for (double i = 10; i < 1e50; i = i * 10) {
			if (value / i < 1 || value == 0) {											
				break;
			} else {
				
				DIFF++;
			}
		}
		return DIFF;
	}
	/**
	 * Calculate the power of a number
	 * @param val
	 * @param num
	 * @return
	 */
	public static BIGNUM pow(BIGNUM val, int num) {
		BIGNUM result = val;
		if (num == 1) {
			return result;
		}
		if (num == 0) {
			if (val.NUM == 0) {
				return null;
			}
			return new BIGNUM(1.0);
		}
		for (int i = 1; i < num; i++) {
			result = multiply(result, val);
		}
		return result;
	}
	/**
	 * ADD two big num
	 * If the factor between the two number is greater than 52 then just return the larger number
	 * @param bn1
	 * @param bn2
	 * @return
	 */
	public static BIGNUM add(BIGNUM bn1, BIGNUM bn2) {
		int diff = 0;
		double COMB = 0;
		double newE = 0;
		if (bn1.E > bn2.E) {
			
			double factor = bn1.E - bn2.E;
			if (factor > 52 ) {
				return bn1;
			}
			COMB = bn1.NUM + bn2.NUM / Math.pow(10, factor);
			newE = bn1.E;
		} else {
			double factor = bn2.E - bn1.E;
			if (factor > 52 ) {
				return bn2;
			}
			COMB = bn2.NUM + bn1.NUM / Math.pow(10, factor);
			newE = bn2.E;
		}
		
		BIGNUM result = new BIGNUM(COMB, newE);
		return result;
	}
	
	/**
	 * Function for testing the BIGNUM function
	 */
	public static void TESTBIGNUM() {
		
		BIGNUM num = new BIGNUM(0);
		System.out.println("Test 0 should be 0 and 0");
		System.out.println(num.NUM + "\t" + num.E);
		
		num = new BIGNUM(1);
		System.out.println("Test 1 should be 1.0 and 0");
		System.out.println(num.NUM + "\t" + num.E);
		
		num = new BIGNUM(10);
		System.out.println("Test 10 should be 1.0 and 1");
		System.out.println(num.NUM + "\t" + num.E);
		
		num = new BIGNUM(5);
		System.out.println("Test 5 should be 5 and 0");
		System.out.println(num.NUM + "\t" + num.E);
		
		num = new BIGNUM(51);
		System.out.println("Test 51 should be 5.1 and 1");
		System.out.println(num.NUM + "\t" + num.E);
		
		num = new BIGNUM(9351);
		System.out.println("Test 9351 should be 9.351 and 3");
		System.out.println(num.NUM + "\t" + num.E);
		
		num = new BIGNUM(0.5);
		System.out.println("Test 0.5 should be 5 and -1");
		System.out.println(num.NUM + "\t" + num.E);
		
		num = new BIGNUM(0.145);
		System.out.println("Test 0.145 should be 1.45 and -1");
		System.out.println(num.NUM + "\t" + num.E);
		
		num = new BIGNUM(0.000008145);
		System.out.println("Test 0.000008145 should be 8.145 and -6");
		System.out.println(num.NUM + "\t" + num.E);
		
		// try multiplication function
		BIGNUM num1 = new BIGNUM(4);		
		BIGNUM num2 = new BIGNUM(5);
		
		BIGNUM multiply = multiply(num1, num2);
		System.out.println("Multiply 4 and 5 should get 2.0 and 1");
		System.out.println(multiply.NUM + "\t" + multiply.E);
		
		// try multiplication function
		num1 = new BIGNUM(100.5);
		num2 = new BIGNUM(401);
		
		multiply = multiply(num1, num2);
		System.out.println("Multiply 100.5 and 401 should get 4.03005 and 4");
		System.out.println(multiply.NUM + "\t" + multiply.E);
		
		
		// test multiplication case 2
		num1 = new BIGNUM(100878098.5);
		num2 = new BIGNUM(0.000000413124);		
		multiply = multiply(num1, num2);
		System.out.println("Multiply 100878098.5 and 0.000000413124 should get 4.03005 and 4");
		System.out.println(multiply.NUM + "\t" + multiply.E);
		
		// test addition 1
		num1 = new BIGNUM(1);
		num2 = new BIGNUM(2);		
		BIGNUM add = add(num1, num2);
		System.out.println("Add 1 and 2 should get 3 and 0");
		System.out.println(add.NUM + "\t" + add.E);

		// test addition 2
		num1 = new BIGNUM(10);
		num2 = new BIGNUM(2);		
		add = add(num1, num2);
		System.out.println("Add 10 and 2 should get 1.2 and 1");
		System.out.println(add.NUM + "\t" + add.E);
		
		// test addition 3
		num1 = new BIGNUM(0.004);
		num2 = new BIGNUM(2234);		
		add = add(num1, num2);
		System.out.println("Add 0.004 and 2234 should get 2.234004 and 0");
		System.out.println(add.NUM + "\t" + add.E);
		
		// test addition 4
		num1 = new BIGNUM(2234);
		num2 = new BIGNUM(0.0045);		
		add = add(num1, num2);
		System.out.println("Add 2234 and 0.0045should get 2.234004 and 0");
		System.out.println(add.NUM + "\t" + add.E);
		
		System.exit(0);
	}

	/**
	 * A Class for storing compound information
	 * @author tshaw
	 * 
	 */
	static class COMPOUND {
		public String TYPE = "";
		public int TYPE1 = 0; // isotope cumulative types
		public int TYPE2 = 0; // isotope cumulative types
		public int TYPE3 = 0; // isotope cumulative types
		public int TYPE4 = 0; // isotope cumulative types
		public int TYPE5 = 0; // isotope cumulative types
		public BIGNUM MASS = new BIGNUM(0.0);
		public BIGNUM PROB = new BIGNUM(0.0);
		public BIGNUM FREQ = new BIGNUM(0.0);
		public BIGNUM TOTALPROB = new BIGNUM(0.0);		
		public int NUM_ISOTOPE = 0; // the number of isotopes
		public int SIZE = 0;
		public void SETTYPE(int ID) {
			if (ID == 0) {
				TYPE1 = 1;
			}
			if (ID == 1) {
				TYPE2 = 1;
			}
			if (ID == 2) {
				TYPE3 = 1;
			}
			if (ID == 3) {
				TYPE4 = 1;
			}
			if (ID == 4) {
				TYPE5 = 1;
			}
			
		}
	}
	/**
	 * Parsing the string and retrieve the number of elements
	 * @param str
	 * @param query_element
	 * @return
	 */
	public static int retrieve_num_element(String str, String query_element) {
		int total = 0;
		int str_index = str.indexOf(query_element);
			
		for (int i = 0; i < str.length(); i++) {
			if (str.substring(i, i + 1).equals(query_element)) {
				int value = 0;
				String num_str = "";
				for (int j = i + 1; j < str.length(); j++) {
					
					if (str.substring(j, j + 1).matches("[0-9]")) {
						num_str += str.substring(j, j + 1);
						value = new Integer(num_str);
					} else {
						
						if (value == 0) {
							value = 1;
						}
						break;
					}
				}
				if (i == str.length() - 1) {
					value = 1;
				}
				total += value;
			}
		}

		return total;
	}
	/**
	 * A Class for storing large numbers;
	 * @author tshaw
	 *
	 */
	static class BIGNUM {
		
		public double NUM = 0;
		public double E = 0;
		
		public String getString() {
			return (new Double(NUM).toString() + "E" + new Integer((int)E).toString());
		}
		public double toDouble() {
			return new Double(this.getString());
		}
		public BIGNUM(double newNUM, double newE) {
			NUM = newNUM;
			E = newE;
			double value = newNUM;
			
			if (value >= 1) {
				
				for (double i = 10; i < 1e50; i = i * 10) {
					if (value / i < 1 || value == 0) {											
						break;
					} else {
						NUM = value / i;
						E++;
					}
				}
			} else {
				for (double i = 10; i < 1e50; i = i * 10) {
					if (value * i > 10 || value == 0) {
						break;
					} else {
						NUM = value * i;
						E--;
					}
				}
			}
		}
		public BIGNUM(String value2) {
			//System.out.println(value);
			if (value2.contains("E")) {
				String[] split = value2.split("E");
				NUM = new Double(split[0]);
				E = new Double(split[1]);
				double value = new Double(split[0]);
				
				if (value >= 1) {
					
					for (double i = 10; i < 1e50; i = i * 10) {
						if (value / i < 1 || value == 0) {											
							break;
						} else {
							NUM = value / i;
							E++;
						}
					}
				} else {
					for (double i = 10; i < 1e50; i = i * 10) {
						if (value * i > 10 || value == 0) {
							break;
						} else {
							NUM = value * i;
							E--;
						}
					}
				}
			} else {
				double value = new Double(value2);
				NUM = new Double(value);
				if (value >= 1) {
					E = 0;
					for (double i = 10; i < 1e50; i = i * 10) {
						if (value / i < 1 || value == 0) {											
							break;
						} else {
							NUM = value / i;
							E++;
						}
					}
				} else {
					for (double i = 10; i < 1e50; i = i * 10) {
						if (value * i > 10 || value == 0) {
							break;
						} else {
							NUM = value * i;
							E--;
						}
					}
				}
			}
		}
		public BIGNUM(double value) {
			//System.out.println(value);
			NUM = value;
			if (value >= 1) {
				E = 0;
				for (double i = 10; i < 1e50; i = i * 10) {
					if (value / i < 1 || value == 0) {											
						break;
					} else {
						NUM = value / i;
						E++;
					}
				}
			} else {
				for (double i = 10; i < 1e50; i = i * 10) {
					if (value * i > 10 || value == 0) {
						break;
					} else {
						NUM = value * i;
						E--;
					}
				}
			}
		}
	}
	public static boolean checkInvalidCharacter(String str) {
		if (str.contains("Cl")) {
			return false;
		}
		if (str.contains("Mg")) {
			return false;
		}
		if (str.contains("V")) {
			return false;
		}
		if (str.contains("He")) {
			return false;
		}
		if (str.contains("B")) {
			return false;
		}
		if (str.contains("Be")) {
			return false;
		}
		if (str.contains("Si")) {
			return false;
		}
		if (str.contains("Na")) {
			return false;
		}
		if (str.contains("F")) {
			return false;
		}
		if (str.contains("K")) {
			return false;
		}
		if (str.contains("Ar")) {
			return false;
		}
		if (str.contains("Ca")) {
			return false;
		}
		if (str.contains("Ti")) {
			return false;
		}
		if (str.contains("Cr")) {
			return false;
		}
		if (str.contains("Mn")) {
			return false;
		}
		if (str.contains("Fe")) {
			return false;
		}
		if (str.contains("Ni")) {
			return false;
		}
		if (str.contains("Co")) {
			return false;
		}
		if (str.contains("Cu")) {
			return false;
		}
		if (str.contains("Zn")) {
			return false;
		}
		if (str.contains("Ge")) {
			return false;
		}
		if (str.contains("As")) {
			return false;
		}
		if (str.contains("Br")) {
			return false;
		}
		if (str.contains("Se")) {
			return false;
		}
		if (str.contains("Rb")) {
			return false;
		}
		if (str.contains("Sr")) {
			return false;
		}
		if (str.contains("Zr")) {
			return false;
		}
		if (str.contains("Mo")) {
			return false;
		}
		if (str.contains("Ru")) {
			return false;
		}
		if (str.contains("Pd")) {
			return false;
		}
		if (str.contains("Ag")) {
			return false;
		}
		if (str.contains("Cd")) {
			return false;
		}
		if (str.contains("Sb")) {
			return false;
		}
		if (str.contains("I")) {
			return false;
		}
		if (str.contains("Te")) {
			return false;
		}
		if (str.contains("Cs")) {
			return false;
		}
		if (str.contains("Ba")) {
			return false;
		}
		if (str.contains("Ce")) {
			return false;
		}
		if (str.contains("Nd")) {
			return false;
		}
		if (str.contains("Gd")) {
			return false;
		}
		if (str.contains("Hf")) {
			return false;
		}
		if (str.contains("Ta")) {
			return false;
		}
		if (str.contains("W")) {
			return false;
		}
		if (str.contains("Re")) {
			return false;
		}
		if (str.contains("Pt")) {
			return false;
		}
		if (str.contains("Au")) {
			return false;
		}
		if (str.contains("Hg")) {
			return false;
		}
		if (str.contains("Pb")) {
			return false;
		}
		if (str.contains("Bi")) {
			return false;
		}
		if (str.contains("Th")) {
			return false;
		}
		if (str.contains("ClTl")) {
			return false;
		}
		return true;
	}
}
