package ISOTOPEDISTRIBUTION.MS1;

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

import ISOTOPEDISTRIBUTION.IsotopeCalculator;
import MISC.ToolBox;


/**
 * calculates the isotope pattern
 * @author tshaw
 *
 */
public class CalculateDegradedPeptideIsotopePattern {

	/*public static double A = -1;
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
	*/
	public static double C12 = 0.9893;
	public static double N14 = 0.99632;
	public static double H1 = 0.999885;

	public static Random rand = new Random();
	
	public static void execute(String[] args) {
		
		try {
			

			
			//String outputFile = "C:\\Users\\tshaw\\Desktop\\PROTEOMICS\\PeptideSimulation\\simulation.txt";
			int i = 0;
			String inputFile = args[0];
			String outputFile = args[1];
			String type = args[2];
			FileWriter fwriter = new FileWriter(outputFile);
			BufferedWriter out = new BufferedWriter(fwriter);
			
			LinkedList list = readPeptideList(inputFile);
			Iterator itr = list.iterator();
			while (itr.hasNext()) {
				
			
					//String peptideFormula = generateRandPeptide(i);
				String peptideFormula = (String)itr.next();
				i++;
				String chemFormula = peptide2formula(peptideFormula);
				int number_of_carbons = ToolBox.retrieve_num_element(chemFormula, "C");
				int number_of_oxygen = ToolBox.retrieve_num_element(chemFormula, "O");
				int number_of_nitrogen = ToolBox.retrieve_num_element(chemFormula, "N");
				//System.out.println(peptideFormula + "\t" + chemFormula);
				
				if (type.equals("TMT")) {
					chemFormula += "C8H20Ci4N1Nm1O2";
				}
				IsotopePattern pattern = IsotopeCalculator.calculate_pattern(chemFormula, 50, C12, H1, N14, 1);
				
				
				double highestMass = getMassOfHighestPeak(pattern);
				double firstMass = getFirstMass(pattern);
				
				double highestIntensity = getIntensityOfHighestPeak(pattern);
				double firstIntensity = getFirstIntensity(pattern);
				
				/*for (int k = 0; k < pattern.getNumberOfIsotopes(); k++) {
					IsotopeContainer container = pattern.getIsotope(k);
					
					out2.write(container.getMass() + "\t" + container.getIntensity() + "\t" + new Double(highestMass - firstMass).intValue() + "\n");
					//System.out.println(container.getMass() + "\t" + container.getIntensity());
				}*/
				
				
				String stuff = i + "\t" + firstIntensity + "\t" + highestIntensity + "\t" + (highestIntensity / firstIntensity) 
						+ "\t" + firstMass + "\t" + highestMass + "\t" + new Double(highestMass - firstMass).intValue() + "\t" + peptideFormula + "\t" 
						+ chemFormula + "\t" + number_of_carbons + "\t" + number_of_oxygen 
						+ "\t" + number_of_nitrogen	+ "\t" + (number_of_carbons + number_of_oxygen) + "\t" 
						+ (number_of_carbons + number_of_nitrogen) + "\t" + (number_of_carbons + number_of_oxygen + number_of_nitrogen) + "\t" + getStrIsotope(pattern);
				//System.out.println(firstMass + "\t" + highestMass + "\t" + (highestMass - firstMass) + "\t" + peptideFormula + "\t" 
				//		+ chemFormula + "\t" + number_of_carbons + "\t" + number_of_oxygen 
				//		+ "\t" + number_of_nitrogen	+ "\t" + getStrIsotope(pattern));
				out.write(stuff + "\n");
				if (i % 1000 == 0) {
					out.flush();
				}
			}
			out.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public static void main(String[] args) {
		
		try {
			
			//String outputFile = "C:\\Users\\tshaw\\Desktop\\PROTEOMICS\\PeptideSimulation\\simulation.txt";
			int i = 0;
			String inputFile = args[0];
			String outputFile = args[1];
			FileWriter fwriter = new FileWriter(outputFile);
			BufferedWriter out = new BufferedWriter(fwriter);
			
			LinkedList list = readPeptideList(inputFile);
			Iterator itr = list.iterator();
			while (itr.hasNext()) {
				
			
					//String peptideFormula = generateRandPeptide(i);
				String peptideFormula = (String)itr.next();
				i++;
				String chemFormula = peptide2formula(peptideFormula);
				int number_of_carbons = ToolBox.retrieve_num_element(chemFormula, "C");
				int number_of_oxygen = ToolBox.retrieve_num_element(chemFormula, "O");
				int number_of_nitrogen = ToolBox.retrieve_num_element(chemFormula, "N");
				//System.out.println(peptideFormula + "\t" + chemFormula);
				IsotopePattern pattern = IsotopeCalculator.calculate_pattern(chemFormula, 50, C12, H1, N14, 1);
				
				
				double highestMass = getMassOfHighestPeak(pattern);
				double firstMass = getFirstMass(pattern);
				
				double highestIntensity = getIntensityOfHighestPeak(pattern);
				double firstIntensity = getFirstIntensity(pattern);
				
				/*for (int k = 0; k < pattern.getNumberOfIsotopes(); k++) {
					IsotopeContainer container = pattern.getIsotope(k);
					
					out2.write(container.getMass() + "\t" + container.getIntensity() + "\t" + new Double(highestMass - firstMass).intValue() + "\n");
					//System.out.println(container.getMass() + "\t" + container.getIntensity());
				}*/
				
				
				String stuff = i + "\t" + firstIntensity + "\t" + highestIntensity + "\t" + (highestIntensity / firstIntensity) 
						+ "\t" + firstMass + "\t" + highestMass + "\t" + new Double(highestMass - firstMass).intValue() + "\t" + peptideFormula + "\t" 
						+ chemFormula + "\t" + number_of_carbons + "\t" + number_of_oxygen 
						+ "\t" + number_of_nitrogen	+ "\t" + (number_of_carbons + number_of_oxygen) + "\t" 
						+ (number_of_carbons + number_of_nitrogen) + "\t" + (number_of_carbons + number_of_oxygen + number_of_nitrogen) + "\t" + getStrIsotope(pattern);
				//System.out.println(firstMass + "\t" + highestMass + "\t" + (highestMass - firstMass) + "\t" + peptideFormula + "\t" 
				//		+ chemFormula + "\t" + number_of_carbons + "\t" + number_of_oxygen 
				//		+ "\t" + number_of_nitrogen	+ "\t" + getStrIsotope(pattern));
				out.write(stuff + "\n");
				if (i % 1000 == 0) {
					out.flush();
				}
			}
			out.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static LinkedList readPeptideList(String fileName) {
		
		LinkedList list = new LinkedList();
		try {
			String seq = "";
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				list.add(str);
			}
			in.close();
			list.add(seq);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return list;
	}
	
	public static double getFirstMass(IsotopePattern pattern) {
		IsotopeContainer container = pattern.getIsotope(0);
		return container.getMass();
	}
	
	public static double getMassOfHighestPeak(IsotopePattern pattern) {
		String returnStr = pattern.getNumberOfIsotopes() + "";
		double max_intensity = -1;
		double max_mass = -1;
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
			IsotopeContainer container = pattern.getIsotope(i);
			if (max_intensity < container.getIntensity()) {
				max_intensity = container.getIntensity();
				max_mass = container.getMass();
			}
			//System.out.println(container.getMass() + "\t" + container.getIntensity());				
		}
		return max_mass;
	}
	
	public static double getFirstIntensity(IsotopePattern pattern) {
		IsotopeContainer container = pattern.getIsotope(0);
		return container.getIntensity();
	}
	
	public static double getIntensityOfHighestPeak(IsotopePattern pattern) {
		String returnStr = pattern.getNumberOfIsotopes() + "";
		double max_intensity = -1;
		double max_mass = -1;
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
			IsotopeContainer container = pattern.getIsotope(i);
			if (max_intensity < container.getIntensity()) {
				max_intensity = container.getIntensity();
				max_mass = container.getMass();
			}
			//System.out.println(container.getMass() + "\t" + container.getIntensity());				
		}
		return max_intensity;
	}
	public static String getStrIsotope(IsotopePattern pattern) {
		int countNumIsotopes = 0;
		double cutoff = 0.00001;
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
			IsotopeContainer container = pattern.getIsotope(i);
			if (container.getIntensity() > cutoff) {
				countNumIsotopes++;
			}
		}
		String returnStr = countNumIsotopes + "";
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
			IsotopeContainer container = pattern.getIsotope(i);
			if (container.getIntensity() > cutoff) {
				if (returnStr.equals("")) {
					returnStr += container.getMass() + ":" + container.getIntensity();
				} else {
					returnStr += "," + container.getMass() + ":" + container.getIntensity();
				}
			}
			//System.out.println(container.getMass() + "\t" + container.getIntensity());				
		}
		return returnStr;
	}
	public static String peptide2formula(String str) {
		str = str.replaceAll("FORMULA",  "").trim();
		str = str.replaceAll("=", "").trim();
		MOLECULE mol = new MOLECULE();
		// Grab Alanine A
		int num = ToolBox.retrieve_num_element(str, "A");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("A"));
		}
		// Grab Arginine R
		num = ToolBox.retrieve_num_element(str, "R");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("R"));
		}
		// Grab Asaparagine N
		num = ToolBox.retrieve_num_element(str, "N");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("N"));
		}
		// Grab Aspartic acid D
		num = ToolBox.retrieve_num_element(str, "D");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("D"));
		}
		// Grab Cysteine C
		num = ToolBox.retrieve_num_element(str, "C");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("C"));
		}
		// Grab Glutamine Q
		num = ToolBox.retrieve_num_element(str, "Q");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("Q"));
		}
		// Grab Glutamic acid E
		num = ToolBox.retrieve_num_element(str, "E");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("E"));
		}
		// Grab Glycine G
		num = ToolBox.retrieve_num_element(str, "G");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("G"));
		}
		// Grab Histidine H
		num = ToolBox.retrieve_num_element(str, "H");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("H"));
		}
		// Grab Isoleucine I
		num = ToolBox.retrieve_num_element(str, "I");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("I"));
		}
		// Grab Leucine L
		num = ToolBox.retrieve_num_element(str, "L");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("L"));
		}
		// Grab Lysine K
		num = ToolBox.retrieve_num_element(str, "K");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("K"));
		}
		// Grab Methionine M
		num = ToolBox.retrieve_num_element(str, "M");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("M"));
		}
		// Grab Phenylalanine F
		num = ToolBox.retrieve_num_element(str, "F");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("F"));
		}
		// Grab Proline P
		num = ToolBox.retrieve_num_element(str, "P");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("P"));
		}
		// Grab Serine S
		num = ToolBox.retrieve_num_element(str, "S");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("S"));
		}
		// Grab Threonine T
		num = ToolBox.retrieve_num_element(str, "T");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("T"));
		}
		// Grab Tryptophan W
		num = ToolBox.retrieve_num_element(str, "W");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("W"));
		}
		// Grab Tyrosine Y
		num = ToolBox.retrieve_num_element(str, "Y");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("Y"));
		}
		// Grab Valine V
		num = ToolBox.retrieve_num_element(str, "V");
		for (int i = 0; i < num; i++) {
			mol = concatenate_molecule(mol, ConvertAAtoFormula("V"));
		}
		
		
		// add an additional of OH group to the concatenated residues.
		mol.O++;
		mol.H++;
		
		
		String finalStr = "";
		if (mol.C > 0) {
			finalStr += "C" + mol.C;
		}
		if (mol.H > 0) {
			finalStr += "H" + (new Integer(mol.H) + 1);
		}
		if (mol.N > 0) {
			finalStr += "N" + mol.N;
		}
		if (mol.O > 0) {
			finalStr += "O" + mol.O;
		}
		if (mol.S > 0) {
			finalStr += "S" + mol.S;
		}
		return finalStr;
		
		/*C_num = mol.C;
		H_num = mol.H;
		N_num = mol.N;
		O_num = mol.O;
		S_num = mol.S;*/
		
	}
	public static MOLECULE concatenate_molecule(MOLECULE mol1, MOLECULE mol2) {
		MOLECULE final_mol = new MOLECULE();
		final_mol.C = mol1.C + mol2.C;
		final_mol.H = mol1.H + mol2.H;
		final_mol.N = mol1.N + mol2.N;
		final_mol.S = mol1.S + mol2.S;
		final_mol.O = mol1.O + mol2.O;
		return final_mol;
	}
	/**
	 * data class for molecule/peptide
	 */
	static class MOLECULE {
		public int C = 0;
		public int H = 0;
		public int N = 0;
		public int S = 0;
		public int O = 0;
		public String getFormula() {
			String formula = "";
			if (C > 0) {
				formula += "C" + C;
			}
			if (H > 0) {
				formula += "H" + H;
			}
			if (N > 0) {
				formula += "N" + N;
			}
			if (S > 0) {
				formula += "S" + S;
			}
			if (O > 0) {
				formula += "O" + O;
			}			
			return formula; 
		}
	}
	/**
	 * Formula conversion based on 
	 * http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/
	 * @param aa
	 * @return
	 */
	public static MOLECULE ConvertAAtoFormula(String aa) {
		MOLECULE mol = new MOLECULE();
		if (aa.equals("A")) {
			//C3H7NO2
			mol.C += 3;
			mol.H += 5;
			mol.N += 1;
			mol.O += 1;
		}
		if (aa.equals("R")) {
			//C6H14N4O2
			mol.C += 6;
			mol.H += 12;
			mol.N += 4;
			mol.O += 1;
		}
		if (aa.equals("N")) {
			//C4H8N2O3
			mol.C += 4;
			mol.H += 6;
			mol.N += 2;
			mol.O += 2;
		}
		if (aa.equals("D")) {
			//C4H7NO4
			mol.C += 4;
			mol.H += 5;
			mol.N += 1;
			mol.O += 3;
		}
		if (aa.equals("C")) {
			//C3H7NO2S
			mol.C += 3;
			mol.H += 5;
			mol.N += 1;
			mol.O += 1;
			mol.S += 1;
		}
		if (aa.equals("Q")) {
			//C5H10N2O3
			mol.C += 5;
			mol.H += 8;
			mol.N += 2;
			mol.O += 2;			
		}
		if (aa.equals("E")) {
			//C5H9NO4
			mol.C += 5;
			mol.H += 7;
			mol.N += 1;
			mol.O += 3;			
		}
		if (aa.equals("G")) {
			//C2H5NO2
			mol.C += 2;
			mol.H += 3;
			mol.N += 1;
			mol.O += 1;			
		}
		if (aa.equals("H")) {
			//C6H9N3O2
			mol.C += 6;
			mol.H += 7;
			mol.N += 3;
			mol.O += 1;			
		}
		if (aa.equals("I") || aa.equals("L")) {
			//C6H13NO2
			mol.C += 6;
			mol.H += 11;
			mol.N += 1;
			mol.O += 1;			
		}
		if (aa.equals("K")) {
			//C6H14N2O2
			mol.C += 6;
			mol.H += 12;
			mol.N += 2;
			mol.O += 1;			
		}
		if (aa.equals("M")) {
			//C5H11NO2S
			mol.C += 5;
			mol.H += 9;
			mol.N += 1;
			mol.O += 1;
			mol.S += 1;
		}	
		if (aa.equals("F")) {
			//C9H11NO2
			mol.C += 9;
			mol.H += 9;
			mol.N += 1;
			mol.O += 1;			
		}
		if (aa.equals("P")) {
			//C5H9NO2
			mol.C += 5;
			mol.H += 7;
			mol.N += 1;
			mol.O += 1;			
		}
		if (aa.equals("S")) {
			//C3H7NO3
			mol.C += 3;
			mol.H += 5;
			mol.N += 1;
			mol.O += 2;			
		}
		if (aa.equals("T")) {
			//C4H9NO3
			mol.C += 4;
			mol.H += 7;
			mol.N += 1;
			mol.O += 2;			
		}
		if (aa.equals("W")) {
			//C11H12N2O2
			mol.C += 11;
			mol.H += 10;
			mol.N += 2;
			mol.O += 1;			
		}
		if (aa.equals("Y")) {
			//C9H11NO3
			mol.C += 9;
			mol.H += 9;
			mol.N += 1;
			mol.O += 2;			
		}
		if (aa.equals("V")) {
			//C5H11NO2
			mol.C += 5;
			mol.H += 9;
			mol.N += 1;
			mol.O += 1;			
		}
		return mol;
	}
	public static String generateRandPeptide(int size) {
		String result = "";
		HashMap map = new HashMap();
		for (int i = 0; i < size; i++) {
					
			int n = rand.nextInt(20);
			//System.out.println(n);
			String aa = convertInt2AA(n);
			if (map.containsKey(aa)) {
				int num = (Integer)map.get(aa);
				num++;
				map.put(aa, num);
			} else {
				map.put(aa, 1);
			}
		}
		Iterator itr = map.keySet().iterator();
		while (itr.hasNext()) {
			String key = (String)itr.next();
			int num = (Integer)map.get(key);
			result += key + num;
		}
		return result;
	}
	public static String convertInt2AA(int num) {
		
		String result = "";
		switch (num) {
		case 0:
			result = "R";
			break;
		case 1:
			result = "H";
			break;
		case 2:
			result = "K";
			break;
		case 3:
			result = "D";
			break;
		case 4:
			result = "E";
			break;
		case 5:
			result = "S";
			break;
		case 6:
			result = "T";
			break;
		case 7:
			result = "N";
			break;
		case 8:
			result = "Q";
			break;
		case 9:
			result = "C";
			break;
		case 10:
			result = "G";
			break;
		case 11:
			result = "P";
			break;
		case 12:
			result = "A";
			break;
		case 13:
			result = "V";
			break;
		case 14:
			result = "I";
			break;
		case 15:
			result = "L";
			break;
		case 16:
			result = "M";
			break;
		case 17:
			result = "F";
			break;
		case 18:
			result = "Y";
			break;
		case 19:
			result = "W";
			break;
		case 20:
			result = "U";
			break;
		}
		return result;
	}
	
	public static BIGNUM cutoff = new BIGNUM(1, -50);

	/**
	 * Calculate the pattern for a formula
	 * @param C12
	 * @param N14
	 * @param H1
	 * @return
	 */
	/*public static IsotopePattern calculate_pattern(String formula, int charge) {
		
		C_num = retrieve_num_element(formula, "C");
		H_num = retrieve_num_element(formula, "H");
		N_num = retrieve_num_element(formula, "N");
		O_num = retrieve_num_element(formula, "O");
		S_num = retrieve_num_element(formula, "S");
		P_num = retrieve_num_element(formula, "P");
		
		double C12 = 0.9893;
		double N14 = 0.99632;
		double H1 = 0.999885;
		
		LinkedList ini_list = create_ini_list(C12, N14, H1);
		
		LinkedList listC = preload_ini_file(ini_list, "C");
		LinkedList listH = preload_ini_file(ini_list, "H");
		LinkedList listN = preload_ini_file(ini_list, "N");
		LinkedList listS = preload_ini_file(ini_list, "S");
		LinkedList listO = preload_ini_file(ini_list, "O");
		LinkedList listP = preload_ini_file(ini_list, "P");
		// possibly the following line might be causing trouble
		
		boolean tolerance_type = true;
		double ppm_val = 50;
		
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
	}*/
	
	/*public static LinkedList create_ini_list(double C12, double N14, double H1) {
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
	}*/
	
	/**
	 * Subtract two bignum
	 * subtract bn2 from bn1
	 */
	public static BIGNUM subtract(BIGNUM bn1, BIGNUM bn2) {
		BIGNUM newbn = new BIGNUM(-bn2.NUM, bn2.E);
		return add(bn1, newbn);
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
				//mass = (mass / charge + charge * 1.007825);
				mass = (mass / charge + 1.007825);
				pattern.addIsotope(new IsotopeContainer(mass, intensity));
				//System.out.println(mass + "\t" + intensity);
			}
			count++;
		}
		return pattern;
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
	/*public static int retrieve_num_element(String str, String query_element) {
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
	}*/
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
}
