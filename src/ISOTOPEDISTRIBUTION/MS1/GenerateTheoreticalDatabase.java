package ISOTOPEDISTRIBUTION.MS1;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;


/**
 * Parse through the degraded peptide isotope to generate 
 * @author tshaw
 *
 */
public class GenerateTheoreticalDatabase {

	public static void main(String[] args) {
		try {
			
		
			THEORETICAL_DATABASE[] database_M0 = new THEORETICAL_DATABASE[650000];
			THEORETICAL_DATABASE[] database_M1 = new THEORETICAL_DATABASE[650000];
			THEORETICAL_DATABASE[] database_M2 = new THEORETICAL_DATABASE[650000];
			THEORETICAL_DATABASE[] database_M3 = new THEORETICAL_DATABASE[650000];
			THEORETICAL_DATABASE[] database_M4 = new THEORETICAL_DATABASE[650000];
			THEORETICAL_DATABASE[] database_M5 = new THEORETICAL_DATABASE[650000];
			
			for (int i = 0; i < 650000; i++) {
				database_M0[i] = new THEORETICAL_DATABASE();
				database_M1[i] = new THEORETICAL_DATABASE();
				database_M2[i] = new THEORETICAL_DATABASE();
				database_M3[i] = new THEORETICAL_DATABASE();
				database_M4[i] = new THEORETICAL_DATABASE();
				database_M5[i] = new THEORETICAL_DATABASE();
			}
			double mass = -1;
			int count = 0;
			String fileName = args[0];
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String peptide = split[7];
				String isotope = split[15];
				double highest_intensity = getHIntensity(isotope);
				double highest_mass = getHImass(isotope);
				int DBindex = new Double(highest_mass * 100).intValue();
				int position = getHMindex(isotope);
				if (mass < highest_mass) {
					mass = highest_mass;
				}
				count++;
				
				database_M0[DBindex].MASS.add(highest_mass);
				database_M0[DBindex].PEAKPOSITION = position;
				double m0_mass = getMMass(isotope, 0);
				double m1_mass = getMMass(isotope, 1);
				double m2_mass = getMMass(isotope, 2);
				double m3_mass = getMMass(isotope, 3);
				double m4_mass = getMMass(isotope, 4);
				double m5_mass = getMMass(isotope, 5);
				double m0_intensity = getMIntensity(isotope, 0);
				double m1_intensity = getMIntensity(isotope, 1);
				double m2_intensity = getMIntensity(isotope, 2);
				double m3_intensity = getMIntensity(isotope, 3);
				double m4_intensity = getMIntensity(isotope, 4);
				double m5_intensity = getMIntensity(isotope, 5);
				if (m0_mass != Double.NaN) {
					database_M0[DBindex].M0_I.add(m0_intensity);
					database_M0[DBindex].M0_M.add(m0_mass);
				}
				if (m1_mass != Double.NaN) {
					database_M0[DBindex].M1_I.add(m1_intensity);
					database_M0[DBindex].M1_M.add(m1_mass);
				}
				if (m2_mass != Double.NaN) {
					database_M0[DBindex].M2_I.add(m2_intensity);
					database_M0[DBindex].M2_M.add(m2_mass);
				}
				if (m3_mass != Double.NaN) {
					database_M0[DBindex].M3_I.add(m3_intensity);
					database_M0[DBindex].M3_M.add(m3_mass);
				}
				if (m4_mass != Double.NaN) {
					database_M0[DBindex].M4_I.add(m4_intensity);
					database_M0[DBindex].M4_M.add(m4_mass);
				}
				if (m5_mass != Double.NaN) {
					database_M0[DBindex].M5_I.add(m5_intensity);
					database_M0[DBindex].M5_M.add(m5_mass);
				}
			}
			in.close();
			updateMeanSD(database_M0);
			//System.out.println(count + "\t" + mass);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void execute(String[] args) {
		try {
			int totalSize = 6500000;
			
			THEORETICAL_DATABASE[] database_M0 = new THEORETICAL_DATABASE[totalSize];
			THEORETICAL_DATABASE[] database_M1 = new THEORETICAL_DATABASE[totalSize];
			THEORETICAL_DATABASE[] database_M2 = new THEORETICAL_DATABASE[totalSize];
			THEORETICAL_DATABASE[] database_M3 = new THEORETICAL_DATABASE[totalSize];
			THEORETICAL_DATABASE[] database_M4 = new THEORETICAL_DATABASE[totalSize];
			THEORETICAL_DATABASE[] database_M5 = new THEORETICAL_DATABASE[totalSize];
			
			for (int i = 0; i < totalSize; i++) {
				database_M0[i] = new THEORETICAL_DATABASE();
				database_M1[i] = new THEORETICAL_DATABASE();
				database_M2[i] = new THEORETICAL_DATABASE();
				database_M3[i] = new THEORETICAL_DATABASE();
				database_M4[i] = new THEORETICAL_DATABASE();
				database_M5[i] = new THEORETICAL_DATABASE();
			}
			double mass = -1;
			int count = 0;
			String fileName = args[0];
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String peptide = split[7];
				String isotope = split[15];
				double highest_intensity = getHIntensity(isotope);
				double highest_mass = getHImass(isotope);
				int DBindex = new Double(highest_mass * 1000).intValue();
				int position = getHMindex(isotope);
				if (mass < highest_mass) {
					mass = highest_mass;
				}
				count++;
				//System.out.println(highest_mass + "\t" + position);
				if (position == 0) {
					updateDatabase(database_M0, DBindex, highest_mass, position, isotope);
				} else if (position == 1) {
					updateDatabase(database_M1, DBindex, highest_mass, position, isotope);
				} else if (position == 2) {
					updateDatabase(database_M2, DBindex, highest_mass, position, isotope);
				} else if (position == 3) {
					updateDatabase(database_M3, DBindex, highest_mass, position, isotope);
				} else if (position == 4) {
					updateDatabase(database_M4, DBindex, highest_mass, position, isotope);
				} else if (position == 5) {
					updateDatabase(database_M5, DBindex, highest_mass, position, isotope);
				}
				//(THEORETICAL_DATABASE[] database_charge1, int DBindex, double highest_mass, int position, String isotope) {
			}
			in.close();
			updateMeanSD(database_M0);
			updateMeanSD(database_M1);
			updateMeanSD(database_M2);
			updateMeanSD(database_M3);
			updateMeanSD(database_M4);
			updateMeanSD(database_M5);
			//System.out.println(cont + "\t" + mass);
			
			String outputFile = args[1];
			FileWriter fwriter = new FileWriter(outputFile);
			BufferedWriter out = new BufferedWriter(fwriter);
			
			for (int i = 0; i < totalSize; i++) {
				if (database_M0[i].PEAKPOSITION >= 0){ 
					out.write(database_M0[i].FINAL_MASS + "\t1\t" + database_M0[i].PEAKPOSITION  
							+ "\t" + database_M0[i].TOTAL + "\t" + database_M0[i].M0_I_avg
							+ "\t" + database_M0[i].M1_I_avg + "\t" + database_M0[i].M2_I_avg + "\t" + database_M0[i].M3_I_avg
							+ "\t" + database_M0[i].M4_I_avg + "\t" + database_M0[i].M5_I_avg + "\t" + database_M0[i].M0_I_stdev
							+ "\t" + database_M0[i].M1_I_stdev + "\t" + database_M0[i].M2_I_stdev + "\t" + database_M0[i].M3_I_stdev
							+ "\t" + database_M0[i].M4_I_stdev + "\t" + database_M0[i].M5_I_stdev + "\t" + database_M0[i].M0_M_avg
							+ "\t" + database_M0[i].M1_M_avg + "\t" + database_M0[i].M2_M_avg + "\t" + database_M0[i].M3_M_avg
							+ "\t" + database_M0[i].M4_M_avg + "\t" + database_M0[i].M5_M_avg + "\t" + database_M0[i].M0_M_stdev
							+ "\t" + database_M0[i].M1_M_stdev + "\t" + database_M0[i].M2_M_stdev + "\t" + database_M0[i].M3_M_stdev
							+ "\t" + database_M0[i].M4_M_stdev + "\t" + database_M0[i].M5_M_stdev + "\n");					
				}
				if (database_M1[i].PEAKPOSITION >= 0){ 
					out.write(database_M1[i].FINAL_MASS + "\t1\t" + database_M1[i].PEAKPOSITION  
							+ "\t" + database_M1[i].TOTAL + "\t" + database_M1[i].M0_I_avg
							+ "\t" + database_M1[i].M1_I_avg + "\t" + database_M1[i].M2_I_avg + "\t" + database_M1[i].M3_I_avg
							+ "\t" + database_M1[i].M4_I_avg + "\t" + database_M1[i].M5_I_avg + "\t" + database_M1[i].M0_I_stdev
							+ "\t" + database_M1[i].M1_I_stdev + "\t" + database_M1[i].M2_I_stdev + "\t" + database_M1[i].M3_I_stdev
							+ "\t" + database_M1[i].M4_I_stdev + "\t" + database_M1[i].M5_I_stdev + "\t" + database_M1[i].M0_M_avg
							+ "\t" + database_M1[i].M1_M_avg + "\t" + database_M1[i].M2_M_avg + "\t" + database_M1[i].M3_M_avg
							+ "\t" + database_M1[i].M4_M_avg + "\t" + database_M1[i].M5_M_avg + "\t" + database_M1[i].M0_M_stdev
							+ "\t" + database_M1[i].M1_M_stdev + "\t" + database_M1[i].M2_M_stdev + "\t" + database_M1[i].M3_M_stdev
							+ "\t" + database_M1[i].M4_M_stdev + "\t" + database_M1[i].M5_M_stdev + "\n");					
				}
				if (database_M2[i].PEAKPOSITION >= 0){ 
					out.write(database_M2[i].FINAL_MASS + "\t1\t" + database_M2[i].PEAKPOSITION  
							+ "\t" + database_M2[i].TOTAL + "\t" + database_M2[i].M0_I_avg
							+ "\t" + database_M2[i].M1_I_avg + "\t" + database_M2[i].M2_I_avg + "\t" + database_M2[i].M3_I_avg
							+ "\t" + database_M2[i].M4_I_avg + "\t" + database_M2[i].M5_I_avg + "\t" + database_M2[i].M0_I_stdev
							+ "\t" + database_M2[i].M1_I_stdev + "\t" + database_M2[i].M2_I_stdev + "\t" + database_M2[i].M3_I_stdev
							+ "\t" + database_M2[i].M4_I_stdev + "\t" + database_M2[i].M5_I_stdev + "\t" + database_M2[i].M0_M_avg
							+ "\t" + database_M2[i].M1_M_avg + "\t" + database_M2[i].M2_M_avg + "\t" + database_M2[i].M3_M_avg
							+ "\t" + database_M2[i].M4_M_avg + "\t" + database_M2[i].M5_M_avg + "\t" + database_M2[i].M0_M_stdev
							+ "\t" + database_M2[i].M1_M_stdev + "\t" + database_M2[i].M2_M_stdev + "\t" + database_M2[i].M3_M_stdev
							+ "\t" + database_M2[i].M4_M_stdev + "\t" + database_M2[i].M5_M_stdev + "\n");					
				}
				if (database_M3[i].PEAKPOSITION >= 0){ 
					out.write(database_M3[i].FINAL_MASS + "\t1\t" + database_M3[i].PEAKPOSITION  
							+ "\t" + database_M3[i].TOTAL + "\t" + database_M3[i].M0_I_avg
							+ "\t" + database_M3[i].M1_I_avg + "\t" + database_M3[i].M2_I_avg + "\t" + database_M3[i].M3_I_avg
							+ "\t" + database_M3[i].M4_I_avg + "\t" + database_M3[i].M5_I_avg + "\t" + database_M3[i].M0_I_stdev
							+ "\t" + database_M3[i].M1_I_stdev + "\t" + database_M3[i].M2_I_stdev + "\t" + database_M3[i].M3_I_stdev
							+ "\t" + database_M3[i].M4_I_stdev + "\t" + database_M3[i].M5_I_stdev + "\t" + database_M3[i].M0_M_avg
							+ "\t" + database_M3[i].M1_M_avg + "\t" + database_M3[i].M2_M_avg + "\t" + database_M3[i].M3_M_avg
							+ "\t" + database_M3[i].M4_M_avg + "\t" + database_M3[i].M5_M_avg + "\t" + database_M3[i].M0_M_stdev
							+ "\t" + database_M3[i].M1_M_stdev + "\t" + database_M3[i].M2_M_stdev + "\t" + database_M3[i].M3_M_stdev
							+ "\t" + database_M3[i].M4_M_stdev + "\t" + database_M3[i].M5_M_stdev + "\n");					
				}
				if (database_M4[i].PEAKPOSITION >= 0){ 
					out.write(database_M4[i].FINAL_MASS + "\t1\t" + database_M4[i].PEAKPOSITION  
							+ "\t" + database_M4[i].TOTAL + "\t" + database_M4[i].M0_I_avg
							+ "\t" + database_M4[i].M1_I_avg + "\t" + database_M4[i].M2_I_avg + "\t" + database_M4[i].M3_I_avg
							+ "\t" + database_M4[i].M4_I_avg + "\t" + database_M4[i].M5_I_avg + "\t" + database_M4[i].M0_I_stdev
							+ "\t" + database_M4[i].M1_I_stdev + "\t" + database_M4[i].M2_I_stdev + "\t" + database_M4[i].M3_I_stdev
							+ "\t" + database_M4[i].M4_I_stdev + "\t" + database_M4[i].M5_I_stdev + "\t" + database_M4[i].M0_M_avg
							+ "\t" + database_M4[i].M1_M_avg + "\t" + database_M4[i].M2_M_avg + "\t" + database_M4[i].M3_M_avg
							+ "\t" + database_M4[i].M4_M_avg + "\t" + database_M4[i].M5_M_avg + "\t" + database_M4[i].M0_M_stdev
							+ "\t" + database_M4[i].M1_M_stdev + "\t" + database_M4[i].M2_M_stdev + "\t" + database_M4[i].M3_M_stdev
							+ "\t" + database_M4[i].M4_M_stdev + "\t" + database_M4[i].M5_M_stdev + "\n");					
				}
				if (database_M5[i].PEAKPOSITION >= 0){ 
					out.write(database_M5[i].FINAL_MASS + "\t1\t" + database_M5[i].PEAKPOSITION  
							+ "\t" + database_M5[i].TOTAL + "\t" + database_M5[i].M0_I_avg
							+ "\t" + database_M5[i].M1_I_avg + "\t" + database_M5[i].M2_I_avg + "\t" + database_M5[i].M3_I_avg
							+ "\t" + database_M5[i].M4_I_avg + "\t" + database_M5[i].M5_I_avg + "\t" + database_M5[i].M0_I_stdev
							+ "\t" + database_M5[i].M1_I_stdev + "\t" + database_M5[i].M2_I_stdev + "\t" + database_M5[i].M3_I_stdev
							+ "\t" + database_M5[i].M4_I_stdev + "\t" + database_M5[i].M5_I_stdev + "\t" + database_M5[i].M0_M_avg
							+ "\t" + database_M5[i].M1_M_avg + "\t" + database_M5[i].M2_M_avg + "\t" + database_M5[i].M3_M_avg
							+ "\t" + database_M5[i].M4_M_avg + "\t" + database_M5[i].M5_M_avg + "\t" + database_M5[i].M0_M_stdev
							+ "\t" + database_M5[i].M1_M_stdev + "\t" + database_M5[i].M2_M_stdev + "\t" + database_M5[i].M3_M_stdev
							+ "\t" + database_M5[i].M4_M_stdev + "\t" + database_M5[i].M5_M_stdev + "\n");					
				}

			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void updateDatabase(THEORETICAL_DATABASE[] database_charge1, int DBindex, double highest_mass, int position, String isotope) {
		database_charge1[DBindex].TOTAL++;
		database_charge1[DBindex].MASS.add(highest_mass);
		database_charge1[DBindex].PEAKPOSITION = position;
		double m0_mass = getMMass(isotope, 0);
		double m1_mass = getMMass(isotope, 1);
		double m2_mass = getMMass(isotope, 2);
		double m3_mass = getMMass(isotope, 3);
		double m4_mass = getMMass(isotope, 4);
		double m5_mass = getMMass(isotope, 5);
		double m0_intensity = getMIntensity(isotope, 0);
		double m1_intensity = getMIntensity(isotope, 1);
		double m2_intensity = getMIntensity(isotope, 2);
		double m3_intensity = getMIntensity(isotope, 3);
		double m4_intensity = getMIntensity(isotope, 4);
		double m5_intensity = getMIntensity(isotope, 5);
		if (m0_mass != Double.NaN) {
			database_charge1[DBindex].M0_I.add(m0_intensity);
			database_charge1[DBindex].M0_M.add(m0_mass);
		}
		if (m1_mass != Double.NaN) {
			database_charge1[DBindex].M1_I.add(m1_intensity);
			database_charge1[DBindex].M1_M.add(m1_mass);
		}
		if (m2_mass != Double.NaN) {
			database_charge1[DBindex].M2_I.add(m2_intensity);
			database_charge1[DBindex].M2_M.add(m2_mass);
		}
		if (m3_mass != Double.NaN) {
			database_charge1[DBindex].M3_I.add(m3_intensity);
			database_charge1[DBindex].M3_M.add(m3_mass);
		}
		if (m4_mass != Double.NaN) {
			database_charge1[DBindex].M4_I.add(m4_intensity);
			database_charge1[DBindex].M4_M.add(m4_mass);
		}
		if (m5_mass != Double.NaN) {
			database_charge1[DBindex].M5_I.add(m5_intensity);
			database_charge1[DBindex].M5_M.add(m5_mass);
		}
	}
	public static void updateMeanSD(THEORETICAL_DATABASE[] database_charge) {
		for (int i = 0; i < database_charge.length; i++) {
			StandardDeviation sd = new StandardDeviation();
			database_charge[i].M0_I_stdev = sd.evaluate(list2double(database_charge[i].M0_I));
			database_charge[i].M0_I_avg = mean(database_charge[i].M0_I);
			database_charge[i].M1_I_stdev = sd.evaluate(list2double(database_charge[i].M1_I));
			database_charge[i].M1_I_avg = mean(database_charge[i].M1_I);
			database_charge[i].M2_I_stdev = sd.evaluate(list2double(database_charge[i].M2_I));
			database_charge[i].M2_I_avg = mean(database_charge[i].M2_I);
			database_charge[i].M3_I_stdev = sd.evaluate(list2double(database_charge[i].M3_I));
			database_charge[i].M3_I_avg = mean(database_charge[i].M3_I);
			database_charge[i].M4_I_stdev = sd.evaluate(list2double(database_charge[i].M4_I));
			database_charge[i].M4_I_avg = mean(database_charge[i].M4_I);
			database_charge[i].M5_I_stdev = sd.evaluate(list2double(database_charge[i].M5_I));
			database_charge[i].M5_I_avg = mean(database_charge[i].M5_I);
			
			database_charge[i].M0_M_stdev = sd.evaluate(list2double(database_charge[i].M0_M));
			database_charge[i].M0_M_avg = mean(database_charge[i].M0_M);
			database_charge[i].M1_M_stdev = sd.evaluate(list2double(database_charge[i].M1_M));
			database_charge[i].M1_M_avg = mean(database_charge[i].M1_M);
			database_charge[i].M2_M_stdev = sd.evaluate(list2double(database_charge[i].M2_M));
			database_charge[i].M2_M_avg = mean(database_charge[i].M2_M);
			database_charge[i].M3_M_stdev = sd.evaluate(list2double(database_charge[i].M3_M));
			database_charge[i].M3_M_avg = mean(database_charge[i].M3_M);
			database_charge[i].M4_M_stdev = sd.evaluate(list2double(database_charge[i].M4_M));
			database_charge[i].M4_M_avg = mean(database_charge[i].M4_M);
			database_charge[i].M5_M_stdev = sd.evaluate(list2double(database_charge[i].M5_M));
			database_charge[i].M5_M_avg = mean(database_charge[i].M5_M);
			
			database_charge[i].FINAL_MASS = mean(database_charge[i].MASS);
		}
		
		
	}
	public static double mean(LinkedList list) {
		double total = 0;
		Iterator itr = list.iterator();
		while (itr.hasNext()) {
			total += (Double)itr.next();			
		}
		return total / list.size();
	}
	
	public static double[] list2double(LinkedList list) {
		double[] result = new double[list.size()];
		int n = 0;
		Iterator itr = list.iterator();
		while (itr.hasNext()) {
			double val = (Double)itr.next();
			result[n] = val;
			n++;
		}
		return result;
	}
	public static double getMMass(String isotope, int j) {
		String[] split = isotope.split(",");
		int index = -1;
		double return_mass = -1;
		double highestIntensity = -1;
		for (int i = 1; i < split.length; i++) {
			if (i == j + 1) {
				double mass = new Double(split[i].split(":")[0]);
				double intensity = new Double(split[i].split(":")[1]);
				return mass;				
			}			
		}
		return Double.NaN;
	}
	public static double getMIntensity(String isotope, int j) {
		String[] split = isotope.split(",");
		int index = -1;
		double return_mass = -1;
		double highestIntensity = -1;
		for (int i = 1; i < split.length; i++) {
			if (i == j + 1) {
				double mass = new Double(split[i].split(":")[0]);
				double intensity = new Double(split[i].split(":")[1]);
				return intensity;
				
			}
			
		}
		return Double.NaN;
	}
	
	public static class THEORETICAL_DATABASE {
		LinkedList MASS = new LinkedList();
		int PEAKPOSITION = -1;
		int TOTAL = 0;
		LinkedList M0_M = new LinkedList();
		LinkedList M1_M = new LinkedList();
		LinkedList M2_M = new LinkedList();
		LinkedList M3_M = new LinkedList();
		LinkedList M4_M = new LinkedList();
		LinkedList M5_M = new LinkedList();
		
		double FINAL_MASS = -1;
		
		double M0_I_avg = -1;
		double M1_I_avg = -1;
		double M2_I_avg = -1;
		double M3_I_avg = -1;
		double M4_I_avg = -1;
		double M5_I_avg = -1;
		
		double M0_I_stdev= -1;
		double M1_I_stdev = -1;
		double M2_I_stdev = -1;
		double M3_I_stdev = -1;
		double M4_I_stdev = -1;
		double M5_I_stdev = -1;
		
		double M0_M_avg = -1;
		double M1_M_avg = -1;
		double M2_M_avg = -1;
		double M3_M_avg = -1;
		double M4_M_avg = -1;
		double M5_M_avg = -1;
		
		double M0_M_stdev= -1;
		double M1_M_stdev = -1;
		double M2_M_stdev = -1;
		double M3_M_stdev = -1;
		double M4_M_stdev = -1;
		double M5_M_stdev = -1;
		
		LinkedList M0_I = new LinkedList();
		LinkedList M1_I = new LinkedList();
		LinkedList M2_I = new LinkedList();
		LinkedList M3_I = new LinkedList();
		LinkedList M4_I = new LinkedList();
		LinkedList M5_I = new LinkedList();
	}
	public static double getHIntensity(String isotope) {
		String[] split = isotope.split(",");
		int index = -1;
		double return_mass = -1;
		double highestIntensity = -1;
		for (int i = 1; i < split.length; i++) {
			double mass = new Double(split[i].split(":")[0]);
			double intensity = new Double(split[i].split(":")[1]);
			if (highestIntensity < intensity) {
				highestIntensity = intensity;
				index = i - 1;
				
			}
		}
		return highestIntensity;
	}
	public static double getHImass(String isotope) {
		String[] split = isotope.split(",");
		int index = -1;
		double return_mass = -1;
		double highestIntensity = -1;
		for (int i = 1; i < split.length; i++) {
			double mass = new Double(split[i].split(":")[0]);
			double intensity = new Double(split[i].split(":")[1]);
			if (highestIntensity < intensity) {
				highestIntensity = intensity;
				index = i - 1;
				return_mass = mass;
			}
		}
		return return_mass;
	}
	public static int getHMindex(String isotope) {
		String[] split = isotope.split(",");
		int index = -1;
		double highestIntensity = -1;
		for (int i = 1; i < split.length; i++) {
			double mass = new Double(split[i].split(":")[0]);
			double intensity = new Double(split[i].split(":")[1]);
			if (highestIntensity < intensity) {
				highestIntensity = intensity;
				index = i - 1;
			}
		}
		return index;
	}
}
