package ISOTOPEDISTRIBUTION.MS1;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import MISC.ToolBox;

public class GenerateStatisticalMatrix {

	public static void main(String[] args) {
		
	}
	public static void execute(String[] args) {
		computeStatistics(args[0], args[2], args[3], new Integer(args[1]));
	}
	
	public static void computeStatistics(String inputFile, String outputFile, String outputFile2, int resolution) {
		
		try {
			
			int totalFormulas = 0;
			int total = 0;
			HashMap map_0 = new HashMap();
			HashMap map_1 = new HashMap();
			HashMap map_2 = new HashMap();
			HashMap map_3 = new HashMap();
			HashMap map_4 = new HashMap();
			HashMap map_5 = new HashMap();
			
			HashMap map0_distribution = new HashMap();
			HashMap map1_distribution = new HashMap();
			HashMap map2_distribution = new HashMap();
			HashMap map3_distribution = new HashMap();
			HashMap map4_distribution = new HashMap();
			HashMap map5_distribution = new HashMap();
			
			int[] histogram = new int[6500];
			for (int i = 0; i < histogram.length; i++) {
				histogram[i] = 0;
			}
			
			//int buffer = 1000;
			int buffer = resolution;
			double[] histogram_score = new double[6500 * buffer];
			for (int i = 0; i < histogram_score.length; i++) {
				histogram_score[i] = 0;
			}
			//String outputFile = "C:\\Users\\tshaw\\Desktop\\PROTEOMICS\\Software\\Decharging\\Summarized_Human_Overlap_001.txt"; //args[0];
			//String outputFile = "C:\\Users\\tshaw\\Desktop\\PROTEOMICS\\Software\\Decharging\\Summarized_Human_Overlap_01.txt"; //args[0];
			FileWriter fwriter = new FileWriter(outputFile);
			BufferedWriter out = new BufferedWriter(fwriter);
			
			////args[0];
			//String outputFile2 = "C:\\Users\\tshaw\\Desktop\\PROTEOMICS\\Software\\Decharging\\Score_Similarity_01.txt"; //args[0];
			//String outputFile2 = "C:\\Users\\tshaw\\Desktop\\PROTEOMICS\\Software\\Decharging\\Score_Similarity_001.txt"; //args[0];
			FileWriter fwriter2 = new FileWriter(outputFile2);
			BufferedWriter out2 = new BufferedWriter(fwriter2);
			
			//String fileName = "C:\\Users\\tshaw\\Desktop\\PROTEOMICS\\Software\\Decharging\\Summarized_Human_001.result"; //args[0];
			//String fileName = "C:\\Users\\tshaw\\Desktop\\PROTEOMICS\\Software\\Decharging\\Summarized_Human.result"; //args[0];
			String fileName = inputFile;
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				double val = new Double(split[0]);
				val = new Double(new Double(val * buffer).intValue()) / buffer;
				int count = new Integer(split[3]);
				String MS1_I = split[4];
				String MS2_I = split[5];
				String MS3_I = split[6];
				String MS4_I = split[7];
				String MS5_I = split[8];
				String MS6_I = split[9];
				
				String MS1_M = split[16];
				String MS2_M = split[17];
				String MS3_M = split[18];
				String MS4_M = split[19];
				String MS5_M = split[20];
				String MS6_M = split[21];
				
				LinkedList spectrum = new LinkedList();
				if (!MS1_I.equals("NaN")) {
					spectrum.add(MS1_M + "\t" + MS1_I);
				}
				if (!MS2_I.equals("NaN")) {
					spectrum.add(MS2_M + "\t" + MS2_I);
				}
				if (!MS3_I.equals("NaN")) {
					spectrum.add(MS3_M + "\t" + MS3_I);
				}
				if (!MS4_I.equals("NaN")) {
					spectrum.add(MS4_M + "\t" + MS4_I);
				}
				if (!MS5_I.equals("NaN")) {
					spectrum.add(MS5_M + "\t" + MS5_I);
				}
				if (!MS6_I.equals("NaN")) {
					spectrum.add(MS6_M + "\t" + MS6_I);
				}
				
				totalFormulas += count;
				total++;
				int M = new Integer(split[2]);
				if (M == 0) {
					map_0.put(val, count);
					map0_distribution.put(val,  spectrum);
				} else if (M == 1) {
					map_1.put(val, count);
					map1_distribution.put(val,  spectrum);
				} else if (M == 2) {
					map_2.put(val, count);
					map2_distribution.put(val,  spectrum);
				} else if (M == 3) {
					map_3.put(val, count);
					map3_distribution.put(val,  spectrum);
				} else if (M == 4) {
					map_4.put(val, count);
					map4_distribution.put(val,  spectrum);
				} else if (M == 5) {
					map_5.put(val, count);
					map5_distribution.put(val,  spectrum);
				}
			}
			in.close();
			
			int overlap_0_1 = 0;
			int overlap_1_2 = 0;
			int overlap_2_3 = 0;
			int overlap_3_4 = 0;
			int overlap_4_5 = 0;
			
			int count01 = 0;
			int count12 = 0;
			int count23 = 0;
			int count34 = 0;
			int count45 = 0;
			Iterator itr = map_0.keySet().iterator();
			while (itr.hasNext()) {
				double val = (Double)itr.next();
				
				if (map_1.containsKey(val)) {
					int count0 = (Integer)map_0.get(val);
					int count1 = (Integer)map_1.get(val);
					overlap_0_1++;
					
					count01 += count0 + count1;
					int index = new Double(val).intValue();
					histogram[index] += (count0 + count1);
					
					index = new Double(val * buffer).intValue();
					LinkedList isotope0 = (LinkedList)map0_distribution.get(val);
					LinkedList isotope1 = (LinkedList)map1_distribution.get(val);
					
					histogram_score[index] = ToolBox.comparePeaks(isotope0, isotope1); 
					//out.write(val + "\t" + (count0 + count1) + "\n");
				}
			}
			itr = map_1.keySet().iterator();
			while (itr.hasNext()) {
				double val = (Double)itr.next();
				if (map_2.containsKey(val)) {
					int count1 = (Integer)map_1.get(val);
					int count2 = (Integer)map_2.get(val);
					overlap_1_2++;
					count12 += count1 + count2;
					int index = new Double(val).intValue();
					histogram[index] += (count1 + count2);
					//out.write(val + "\t" + (count1 + count2) + "\n");
					
					index = new Double(val * buffer).intValue();
					LinkedList isotope1 = (LinkedList)map1_distribution.get(val);
					LinkedList isotope2 = (LinkedList)map2_distribution.get(val);
					
					histogram_score[index] = ToolBox.comparePeaks(isotope1, isotope2);
				}
			}
			itr = map_2.keySet().iterator();
			while (itr.hasNext()) {
				double val = (Double)itr.next();
				if (map_3.containsKey(val)) {
					int count2 = (Integer)map_2.get(val);
					int count3 = (Integer)map_3.get(val);
					overlap_2_3++;
					count23 += count2 + count3;
					int index = new Double(val).intValue();
					histogram[index] += (count2 + count3);
					//out.write(val + "\t" + (count2 + count3) + "\n");
					
					index = new Double(val * buffer).intValue();
					LinkedList isotope2 = (LinkedList)map2_distribution.get(val);
					LinkedList isotope3 = (LinkedList)map3_distribution.get(val);
					
					histogram_score[index] = ToolBox.comparePeaks(isotope2, isotope3);
				}
			}
			itr = map_3.keySet().iterator();
			while (itr.hasNext()) {
				double val = (Double)itr.next();
				if (map_4.containsKey(val)) {
					int count3 = (Integer)map_3.get(val);
					int count4 = (Integer)map_4.get(val);
					overlap_3_4++;
					count34 += count3 + count4;
					int index = new Double(val).intValue();
					histogram[index] += (count3 + count4);
					//out.write(val + "\t" + (count3 + count4) + "\n");
					
					index = new Double(val * buffer).intValue();
					LinkedList isotope3 = (LinkedList)map3_distribution.get(val);
					LinkedList isotope4 = (LinkedList)map4_distribution.get(val);
					
					histogram_score[index] = ToolBox.comparePeaks(isotope3, isotope4);
				}
			}
			itr = map_4.keySet().iterator();
			while (itr.hasNext()) {
				double val = (Double)itr.next();
				if (map_5.containsKey(val)) {
					int count4 = (Integer)map_4.get(val);
					int count5 = (Integer)map_5.get(val);
					overlap_4_5++;
					count45 += count4 + count5;
					int index = new Double(val).intValue();
					histogram[index] += (count4 + count5);
					//out.write(val + "\t" + (count4 + count5) + "\n");
					
					index = new Double(val * buffer).intValue();
					LinkedList isotope4 = (LinkedList)map4_distribution.get(val);
					LinkedList isotope5 = (LinkedList)map5_distribution.get(val);
					
					histogram_score[index] = ToolBox.comparePeaks(isotope4, isotope5);
				}
			}
			
			for (int i = 0; i < histogram.length; i++) {
				out.write(i + "\t" + histogram[i] + "\n");
			}
			for (int i = 0; i < histogram_score.length; i++) {
				if (histogram_score[i] > 0) {
					out2.write(i + "\t" + histogram_score[i] + "\n");
				}
			}

			System.out.println("Total number of entries: " + total);			
			System.out.println("Number of entries sharing same mass index between M and M+1: " + overlap_0_1);			
			System.out.println("Number of entries sharing same mass index between M+1 and M+2: " + overlap_1_2);			
			System.out.println("Number of entries sharing same mass index between M+2 and M+3: " + overlap_2_3);			
			System.out.println("Number of entries sharing same mass index between M+3 and M+4: " + overlap_3_4);			
			System.out.println("Number of entries sharing same mass index between M+4 and M+5: " + overlap_4_5);						
			System.out.println("Total number of peptides: " + totalFormulas);
			System.out.println("Number peptides sharing same mass index between M and M+1 formulas: " + count01);
			System.out.println("Number peptides sharing same mass index between M+1 and M+2 formulas: " + count12);
			System.out.println("Number peptides sharing same mass index between M+2 and M+3 formulas: " + count23);
			System.out.println("Number peptides sharing same mass index between M+3 and M+4 formulas: " + count34);
			System.out.println("Number peptides sharing same mass index between M+4 and M+5 formulas: " + count45);
			
			out.close();
			out2.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
